import pathlib
import logging
from queue import Queue
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict, Counter
from typing import NamedTuple
from bisect import bisect_left, bisect_right
import threading

from ._utils import QueueTee, CandidateAlignment, QueueConsumer
from .bwamem import BWAMEM

from tqdm.auto import tqdm
from geffa.geffa import Seq, GffFile, GeneNode
import pandas as pd

logger = logging.getLogger('slapspan')

def do_alignment(bwa: BWAMEM, readfiles: list[pathlib.Path | str], bwa_threads: int | None=None, this_logger: logging.Logger=logger) -> tuple[Queue[CandidateAlignment], threading.Thread]:
    # This queue will hold the pre-filtered softclipped sequences from BWA-MEM's alignment.
    queue = Queue[CandidateAlignment](maxsize=100000)
    # Start aligning
    bwa_thread = bwa.align(readfiles, queue, slsequence=sl_sequence, threads=bwa_threads, this_logger=this_logger, filtering='allmatches')
    # Return the queue and the thread (for future joining when it is done)
    return queue, bwa_thread

class Site(NamedTuple):
    position: int
    usage: int
    gene: str

class SiteCollection:
    def __init__(self, sites: list[Site]):
        self._sites = sorted(sites, key=lambda s: s.position)

    def __getitem__(self, positionrange: slice | int):
        start = bisect_left(self._sites, positionrange.start, key=lambda s: s.position)
        end = bisect_right(self._sites, positionrange.stop, key=lambda s: s.position)
        return SiteCollection(self._sites[start:end])

    def total_usage(self):
        totals = Counter()
        for s in self._sites:
            totals[s.gene] += s.usage
        return dict(totals)
    
    def __len__(self):
        return len(self._sites)

def process_reads(slas_pas_file: pathlib.Path, reference_genome: pathlib.Path, rnaseq_reads: list[pathlib.Path], sl_sequence: Seq | None = None):
    # Open GFF file containing SLAS and PAS features (from a previous slapquant run)
    gff = GffFile(slas_pas_file, ignore_unknown_feature_types=True)
    SLAS = {seqreg.name: SiteCollection([Site(node.start, int(node.attributes['usage']), node.parents[0].attributes['ID'] if node.parents else None) for node in seqreg.node_registry.values() if node.type == 'SLAS']) for seqreg in gff.sequence_regions.values()}
    PAS = {seqreg.name: SiteCollection([Site(node.start, int(node.attributes['usage']), node.parents[0].attributes['ID'] if node.parents else None) for node in seqreg.node_registry.values() if node.type == 'PAS']) for seqreg in gff.sequence_regions.values()}

    # Index the genome
    bwa = BWAMEM(reference_genome)
    n_cpus = multiprocessing.cpu_count()
    # Running the alignment in parallel with 8 cores allocated to each worker seems to be most efficient.
    n_workers = min(n_cpus // 8, len(rnaseq_reads))
    n_bwa_threads = n_cpus // n_workers
    with ProcessPoolExecutor(n_workers, initializer=_init_worker, initargs=(bwa, n_bwa_threads, sl_sequence, SLAS, PAS)) as executor:
        sites_iterator = executor.map(_process_read_file, rnaseq_reads, chunksize=1)
        if not logger.isEnabledFor(logging.INFO):
            spans = list(tqdm(sites_iterator, total=len(rnaseq_reads), desc="Read files"))
        else:
            spans = list(sites_iterator)

    slas_usage = Counter[tuple[str, str]]()
    pas_usage = Counter[tuple[str, str]]()
    for span in spans:
        slas_usage.update(span.slas_usage)
        pas_usage.update(span.pas_usage)

    slas_spans = pd.Series(slas_usage).rename_axis(["contig", "gene"]).rename("weighted_slas_spans")
    pas_spans = pd.Series(pas_usage).rename_axis(["contig", "gene"]).rename("weighted_pas_spans")

    return pd.concat([slas_spans, pas_spans], axis=1).fillna(0).astype(int)

def _init_worker(_bwa: BWAMEM, _n_threads: int, _sl_sequence: str, _SLAS: dict[str, SiteCollection], _PAS: dict[str, SiteCollection]):
    # This just initialises the parallel worker process with the BWAMEM object and the number of threads the alignment should use.
    global bwa, n_threads, sl_sequence, SLAS, PAS
    bwa = _bwa
    n_threads = _n_threads
    sl_sequence = _sl_sequence
    SLAS = _SLAS
    PAS = _PAS

class ReadFileResults(NamedTuple):
    read_file: pathlib.Path
    slas_usage: Counter[tuple[str, str]]
    pas_usage: Counter[tuple[str, str]]

def count_spans(queue: Queue[CandidateAlignment]):
    global SLAS, PAS
    slas_usage = Counter[tuple[str, str]]()
    pas_usage = Counter[tuple[str, str]]()
    def process(alignment: CandidateAlignment):
        region = alignment.sequence_name
        if region in SLAS:
            start = alignment.position
            end = start + alignment.nr_matched + 1
            for gene, usage in SLAS[alignment.sequence_name][start:end].total_usage().items():
                slas_usage[(region, gene)] += usage
            for gene, usage in PAS[alignment.sequence_name][start:end].total_usage().items():
                pas_usage[(region, gene)] += usage
    
    consumer = QueueConsumer(process, queue)
    consumer.start()

    return (slas_usage, pas_usage), consumer

def _process_read_file(reads: pathlib.Path):
    global bwa, n_threads, sl_sequence
    # This is the main worker function for the parallel processing of alignments.

    # Start thread that calls BWA-MEM and the filtering pipe.
    # Successfully aligned and filtered reads will appear in the `alignments` queue.
    # We need to take care to join the `alignment_thread` when done.
    this_logger = logger.getChild(reads.stem.replace('.fastq', ''))
    this_logger.debug(f'Starting aligner with {n_threads} threads.')
    alignments, alignment_thread = do_alignment(bwa, [reads], bwa_threads=n_threads, this_logger=this_logger)

    threads: list[threading.Thread] = [alignment_thread]

    (slas_usage, pas_usage), count_thread = count_spans(alignments)

    for thread in threads[::-1]:
        thread.join()
    
    this_logger.info(f"Found {len(slas_usage)} spliced leader acceptor sites and {len(pas_usage)} polyadenylation sites spanned by nascent transcripts.")

    return ReadFileResults(reads, slas_usage, pas_usage)

