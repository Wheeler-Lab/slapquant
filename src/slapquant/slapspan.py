from functools import reduce
import operator
import pathlib
import logging
from queue import Queue
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
from typing import NamedTuple
from bisect import bisect_left, bisect_right
import threading

from ._utils import CandidateAlignment, QueueConsumer
from .bwamem import BWAMEM

from tqdm.auto import tqdm
from geffa.geffa import Seq, GffFile, SLASNode, PASNode, CDSNode, MRNANode, GeneNode
import pandas as pd

logger = logging.getLogger('slapspan')


def do_alignment(
    bwa: BWAMEM,
    readfiles: list[pathlib.Path | str],
    bwa_threads: int | None = None,
    this_logger: logging.Logger = logger,
) -> tuple[Queue[CandidateAlignment], threading.Thread]:
    # This queue will hold the pre-filtered softclipped sequences from
    # BWA-MEM's alignment.
    queue = Queue[CandidateAlignment](maxsize=100000)
    # Start aligning
    bwa_thread = bwa.align(
        readfiles,
        queue,
        slsequence=sl_sequence,
        polyA_sequence=polyA_sequence,
        threads=bwa_threads,
        this_logger=this_logger,
        filtering='allmatches',
    )
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
        start = bisect_left(self._sites, positionrange.start,
                            key=lambda s: s.position)
        end = bisect_right(self._sites, positionrange.stop,
                           key=lambda s: s.position)
        return SiteCollection(self._sites[start:end])

    def total_usage(self):
        totals = Counter()
        for s in self._sites:
            totals[s.gene] += s.usage
        return dict(totals)

    def __len__(self):
        return len(self._sites)


class Span(NamedTuple):
    start: int
    end: int
    gene: str


class SpanCollection:
    def __init__(self, spans: list[Span]):
        self._spans = sorted(spans, key=lambda s: s.start)

    def __getitem__(self, positionrange: slice) -> list[Span]:
        if self._spans:
            start = bisect_left(
                self._spans,
                positionrange.start,
                key=lambda s: s.start
            )
            if start > 0:
                end = bisect_right(
                    self._spans,
                    positionrange.stop,
                    key=lambda s: s.end
                )
                candidate_spans = self._spans[start-1:end+1]
                return [span for span in candidate_spans
                        if span.end >= positionrange.stop]
        return []

    def __len__(self):
        return len(self._spans)


def process_reads(
    slas_pas_file: pathlib.Path,
    reference_genome: pathlib.Path,
    rnaseq_reads: list[pathlib.Path],
    sl_sequence: Seq | None = None,
    pa_length: int = 6,
):
    # Open GFF file containing SLAS and PAS features (from a previous
    # slapquant run)
    gff = GffFile(slas_pas_file, ignore_unknown_feature_types=True)
    SLAS = {
        seqreg.name: SiteCollection(
            [
                Site(
                    node.start,
                    int(node.attributes['usage']),
                    node.parents[0].attributes['ID'],
                )
                for node in seqreg.nodes_of_type(SLASNode)
                if node.parents
            ]
        )
        for seqreg in gff.sequence_regions.values()
    }
    PAS = {
        seqreg.name: SiteCollection(
            [
                Site(
                    node.start,
                    int(node.attributes['usage']),
                    node.parents[0].attributes['ID'],
                )
                for node in seqreg.nodes_of_type(PASNode)
                if node.parents
            ]
        )
        for seqreg in gff.sequence_regions.values()
    }
    mRNA = {
        seqreg.name: SpanCollection(
            [
                Span(
                    node.start,
                    node.end,
                    node.parents[0].attributes['ID'],
                )
                for node in seqreg.nodes_of_type(MRNANode)
            ]
        )
        for seqreg in gff.sequence_regions.values()
    }

    polyA_sequence = "".join(["A"] * pa_length)

    # Index the genome
    bwa = BWAMEM(reference_genome)
    n_cpus = multiprocessing.cpu_count()
    # Running the alignment in parallel with 8 cores allocated to each worker
    # seems to be most efficient.
    n_workers = min(n_cpus // 8, len(rnaseq_reads))
    n_bwa_threads = n_cpus // n_workers
    with ProcessPoolExecutor(
        n_workers,
        initializer=_init_worker,
        initargs=(
            bwa, n_bwa_threads, sl_sequence, polyA_sequence, SLAS, PAS, mRNA),
    ) as executor:
        sites_iterator = executor.map(
            _process_read_file, rnaseq_reads, chunksize=1)
        if not logger.isEnabledFor(logging.INFO):
            results = list(tqdm(sites_iterator, total=len(
                rnaseq_reads), desc="Read files"))
        else:
            results = list(sites_iterator)

    total: ReadFileResults = reduce(operator.add, results)

    combined = pd.concat(
        [
            pd.Series(entry, name=label).rename_axis(["contig", "gene"])
            for label, entry in [
                ("SLAS_spans", total.slas_count),
                ("usage_weighted_SLAS_spans", total.slas_usage),
                ("PAS_spans", total.pas_count),
                ("usage_weighted_PAS_spans", total.pas_usage),
                ("aligned_mRNA_reads", total.mRNA_reads),
            ]
        ],
        axis=1,
    )

    combined = combined.join(
        pd.Series(
            {
                (seqreg.name, gene.attributes["ID"]): sum([
                    int(node.attributes["usage"])
                    for node in gene.children_of_type(SLASNode)
                ])
                for seqreg in gff.sequence_regions.values()
                for gene in seqreg.nodes_of_type(GeneNode)
            },
            name="total_SLAS_usage",
        ).rename_axis(["contig", "gene"]),
        how="outer",
    )
    combined = combined.join(
        pd.Series(
            {
                (seqreg.name, gene.attributes["ID"]): sum([
                    int(node.attributes["usage"])
                    for node in gene.children_of_type(PASNode)
                ])
                for seqreg in gff.sequence_regions.values()
                for gene in seqreg.nodes_of_type(GeneNode)
            },
            name="total_PAS_usage",
        ).rename_axis(["contig", "gene"]),
        how="outer",
    )

    coding_genes = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(
            {
                (seqreg.name, cds.parents[0].parents[0].attributes["ID"])
                for seqreg in gff.sequence_regions.values()
                for cds in seqreg.nodes_of_type(CDSNode)
            },
            names=["contig", "gene"]
        )
    )
    return (
        combined
        .join(coding_genes, how="outer").fillna(0).astype(int)
        .loc[:, [
            "SLAS_spans",
            "usage_weighted_SLAS_spans",
            "total_SLAS_usage",
            "PAS_spans",
            "usage_weighted_PAS_spans",
            "total_PAS_usage",
            "aligned_mRNA_reads",
        ]]
    )


def _init_worker(
    _bwa: BWAMEM,
    _n_threads: int,
    _sl_sequence: str,
    _polyA_sequence: str,
    _SLAS: dict[str, SiteCollection],
    _PAS: dict[str, SiteCollection],
    _mRNA: dict[str, SpanCollection],
):
    # This just initialises the parallel worker process with the BWAMEM object
    # and the number of threads the alignment should use.
    global bwa, n_threads, sl_sequence, polyA_sequence, SLAS, PAS, mRNA
    bwa = _bwa
    n_threads = _n_threads
    sl_sequence = _sl_sequence
    polyA_sequence = _polyA_sequence
    SLAS = _SLAS
    PAS = _PAS
    mRNA = _mRNA


class ReadFileResults:
    def __init__(self, read_file: pathlib.Path):
        self.read_files = {read_file}
        self.slas_count = Counter[tuple[str, str]]()
        self.slas_usage = Counter[tuple[str, str]]()
        self.pas_count = Counter[tuple[str, str]]()
        self.pas_usage = Counter[tuple[str, str]]()
        self.mRNA_reads = Counter[tuple[str, str]]()

    def __add__(self: "ReadFileResults", other: "ReadFileResults"):
        if self.read_files.intersection(other.read_files):
            raise ValueError("Cannot combine already combined results!")
        new = ReadFileResults(self.read_files.union(other.read_files))
        for part in [self, other]:
            part.slas_count.update(self.slas_count)
            part.slas_usage.update(self.slas_usage)
            part.pas_count.update(self.pas_count)
            part.pas_usage.update(self.pas_usage)
            part.mRNA_reads.update(self.mRNA_reads)

        return new


def count_spans(results: ReadFileResults, queue: Queue[CandidateAlignment]):
    global SLAS, PAS, mRNA

    def process(alignment: CandidateAlignment):
        region = alignment.sequence_name
        if region in SLAS:
            start = alignment.position
            end = start + alignment.nr_matched + 1
            for gene, usage in (
                SLAS[alignment.sequence_name][start:end].total_usage().items()
            ):
                results.slas_count[(region, gene)] += 1
                results.slas_usage[(region, gene)] += usage
            for gene, usage in (
                PAS[alignment.sequence_name][start:end].total_usage().items()
            ):
                results.pas_count[(region, gene)] += 1
                results.pas_usage[(region, gene)] += usage
            for mrna in mRNA[region][start:end]:
                results.mRNA_reads[(region, mrna.gene)] += 1

    consumer = QueueConsumer(process, queue)
    consumer.start()

    return consumer


def _process_read_file(reads: pathlib.Path):
    global bwa, n_threads, sl_sequence
    # This is the main worker function for the parallel processing of
    # alignments.

    # Start thread that calls BWA-MEM and the filtering pipe.
    # Successfully aligned and filtered reads will appear in the `alignments`
    # queue. We need to take care to join the `alignment_thread` when done.
    this_logger = logger.getChild(reads.stem.replace('.fastq', ''))
    this_logger.debug(f'Starting aligner with {n_threads} threads.')
    alignments, alignment_thread = do_alignment(
        bwa, [reads], bwa_threads=n_threads, this_logger=this_logger)

    threads: list[threading.Thread] = [alignment_thread]

    results = ReadFileResults(reads)

    count_thread = count_spans(results, alignments)
    threads.append(count_thread)

    for thread in threads[::-1]:
        thread.join()

    this_logger.info(
        f"Found {len(results.slas_count)} spliced leader acceptor sites and "
        f"{len(results.pas_count)} polyadenylation sites spanned by nascent "
        "transcripts."
    )

    return results
