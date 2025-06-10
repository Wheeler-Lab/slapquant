from __future__ import annotations
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import itertools
import pathlib
from queue import Queue
import re
import threading
from typing import Counter, Literal, NamedTuple
import multiprocessing

from ._utils import QueueConsumer, QueueTee, FinishingQueue, CandidateAlignment
from .bwamem import BWAMEM

from geffa.geffa import Seq, SLASNode, PASNode
from geffa import GffFile

from tqdm.auto import tqdm

import logging
logger = logging.getLogger('slapquant')


def do_alignment(
    bwa: BWAMEM,
    readfiles: list[pathlib.Path | str],
    bwa_threads: int | None = None,
    this_logger: logging.Logger = logger
) -> tuple[Queue[CandidateAlignment], threading.Thread]:
    # This queue will hold the pre-filtered softclipped sequences from
    # BWA-MEM's alignment.
    queue = Queue[CandidateAlignment](maxsize=100000)
    # Start aligning
    bwa_thread = bwa.align(
        readfiles, queue, threads=bwa_threads, this_logger=this_logger)
    # Return the queue and the thread (for future joining when it is done)
    return queue, bwa_thread


def find_sl_sequence(
    softclipped: FinishingQueue[CandidateAlignment],
    process_mode: Literal["firstn"] | Literal["all"] = "firstn",
):
    # Take the softclipped sequences from the given queue and find the spliced
    # leader sequence. This sequence should be the most abundant softclipped
    # sequence after any poly-A tails.

    # Simple dict that holds the count each sequence occurred with
    counts = Counter[Seq]()
    # Store how many sequences we've processed
    processed: int = [0]
    if process_mode == "firstn":
        # We only want to process the first N sequences, so that the queues
        # don't become deadlocked.
        # Use the queue maxsize as N.
        to_process = softclipped.maxsize - 1
    elif process_mode == "all":
        # In this case, we want to process all alignments. Set to_process to
        # zero, the actual check is for process == "all" below.
        to_process = 0
    else:
        raise ValueError("Unknown process_mode value.")

    def process(alignment: CandidateAlignment):
        if (process_mode == "all") or (processed[0] < to_process):
            processed[0] += 1
            clipped = Seq(alignment.clipped)
            remainder = Seq(alignment.remainder)
            # We don't want sequences shorter than 9 bases
            if (
                len(clipped) < 9 or
                'TTTTTT' in clipped or
                'AAAAAA' in clipped
            ):  # Remove polyA motifs
                return
            if alignment.match_location == 'end':
                clipped = clipped.reverse_complement()
                remainder = remainder.reverse_complement()

            # The candidate SL sequence is the part that is being clipped
            # before the start of the alignment to the genome.
            # The SL sequence contains a "G" at the end. The SL acceptor site
            # is usually the dinucleotide "AG", can be others but always ends
            # with a "G". The aligner will therefore align everything up to
            # and including the "G" that should be part of the SL sequence. We
            # therefore need to recover the missing nucleotide from the
            # remainder (the aligned bit).
            # We also only use 9 nucleotides (the last nine positions) of the
            # SL sequence, this is enough to identify SLAS sites, and avoids
            # unnecessary discarding of reads due to errors in the extended
            # clipped sequence.
            sl_candidate = clipped[-8:] + remainder[0]

            # Count sequence
            # TODO: Could use a Bloom filter to make it more performant.
            counts[sl_candidate] += 1
        else:
            # Now we're ignoring everything that comes down the queue. We need
            # to do this to avoid deadlock.
            # Set the event to say we're finished
            softclipped.finished.set()
    # Process the softclipped sequences
    consumer = QueueConsumer(process, softclipped)
    consumer.start()

    if process_mode == "firstn":
        # Wait until the finished event is set.
        softclipped.finished.wait()

        # Return the most abundant sequence we've found
        # (that's also longer than 8 bases)
        # Also return the thread, needs to be joined later because we still
        # need to empty the queue.
        sl_sequence = counts.most_common(1)[0][0]
        return sl_sequence, consumer
    elif process_mode == "all":
        # Here, we don't wait, but return the full counts object to be
        # processed later.
        return counts, consumer


class FilterPattern:
    def __init__(
        self,
        pattern: str | re.Pattern[str],
        strand: Literal['+'] | Literal['-'],
        match_location: Literal['start'] | Literal['end']
    ):
        # Compile the regex if it's not already done.
        if isinstance(pattern, str):
            pattern = re.compile(pattern)
        self.pattern = pattern
        self.strand = strand
        self.match_location = match_location

        self.counter = 0

    def __del__(self):
        logger.debug(
            f"Filter pattern '{self.pattern.pattern}' strand {self.strand} "
            f"location {self.match_location} matched {self.counter} times")

    def __call__(self, alignment: CandidateAlignment):
        # The SL sequence often (always?) ends in the same nucleotide (often
        # "G") as the SL acceptor site dinucleotide. This means that this "G"
        # is parted of the aligned portion of the sequence and is not retained
        # in the softclipped part. To alleviate this, we match both the
        # softclipped part only, and separately, the softclipped part plus the
        # first nucleotide of the aligned sequence. Alternatively, we could
        # trim a trailing "G" of the given SL sequence, but we don't know if
        # this would hold for all organisms.
        if alignment.strand == '+':
            aligned_sl_nucleotide = alignment.clipped + alignment.remainder[0]
        else:
            aligned_sl_nucleotide = alignment.remainder[-1] + alignment.clipped

        # We want the pattern to match the sequence, but also the strand and
        # the start or end of the alignment.
        match = (
            (alignment.match_location == self.match_location) and
            (alignment.strand == self.strand) and (
                (self.pattern.search(alignment.clipped) is not None) or
                (self.pattern.search(aligned_sl_nucleotide) is not None)
            )
        )
        if match:
            self.counter += 1
        return match


class Site(NamedTuple):
    sequence_name: str
    position: int
    strand: Literal['+'] | Literal['-']

    @staticmethod
    def from_alignment(alignment: CandidateAlignment):
        return Site(
            alignment.sequence_name,
            alignment.position,
            alignment.strand
        )


def filter_alignments_by_clipped_sequence(
    softclipped: Queue[CandidateAlignment],
    patterns: list[FilterPattern]
):
    # Simple counter that holds the usage count of each site.
    positions = Counter[Site]()

    def process(alignment: CandidateAlignment):
        # For each candidate soft clipped sequence we look if it matches one
        # of the given patterns. This not only checks the pattern, but also
        # strand and the respective end at which it was clipped.
        for pattern in patterns:
            if pattern(alignment):
                positions[Site.from_alignment(alignment)] += 1
    consumer = QueueConsumer(process, softclipped)
    consumer.start()

    return positions, consumer


def create_gff(
    reference_fasta: pathlib.Path,
    SLAS_counter: Counter[Site],
    PAS_counter: Counter[Site]
):
    # Generate GFF to output the SLAS and PAS sites.
    gff = GffFile(fasta_file=reference_fasta)
    SLASIDcounter = itertools.count(1)
    for (sequence_name, position, strand), usage in SLAS_counter.items():
        SLASNode(
            -1,
            gff.sequence_regions[sequence_name],
            'slapquant',
            'SLAS',
            position - 1,
            position - 1,
            '.',
            strand,
            '.',
            f'ID={sequence_name}_SLAS_{next(SLASIDcounter)};usage={usage}'
        )
    PASIDcounter = itertools.count(1)
    for (sequence_name, position, strand), usage in PAS_counter.items():
        PASNode(
            -1,
            gff.sequence_regions[sequence_name],
            'slapquant',
            'PAS',
            position,
            position,
            '.',
            strand,
            '.',
            f'ID={sequence_name}_PAS_{next(PASIDcounter)};usage={usage}'
        )
    return gff


def process_reads_slapidentify(
    reference_genome: pathlib.Path,
    rnaseq_reads: list[pathlib.Path],
):
    # Index the genome
    bwa = BWAMEM(reference_genome)
    n_cpus = multiprocessing.cpu_count()
    # Running the alignment in parallel with 8 cores allocated to each worker
    # seems to be most efficient.
    n_workers = min(n_cpus // 8, len(rnaseq_reads))
    n_bwa_threads = min(n_cpus // n_workers, 8)
    results: Counter[Seq] = Counter()
    with ProcessPoolExecutor(
        n_workers,
        initializer=_init_worker,
        initargs=(bwa, n_bwa_threads, None, sl_length, None)
    ) as executor:
        sequence_iterator = executor.map(
            _process_read_file_slidentify, rnaseq_reads, chunksize=1)
        sequence_counter: Counter[Seq]
        for sequence_counter in tqdm(
            sequence_iterator,
            total=len(rnaseq_reads),
            desc="Read files",
            disable=logger.isEnabledFor(logging.INFO),
        ):
            results.update(sequence_counter)

    sl_sequence = results.most_common(1)[0][0]

    return sl_sequence


def process_reads(
    reference_genome: pathlib.Path,
    rnaseq_reads: list[pathlib.Path],
    sl_sequence: Seq | None = None,
    sl_length: int = 9,
    pa_length: int = 6,
    output_type: str = None,
):
    # Index the genome
    bwa = BWAMEM(reference_genome)
    n_cpus = multiprocessing.cpu_count()
    # Running the alignment in parallel with 8 cores allocated to each worker
    # seems to be most efficient.
    n_workers = min(n_cpus // 8, len(rnaseq_reads))
    n_bwa_threads = min(n_cpus // n_workers, 8)
    with ProcessPoolExecutor(
        n_workers,
        initializer=_init_worker,
        initargs=(bwa, n_bwa_threads, sl_sequence, sl_length, pa_length)
    ) as executor:
        sites_iterator = executor.map(
            _process_read_file, rnaseq_reads, chunksize=1)
        if not logger.isEnabledFor(logging.INFO):
            sites = list(tqdm(sites_iterator, total=len(
                rnaseq_reads), desc="Read files"))
        else:
            sites = list(sites_iterator)

    # Now that we've gotten the sites back from the worker processes, we need
    # to check if the same SL sequence was detected for all read files.
    results = defaultdict[Seq, list[ReadFileResults]](list)
    for result in sites:
        results[result.sl_sequence].append(result)
    if len(results) > 1:
        logger.warning('Found multiple different spliced leader sequences:')
        usages = {sl_sequence: sum([result.spliced_leader_sites.total(
        ) for result in entry]) for sl_sequence, entry in results.items()}
        for sl_sequence, readresults in sorted(
            results.items(), key=lambda entry: -usages[entry[0]]
        ):
            found_files = ", ".join([str(r.read_file) for r in readresults])
            usage = sum([
                result.spliced_leader_sites.total() for result in readresults
            ])
            logger.warning(
                f'\t{sl_sequence} found in {found_files}, total usage {usage}'
            )
        sl_sequence = max(usages, key=lambda s: usages[s])
        logger.warning(
            'Reported results will only include spliced leader acceptor sites'
            f'identified from the sequence {sl_sequence} with total usage'
            f'{usages[sl_sequence]}.')
        logger.warning(
            'Consider specifying the spliced leader sequence via the -S'
            'command line option.')
    else:
        sl_sequence = list(results.keys())[0]

    # return the spliced leader sequence if requested, otherwise continue to
    # GFF generation
    if output_type == "sl_sequence":
        return sl_sequence

    # Only use spliced leader acceptor sites from the SL sequence with the
    # highest total usage.
    spliced_leader_sites = Counter[Site]()
    for entry in results[sl_sequence]:
        spliced_leader_sites.update(entry.spliced_leader_sites)

    # Use polyA sites from all read files (since this isn't dependent on
    # successful detection of spliced leader sites)
    polyA_sites = Counter[Site]()
    for entry in (entry for entries in results.values() for entry in entries):
        polyA_sites.update(entry.polyA_sites)

    # Return a GFF object that holds the found sites
    return create_gff(reference_genome, spliced_leader_sites, polyA_sites)


def _init_worker(
    _bwa: BWAMEM,
    _n_threads: int,
    _sl_sequence: Seq | None,
    _sl_length: int,
    _pa_length: int
):
    # This just initialises the parallel worker process with the BWAMEM object
    # and the number of threads the alignment should use.
    global bwa, n_threads, sl_sequence, sl_length, pa_length
    bwa = _bwa
    n_threads = _n_threads
    sl_sequence = _sl_sequence
    sl_length = _sl_length
    pa_length = _pa_length


class ReadFileResults(NamedTuple):
    read_file: pathlib.Path
    sl_sequence: Seq
    spliced_leader_sites: Counter[Site]
    polyA_sites: Counter[Site]


def _process_read_file_slidentify(reads: pathlib.Path):
    global bwa, n_threads
    # This is the main worker function for the parallel processing of
    # alignments - a slimmed down version to only identify the SL
    # sequence.

    # Start thread that calls BWA-MEM and the filtering pipe.
    # Successfully aligned and filtered reads will appear in the `alignments`
    # queue. We need to take care to join the `alignment_thread` when done.
    this_logger = logger.getChild(reads.stem.replace('.fastq', ''))
    this_logger.debug(f'Starting aligner with {n_threads} threads.')
    alignments, alignment_thread = do_alignment(
        bwa, [reads], bwa_threads=n_threads, this_logger=this_logger)

    threads: list[threading.Thread] = [alignment_thread]

    # Find the spliced leader sequence.
    this_logger.info('Starting spliced leader sequence finding.')
    sl_counts, spliced_leader_thread = find_sl_sequence(
        alignments,
        process_mode="all",
    )
    threads.append(spliced_leader_thread)
    for thread in threads[::-1]:
        thread.join()

    return sl_counts


def _process_read_file(reads: pathlib.Path):
    global bwa, n_threads, sl_sequence, sl_length, pa_length
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

    # Take the alignments queue and copy anything that gets dropped into it
    # into two or three separate queues, depending on if we need to find the
    # spliced leader sequence.
    # If we need to find it, one of the queues will take the first N
    # alignments and looks for the appropriate spliced leader sequence.
    # The queue feeding the spliced reads identification step is paused until
    # the spliced leader sequence is identified, so the queue maxsize needs to
    # be at least N, otherwise we run into a deadlock situation. Therefore,
    # `find_sl_sequence` uses the maxsize of the queue it is given as N.
    alignments_tee = QueueTee(
        alignments, copies=3 if sl_sequence is None else 2, maxsize=100000)

    # Start the identification of reads containing the polyadenylation motif
    # The returned `polyA_alignments` will contain a Counter (a special dict
    # that contains the usage of each site) once the `polyA_alignments_thread`
    # is joined.
    polyA_sites, polyA_sites_thread = filter_alignments_by_clipped_sequence(
        alignments_tee.outputs[0],
        [
            # We look for at least `pa_length` As softclipped at the end of a read if
            # we're on the forward strand...
            FilterPattern(f'A{{{pa_length},}}', '+', 'end'),
            # ...and for at least `pa_length` Ts softclipped at the start of a read if
            # we're on the reverse strand.
            FilterPattern(f'T{{{pa_length},}}', '-', 'start'),
        ]
    )
    threads.append(polyA_sites_thread)

    if sl_sequence is None:
        # Find the spliced leader sequence.
        this_logger.info('Starting spliced leader sequence finding.')
        sl_sequence, spliced_leader_thread = find_sl_sequence(
            alignments_tee.outputs[2])
        this_logger.info(f'Found spliced leader sequence {sl_sequence}.')
        threads.append(spliced_leader_thread)
    else:
        this_logger.info(
            f"Spliced leader sequence specified ('{sl_sequence}'), "
            "skipping automatic detection"
        )

    # Only use the `sl_length` ending nucleotides of the spliced leader
    # sequence. This is enough to accurately identify spliced reads, and leads\
    # to less reads being unecessarily discarded.
    sl_sequence = sl_sequence[-sl_length:]

    # Now identify softclipped reads that contain the spliced leader sequence.
    spliced_leader_sites, spliced_leader_sites_thread = (
        filter_alignments_by_clipped_sequence(
            alignments_tee.outputs[1],
            [
                # We look for reads that have a softclipped segment at the
                # start containing the spliced leader sequence if we're on the
                # forward strand...
                FilterPattern(
                    f'{sl_sequence}$',
                    '+',
                    'start'
                ),
                # ... and for reads that have a softclipped segment at the end
                # containing the reverse complement of the spliced leader
                # sequence if we're on the reverse strand.
                FilterPattern(
                    f'^{sl_sequence.reverse_complement()}',
                    '-',
                    'end'
                ),
            ]
        )
    )
    threads.append(spliced_leader_sites_thread)

    # Crucially, we cannot join the spliced_leader_thread before, because it
    # needs to flush items out of the queue until the end.
    for thread in threads[::-1]:
        thread.join()

    this_logger.info(
        f"Found {len(spliced_leader_sites)} unique spliced leader acceptor "
        f" sites and {len(polyA_sites)} unique polyadenylation sites."
    )

    return ReadFileResults(
        reads,
        sl_sequence,
        spliced_leader_sites,
        polyA_sites
    )
