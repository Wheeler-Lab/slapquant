from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor
import itertools
import pathlib
from queue import Queue
import re
import threading
from typing import Counter, Literal, NamedTuple
import os
import tempfile
import subprocess
import multiprocessing

from ._utils import QueueConsumer, QueueTee

from geffa.geffa import Seq, SLASNode, PASNode
from geffa import GffFile

import logging
logger = logging.getLogger('slapquant')

BWA_PATH = os.environ.get('BWA_PATH', 'bwa-mem2')

class CandidateAlignment(NamedTuple):
    sequence_name: str
    position: int
    strand: Literal['+'] | Literal['-']
    match_location: Literal['start'] | Literal['end']
    sequence: str

class BWAMEM:
    def __init__(self, fasta_file: pathlib.Path):
        # Uses a temporary directory to store the FASTA index
        self._working_directory = tempfile.TemporaryDirectory(dir=os.getcwd())
        logger.debug(f'Creating temporary directory to hold bwa-mem index "{self._working_directory.name}".')
        self.fasta_file = pathlib.Path(fasta_file)
        # Create the index.
        self._index_fasta()

    def _index_fasta(self):
        logger.info('Indexing reference genome.')
        tempdir = pathlib.Path(self._working_directory.name)
        self.reference_fasta = tempdir / 'reference.fasta'
        logger.debug('Trying to hard link reference fasta.')
        os.link(self.fasta_file, self.reference_fasta)

        # Run BWA-MEM index in a subprocess.
        logger.debug('Starting indexing...')
        subprocess.run([BWA_PATH, "index", self.reference_fasta], check=True, capture_output=True)
        logger.debug('Indexing finished.')

    def align(self, readfiles: list[pathlib.Path | str], queue: Queue[CandidateAlignment], threads=None, this_logger: logging.Logger=logger):
        if threads is None:
            threads = multiprocessing.cpu_count()-2

        # This starts BWA-MEM in a subprocess. The SAM output is piped into a gawk program which filters
        # the alignments and outputs every matching softclipped alignment in the tab separated format
        # sequence_name position    strand  match_location  clipped_sequence
        # We read this output line-by-line and put these candidate alignments into a queue for further
        # processing. This is done in a separate thread to allow for parallel processing.
        def line_generator():
            for fname in readfiles:
                this_logger.info(f'Starting to align reads.')
                cmd =  [BWA_PATH, "mem", "-t", str(threads), "-Y", self.reference_fasta, fname]
                this_logger.debug(f'Starting BWA-MEM: "{" ".join([str(v) for v in cmd])}')
                bwa = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE)
                awk_cmd = ['gawk', "-f", pathlib.Path(__file__).parent / "__assets__/filtering.awk"]
                this_logger.debug(f'Starting gawk: "{" ".join([str(v) for v in awk_cmd])}')
                awk = subprocess.Popen(awk_cmd, text=True, stdin=bwa.stdout, stdout=subprocess.PIPE)

                this_logger.debug('Processing and queuing incoming aligned reads.')
                for line in iter(awk.stdout.readline, ''):
                    sequence_name, position, strand, match_location, sequence = line.split()
                    queue.put(CandidateAlignment(
                        sequence_name,
                        int(position),
                        strand,
                        match_location,
                        sequence
                    ))
                awk.stdout.close()
                awk_return_code = awk.wait()
                if awk_return_code:
                    raise subprocess.CalledProcessError(awk_return_code, cmd)
                bwa_return_code = bwa.wait()
                if bwa_return_code:
                    raise subprocess.CalledProcessError(bwa_return_code, cmd)

            this_logger.info('Aligning done.')
            queue.put(None)

        # Start the thread, make it a daemon (so that it is killed if the main thread is terminated)
        bwa_thread = threading.Thread(target=line_generator)
        bwa_thread.daemon = True
        bwa_thread.start()

        # Return the thread so it can be joined when it is done.
        return bwa_thread

def do_alignment(bwa: BWAMEM, readfiles: list[pathlib.Path | str], bwa_threads: int | None=None, this_logger: logging.Logger=logger) -> tuple[Queue[CandidateAlignment], threading.Thread]:
    # This queue will hold the pre-filtered softclipped sequences from BWA-MEM's alignment.
    queue = Queue[CandidateAlignment](maxsize=100000)
    # Start aligning
    bwa_thread = bwa.align(readfiles, queue, threads=bwa_threads, this_logger=this_logger)
    # Return the queue and the thread (for future joining when it is done)
    return queue, bwa_thread

def find_sl_sequence(softclipped: Queue[CandidateAlignment]):
    # Take the softclipped sequences from the given queue and find the spliced leader sequence.
    # This sequence should be the most abundant softclipped sequence after any poly-A tails.

    # Simple dict that holds the count each sequence occurred with
    counts = Counter[Seq]()
    # Store how many sequences we've processed
    processed: int = [0]
    # We only want to process the first N sequences, so that the queues don't become deadlocked.
    # Use the queue maxsize as N.
    to_process = softclipped.maxsize
    # We'll raise this event once processed[0] == to_process.
    finished = threading.Event()
    def process(alignment: CandidateAlignment):
        if processed[0] < to_process:
            processed[0] += 1
            sequence = Seq(alignment.sequence)
            if alignment.match_location == 'end':
                sequence = sequence.reverse_complement()
            # We don't want sequences shorter than 8 bases
            if len(sequence) < 8 or 'TTTTTT' in sequence or 'AAAAAA' in sequence:  # Remove polyA motifs
                return
            # Count sequence
            # TODO: Could use a Bloom filter to make it more performant.
            counts[sequence] += 1
        else:
            # Now we're ignoring everything that comes down the queue. We need to do this to avoid deadlock.
            # Set the event to say we're finished
            finished.set()
    # Process the softclipped sequences
    consumer = QueueConsumer(process, softclipped)
    consumer.start()
    # Wait until the finished event is set.
    finished.wait()

    # Return the most abundant sequence we've found (that's also longer than 8 bases)
    # Also return the thread, needs to be joined later because we still need to empty the queue.
    sl_sequence = counts.most_common(1)[0][0]
    return sl_sequence, consumer

class FilterPattern:
    def __init__(self, pattern: str | re.Pattern[str], strand: Literal['+'] | Literal['-'], match_location: Literal['start'] | Literal['end']):
        # Compile the regex if it's not already done.
        if isinstance(pattern, str):
            pattern = re.compile(pattern)
        self.pattern = pattern
        self.strand = strand
        self.match_location = match_location
    
    def __call__(self, alignment: CandidateAlignment):
        # We want the pattern to match the sequence, but also the strand and the start or end of the alignment.
        return (
            (alignment.match_location == self.match_location) and
            (alignment.strand == self.strand) and
            (self.pattern.match(alignment.sequence) is not None)
        )

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

def filter_alignments_by_clipped_sequence(softclipped: Queue[CandidateAlignment], patterns: list[FilterPattern]):
    # Simple counter that holds the usage count of each site.
    positions = Counter[Site]()
    def process(alignment: CandidateAlignment):
        # For each candidate soft clipped sequence we look if it matches one of the given patterns.
        # This not only checks the pattern, but also strand and the respective end at which it was clipped.
        for pattern in patterns:
            if pattern(alignment):
                positions[Site.from_alignment(alignment)] += 1
    consumer = QueueConsumer(process, softclipped)
    consumer.start()

    return positions, consumer

def create_gff(reference_fasta: pathlib.Path, SLAS_counter: Counter[Site], PAS_counter: Counter[Site]):
    # Generate GFF to output the SLAS and PAS sites.
    gff = GffFile(fasta_file=reference_fasta)
    SLASIDcounter = itertools.count(1)
    for (sequence_name, position, strand), usage in SLAS_counter.items():
        SLASNode(-1, gff.sequence_regions[sequence_name], 'slapquant', 'SLAS', position, position+1, '.', strand, '.', f'ID={sequence_name}_SLAS_{next(SLASIDcounter)};usage={usage}')
    PASIDcounter = itertools.count(1)
    for (sequence_name, position, strand), usage in PAS_counter.items():
        PASNode(-1, gff.sequence_regions[sequence_name], 'slapquant', 'PAS', position, position+1, '.', strand, '.', f'ID={sequence_name}_PAS_{next(PASIDcounter)};usage={usage}')
    return gff
    

def process_reads(reference_genome: pathlib.Path, rnaseq_reads: list[pathlib.Path], sl_sequence: Seq | None = None):
    # Index the genome
    bwa = BWAMEM(reference_genome)
    n_cpus = multiprocessing.cpu_count()
    # Running the alignment in parallel with 8 cores allocated to each worker seems to be most efficient.
    n_workers = min(n_cpus // 8, len(rnaseq_reads))
    n_bwa_threads = n_cpus // n_workers
    with ProcessPoolExecutor(n_workers, initializer=_init_worker, initargs=(bwa, n_bwa_threads, sl_sequence)) as executor:
        sites = executor.map(_process_read_file, rnaseq_reads)
    
    # Now that we've gotten the sites back from the worker processes, we need to collect them into their respective counter objects.
    spliced_leader_sites = Counter[Site]()
    polyA_sites = Counter[Site]()
    for spliced, polyA in sites:
        spliced_leader_sites.update(spliced)
        polyA_sites.update(polyA)

    # Return a GFF object that holds the found sites
    return create_gff(reference_genome, spliced_leader_sites, polyA_sites)

def _init_worker(_bwa: BWAMEM, _n_threads: int, _sl_sequence: Seq | None):
    # This just initialises the parallel worker process with the BWAMEM object and the number of threads the alignment should use.
    global bwa, n_threads, sl_sequence
    bwa = _bwa
    n_threads = _n_threads
    sl_sequence = _sl_sequence

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

    # Take the alignments queue and copy anything that gets dropped into it into two or three separate queues,
    # depending on if we need to find the spliced leader sequence.
    # If we need to find it, one of the queues will take the first N alignments and looks for the appropriate spliced leader sequence.
    # The queue feeding the spliced reads identification step is paused until the spliced leader sequence is identified,
    # so the queue maxsize needs to be at least N, otherwise we run into a deadlock situation. Therefore,
    # `find_sl_sequence` uses the maxsize of the queue it is given as N.
    alignments_tee = QueueTee(alignments, copies=3 if sl_sequence is None else 2, maxsize=100000)

    # Start the identification of reads containing the polyadenylation motif
    # The returned `polyA_alignments` will contain a Counter (a special dict that contains the usage of each site)
    # once the `polyA_alignments_thread` is joined.
    polyA_sites, polyA_sites_thread = filter_alignments_by_clipped_sequence(
        alignments_tee.outputs[0],
        [
            # We look for at least 6 As softclipped at the end of a read if we're on the forward strand...
            FilterPattern('^A{6,}', '+', 'end'),
            # ...and for at least 6 Ts softclipped at the start of a read if we're on the reverse strand.
            FilterPattern('T{6,}$', '-', 'start'),
        ]
    )
    threads.append(polyA_sites_thread)

    if sl_sequence is None:
        # Find the spliced leader sequence.
        this_logger.info(f'Starting spliced leader sequence finding.')
        sl_sequence, spliced_leader_thread = find_sl_sequence(alignments_tee.outputs[2])
        # Only use the 8 ending nucleotides of the spliced leader sequence.
        # This is enough to accurately identify spliced reads, and leads to less reads being unecessarily discarded.
        sl_sequence = sl_sequence[-8:]
        this_logger.info(f'Found spliced leader sequence {sl_sequence}.')
        threads.append(spliced_leader_thread)
    else:
        this_logger.info(f"Spliced leader sequence specified ('{sl_sequence}'), skipping automatic detection")

    # Now identify softclipped reads that contain the spliced leader sequence.
    spliced_leader_sites, spliced_leader_sites_thread = filter_alignments_by_clipped_sequence(
        alignments_tee.outputs[1],
        [
            # We look for reads that have a softclipped segment at the start containing the spliced leader sequence if we're on the forward strand...
            FilterPattern(f'{sl_sequence}$', '+', 'start'),
            # ... and for reads that have a softclipped segment at the end containing the reverse complement of the spliced leader sequence if we're on the reverse strand.
            FilterPattern(f'^{sl_sequence.reverse_complement()}', '-', 'end'),
        ]
    )
    threads.append(spliced_leader_sites_thread)

    # Crucially, we cannot join the spliced_leader_thread before, because it needs to flush items out of the queue until the end.
    for thread in threads[::-1]:
        thread.join()
    
    this_logger.info(f"Found {len(spliced_leader_sites)} unique spliced leader acceptor sites and {len(polyA_sites)} unique polyadenylation sites.")

    return spliced_leader_sites, polyA_sites
