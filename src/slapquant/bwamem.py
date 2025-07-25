import pathlib
import logging
import os
import re
import subprocess
from queue import Queue
from typing import Literal
import multiprocessing
import threading
import tempfile

from ._utils import CandidateAlignment

logger = logging.getLogger('bwamem')

BWA_PATH = os.environ.get('BWA_PATH', 'bwa-mem2')

SEQUENCE_RE = re.compile("^[A-Za-z]+$")

ASSETS_PATH = pathlib.Path(__file__).parent / "__assets__/"


class FASTAValidationError(Exception):
    def __init__(self, file: pathlib.Path):
        super().__init__(f"{file} is not a valid FASTA file.")


class BWAMEM:
    def __init__(self, fasta_file: pathlib.Path):
        # Uses a temporary directory to store the FASTA index
        self._working_directory = tempfile.TemporaryDirectory(dir=os.getcwd())
        logger.debug(
            'Creating temporary directory to hold bwa-mem index '
            f'"{self._working_directory.name}".'
        )
        self.fasta_file = pathlib.Path(fasta_file)
        # Validate the FASTA file
        self._validate_fasta()
        # Create the index.
        self._index_fasta()

    def _validate_fasta(self):
        def sequences():
            sequence: str = ""
            header: str | None = None
            with self.fasta_file.open("r") as fasta:
                for line in fasta:
                    if line[0] == ">":
                        if header is not None:
                            yield header, sequence
                        header = line[1:].strip()
                        sequence = ""
                    elif header is None:
                        # First line isn't a header
                        raise FASTAValidationError(self.fasta_file)
                    else:
                        sequence += line.strip()

        logger.info("Validating FASTA file")
        hashes = set()
        for header, sequence in sequences():
            headerhash = hash(header)
            if (
                headerhash in hashes or                 # Duplicated header
                SEQUENCE_RE.match(sequence) is None     # Invalid sequence
            ):
                raise FASTAValidationError(self.fasta_file)

    def _index_fasta(self):
        logger.info('Indexing reference genome.')
        tempdir = pathlib.Path(self._working_directory.name)
        self.reference_fasta = tempdir / 'reference.fasta'
        logger.debug('Trying to hard link reference fasta.')
        os.link(self.fasta_file, self.reference_fasta)

        # Run BWA-MEM index in a subprocess.
        logger.debug('Starting indexing...')
        try:
            subprocess.run([BWA_PATH, "index", self.reference_fasta],
                           check=True, capture_output=True)
        except subprocess.CalledProcessError as error:
            logger.error(
                "BWA indexing failed. Output of BWA index was:\n"
                "---STDOUT---\n"
                f"{error.stdout}"
                "---STDERR---\n"
                f"{error.stderr}"
            )
            raise
        logger.debug('Indexing finished.')

    def align(
        self,
        readfiles: list[pathlib.Path | str],
        queue: Queue[CandidateAlignment],
        slsequence: str | None = None,
        polyA_sequence: str | None = None,
        threads=None,
        this_logger: logging.Logger = logger,
        filtering: Literal['softclipped'] | Literal['allmatches'] = (
            'softclipped'
        ),
    ):
        if threads is None:
            threads = multiprocessing.cpu_count()-2

        if (
            filtering == 'allmatches' and
            slsequence is None and
            polyA_sequence is None
        ):
            raise ValueError(
                "When aligning with 'allmatches', the `slsequence` and"
                "`polyA_sequence` arguments must be given."
            )

        # This starts BWA-MEM in a subprocess. The SAM output is piped into a
        # gawk program which filters the alignments and outputs every matching
        # softclipped alignment in the tab separated format
        # sequence_name position strand  match_location clipped_sequence
        # We read this output line-by-line and put these candidate alignments
        # into a queue for further processing. This is done in a separate
        # thread to allow for parallel processing.
        def line_generator():
            for fname in readfiles:
                this_logger.info('Starting to align reads.')
                cmd = [BWA_PATH, "mem", "-t",
                       str(threads), "-Y", self.reference_fasta, fname]
                this_logger.debug(
                    f'Starting BWA-MEM: "{" ".join([str(v) for v in cmd])}')
                bwa = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    text=True,
                    stderr=subprocess.PIPE
                )
                awk_cmd = [
                    'gawk', "-f",
                    ASSETS_PATH / f"{filtering}.awk"
                ]
                if filtering == 'allmatches':
                    awk_cmd.extend(
                        [
                            '-v', f'slsequence={slsequence}',
                            '-v', f'polyAsequence={polyA_sequence}'
                        ]
                    )
                this_logger.debug(
                    f'Starting gawk: "{" ".join([str(v) for v in awk_cmd])}')

                awk_env = os.environ.copy()
                awk_env["AWKPATH"] = str(ASSETS_PATH)
                awk = subprocess.Popen(
                    awk_cmd,
                    text=True,
                    stdin=bwa.stdout,
                    stdout=subprocess.PIPE,
                    env=awk_env,
                )

                this_logger.debug(
                    'Processing and queuing incoming aligned reads.')
                for line in iter(awk.stdout.readline, ''):
                    queue.put(CandidateAlignment.from_line(line))
                awk.stdout.close()
                awk_return_code = awk.wait()
                if awk_return_code:
                    raise subprocess.CalledProcessError(awk_return_code, cmd)
                bwa_return_code = bwa.wait()
                if bwa_return_code:
                    stdout, stderr = bwa.communicate()
                    logger.error(
                        "BWA mem failed. Output of BWA mem was:\n"
                        "---STDOUT---\n"
                        f"{stdout}"
                        "---STDERR---\n"
                        f"{stderr}"
                    )
                    raise subprocess.CalledProcessError(bwa_return_code, cmd)

            this_logger.info('Aligning done.')
            queue.put(None)

        # Start the thread, make it a daemon (so that it is killed if the main
        # thread is terminated)
        bwa_thread = threading.Thread(target=line_generator)
        bwa_thread.daemon = True
        bwa_thread.start()

        # Return the thread so it can be joined when it is done.
        return bwa_thread
