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
from abc import ABC, abstractmethod

from ._utils import CandidateAlignment

logger = logging.getLogger('base_aligner')

SEQUENCE_RE = re.compile("^[A-Za-z]+$")

ASSETS_PATH = pathlib.Path(__file__).parent / "__assets__/"


class FASTAValidationError(Exception):
    def __init__(self, file: pathlib.Path):
        super().__init__(f"{file} is not a valid FASTA file.")


class BaseAligner(ABC):
    """
    Abstract base aligner that handles temporary directory creation,
    FASTA validation, and common threaded queue alignment processing logic.
    """

    def __init__(self, fasta_file: pathlib.Path):
        # Uses a temporary directory to store the FASTA index
        self._working_directory = tempfile.TemporaryDirectory(dir=os.getcwd())
        logger.debug(
            'Creating temporary directory to hold aligner index '
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

    @abstractmethod
    def _index_fasta(self):
        """Implement to run tool-specific indexing logic."""
        raise NotImplementedError

    @abstractmethod
    def _get_align_cmd(self, fname: pathlib.Path | str, threads: int, **kwargs) -> list[str]:
        """Implement to return tool-specific alignment command list."""
        raise NotImplementedError

    def align(
        self,
        readfiles: list[pathlib.Path | str],
        queue: Queue[CandidateAlignment],
        slsequence: str | None = None,
        polyA_sequence: str | None = None,
        mismatch_proportion: float = 0.0,
        threads: int | None = None,
        this_logger: logging.Logger = logger,
        filtering: Literal['softclipped', 'allmatches'] = 'softclipped',
        **kwargs
    ):
        if threads is None:
            threads = max(1, multiprocessing.cpu_count() - 2)

        if (
            filtering == 'allmatches' and
            slsequence is None and
            polyA_sequence is None
        ):
            raise ValueError(
                "When aligning with 'allmatches', the `slsequence` and"
                "`polyA_sequence` arguments must be given."
            )

        # This starts the aligner in a subprocess. The SAM output is piped into a
        # gawk program which filters the alignments and outputs every matching
        # softclipped alignment in the tab separated format
        # sequence_name position strand  match_location clipped_sequence
        # We read this output line-by-line and put these candidate alignments
        # into a queue for further processing. This is done in a separate
        # thread to allow for parallel processing.
        def line_generator():
            for fname in readfiles:
                this_logger.info('Starting to align reads.')
                
                cmd = self._get_align_cmd(fname=fname, threads=threads, **kwargs)
                
                this_logger.debug(
                    f'Starting aligner: "{" ".join([str(v) for v in cmd])}"'
                )
                aligner_proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    text=True,
                    stderr=subprocess.PIPE
                )
                
                awk_cmd = [
                    'gawk', "-f",
                    str(ASSETS_PATH / f"{filtering}.awk")
                ]
                if filtering == 'allmatches':
                    awk_cmd.extend(
                        [
                            '-v', f'slsequence={slsequence}',
                            '-v', f'polyAsequence={polyA_sequence}',
                            '-v', f'mismatch_prop={mismatch_proportion:.2f}'
                        ]
                    )
                this_logger.debug(
                    f'Starting gawk: "{" ".join([str(v) for v in awk_cmd])}"'
                )

                awk_env = os.environ.copy()
                awk_env["AWKPATH"] = str(ASSETS_PATH)
                awk = subprocess.Popen(
                    awk_cmd,
                    text=True,
                    stdin=aligner_proc.stdout,
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
                aligner_return_code = aligner_proc.wait()
                if aligner_return_code:
                    stdout, stderr = aligner_proc.communicate()
                    logger.error(
                        "Aligner failed. Output of aligner was:\n"
                        "---STDOUT---\n"
                        f"{stdout}"
                        "---STDERR---\n"
                        f"{stderr}"
                    )
                    raise subprocess.CalledProcessError(aligner_return_code, cmd)

            this_logger.info('Aligning done.')
            queue.put(None)

        # Start the thread, make it a daemon (so that it is killed if the main
        # thread is terminated)
        align_thread = threading.Thread(target=line_generator)
        align_thread.daemon = True
        align_thread.start()

        # Return the thread so it can be joined when it is done.
        return align_thread
