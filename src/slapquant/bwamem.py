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

from .base_aligner import BaseAligner

logger = logging.getLogger('bwamem')

BWA_PATH = os.environ.get('BWA_PATH', 'bwa mem')

class BWAMEM(BaseAligner):
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

    def _get_align_cmd(self, fname: pathlib.Path | str, threads: int, **kwargs) -> list[str]:
        return [BWA_PATH, "mem", "-t", str(threads), "-Y", self.reference_fasta, fname]
