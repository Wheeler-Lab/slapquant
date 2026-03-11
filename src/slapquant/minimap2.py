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

logger = logging.getLogger('minimap2')

MINIMAP2_PATH = os.environ.get('MINIMAP2_PATH', 'minimap2')

class Minimap2(BaseAligner):
    def _index_fasta(self):
        logger.info('Indexing reference genome with minimap2.')
        tempdir = pathlib.Path(self._working_directory.name)
        self.reference_fasta = tempdir / 'reference.fasta'
        self.reference_index = tempdir / 'reference.mmi'
        logger.debug('Trying to hard link reference fasta.')
        os.link(self.fasta_file, self.reference_fasta)

        # Run minimap2 index (create .mmi) in a subprocess.
        logger.debug('Starting indexing (minimap2 -d)...')
        try:
            subprocess.run(
                [MINIMAP2_PATH, "-d", str(self.reference_index), str(self.reference_fasta)],
                check=True,
                capture_output=True,
                text=True
            )
        except subprocess.CalledProcessError as error:
            logger.error(
                "minimap2 indexing failed. Output of minimap2 was:\n"
                "---STDOUT---\n"
                f"{error.stdout}"
                "---STDERR---\n"
                f"{error.stderr}"
            )
            raise
        logger.debug('Indexing finished.')

    def _get_align_cmd(self, fname: pathlib.Path | str, threads: int, **kwargs) -> list[str]:
        alignment_mode = kwargs.get('alignment_mode', 'map-ont')
        return [MINIMAP2_PATH, "-t", str(threads), "-a", "--secondary-seq", 
                str(self.reference_index), "-x", alignment_mode, str(fname)]
