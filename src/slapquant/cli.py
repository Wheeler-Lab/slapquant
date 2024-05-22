import logging
import pathlib
import argparse

from slapquant.slapquant import process_reads

def main():
    parser = argparse.ArgumentParser(description="Find spliced leader acceptor and polyadenylation sites by aligning high-quality RNASeq reads to an existing genome. Note that the reads net to be trimmed beforehand.")
    parser.add_argument('reference_genome', type=pathlib.Path, help="""The path to a FASTA file containing the reference genome used to align the RNASeq reads to.""")
    parser.add_argument('rnaseq_reads', nargs='+', type=pathlib.Path, help="""The path(s) to (potentially multiple) FASTQ files containing the RNASeq reads. Note that no special handling for paired-end reads is done.""")
    parser.add_argument('-v', '--verbose', help="""Give more info about the process.""", action="store_const", dest="loglevel", const=logging.INFO, default=logging.WARNING)
    parser.add_argument('-d', '--debug', help="""Debugging info (very verbose)""", action="store_const", dest="loglevel", const=logging.DEBUG)

    args = parser.parse_args()
    logging.basicConfig(filename='/dev/stderr', level=args.loglevel)

    gff = process_reads(args.reference_genome, args.rnaseq_reads)
    gff.save('/dev/stdout')