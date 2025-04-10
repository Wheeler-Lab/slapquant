import logging
import pathlib
import argparse
import logging
logging.basicConfig(filename='/dev/stderr', format='%(levelname)s:\t%(message)s')
logger = logging.getLogger()

from slapquant.slapquant import process_reads
from slapquant.assign import assign_sites, identify_UTRs

from geffa.geffa import Seq

def slapquant_main():
    parser = argparse.ArgumentParser(description="Find spliced leader acceptor and polyadenylation sites by aligning high-quality RNASeq reads to an existing genome. Note that the reads net to be trimmed beforehand.")
    parser.add_argument('reference_genome', type=pathlib.Path, help="""The path to a FASTA file containing the reference genome used to align the RNASeq reads to.""")
    parser.add_argument('rnaseq_reads', nargs='+', type=pathlib.Path, help="""The path(s) to (potentially multiple) FASTQ files containing the RNASeq reads. Note that no special handling for paired-end reads is done.""")
    parser.add_argument('-S', '--spliced-leader-sequence', help="""The spliced leader sequence to look for. If not specified, slapquant will attempt to auto-detect the spliced leader sequence. Your mileage might vary.""", default=None)
    parser.add_argument('-v', '--verbose', help="""Give more info about the process.""", action="store_const", dest="loglevel", const=logging.INFO, default=logging.WARNING)
    parser.add_argument('-d', '--debug', help="""Debugging info (very verbose)""", action="store_const", dest="loglevel", const=logging.DEBUG)

    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    sl_sequence = args.spliced_leader_sequence
    if sl_sequence is not None:
        sl_sequence = Seq(sl_sequence)
    gff = process_reads(args.reference_genome, args.rnaseq_reads, sl_sequence)
    gff.save('/dev/stdout')

def slapidentify_main():
    parser = argparse.ArgumentParser(description="Identify spliced leader sequence from high-quality RNAseq reads. Note that the reads need to be trimmed beforehand.")
    parser.add_argument('reference_genome', type=pathlib.Path, help="""The path to a FASTA file containing the reference genome used to align the RNASeq reads to.""")
    parser.add_argument('rnaseq_reads', nargs='+', type=pathlib.Path, help="""The path(s) to (potentially multiple) FASTQ files containing the RNASeq reads. Note that no special handling for paired-end reads is done.""")
    parser.add_argument('-v', '--verbose', help="""Give more info about the process.""", action="store_const", dest="loglevel", const=logging.INFO, default=logging.WARNING)
    parser.add_argument('-d', '--debug', help="""Debugging info (very verbose)""", action="store_const", dest="loglevel", const=logging.DEBUG)

    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    sl_sequence = process_reads(args.reference_genome, args.rnaseq_reads, None)
    print(sl.sequence, file=open('/dev/stdout', 'w'))

def slapassign_main():
    parser = argparse.ArgumentParser(description="Assign spliced leader accaptor and polyadenylation sites to gene models. The resulting GFF file contains the sites assigned to the given gene models and is printed to standar output.")
    parser.add_argument('gene_models_gff', type=pathlib.Path, help="The path to a GFF file containing the gene models the sites should be matched to. Please note that this must correspond to the original genome slapquant was run on to identify the sites.")
    parser.add_argument('slas_pas_gff', type=pathlib.Path, help="The path to a GFF file containing the spliced leader accaptor and polyadenylation sites identified from a prior slapquant run.")
    parser.add_argument('-v', '--verbose', help="""Give more info about the process.""", action="store_const", dest="loglevel", const=logging.INFO, default=logging.WARNING)

    args = parser.parse_args()
    logger.setLevel(args.loglevel)
    
    gff = assign_sites(args.gene_models_gff, args.slas_pas_gff)
    gff.save('/dev/stdout')

def slaputrs_main():
    parser = argparse.ArgumentParser(description="Call untranslated regions (UTRs) of trans-spliced genes based on spliced leader acceptor sites and polyadenylation sites assigned to genes. This takes the GFF file created by slapassign and identifies UTRs. The resulting GFF file contains has the 3' and 5' UTRs added to all genes for which they could be identified.")
    parser.add_argument('gene_models_slas_pas_gff', nargs='?', type=pathlib.Path, default='/dev/stdin', help="The path to a GFF file containing the gene models and the spliced leader accaptor and polyadenylation sites identified and assigned via a prior slapassign run.")
    parser.add_argument('-v', '--verbose', help="""Give more info about the process.""", action="store_const", dest="loglevel", const=logging.INFO, default=logging.WARNING)

    args = parser.parse_args()
    logger.setLevel(args.loglevel)
    
    gff = identify_UTRs(args.gene_models_slas_pas_gff)
    gff.save('/dev/stdout')
