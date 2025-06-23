import pathlib
import argparse

from slapquant.slapquant import (
    process_reads as slapquant_process_reads,
    process_reads_slapidentify as slapidentify_process_reads,
)
from slapquant.assign import assign_sites, identify_UTRs
from slapquant.slapspan import process_reads as slapspan_process_reads

from geffa.geffa import Seq

import logging
logging.basicConfig(filename='/dev/stderr',
                    format='%(levelname)s:\t%(message)s')
logger = logging.getLogger()


# Predefined spliced leader sequences for some species.
# Listed without trailing "G" due to the way slapquant works.
SL_SEQUENCES = {
    "trypanosomabrucei": "AACGCTATTATTAGAACAGTTTCTGTACTATATT",
    "leishmaniamexicana": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATT",
}
SPECIES_HELP = ", ".join([f'"{species}"' for species in SL_SEQUENCES])


def _decide_sl_sequence(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
    sl_mandatory: bool,
):
    sl_species: str = args.species
    sl_sequence: str | None = args.spliced_leader_sequence

    if sl_mandatory and not (sl_species or sl_sequence):
        parser.error("Either the SL sequence or the species argument must be given.")
    elif sl_species and sl_sequence:
        parser.error(
            "Please provide either a species or a spliced leader sequence, "
            "not both."
        )
    elif sl_species:
        try:
            sl_sequence = SL_SEQUENCES[sl_species]
        except KeyError:
            parser.error(
                f"Unknown species {sl_species} for spliced leader sequence.")
    if sl_sequence is not None:
        return Seq(sl_sequence)

def slapquant_main():
    parser = argparse.ArgumentParser(
        description=(
            "Find spliced leader acceptor and polyadenylation sites by "
            "aligning high-quality RNASeq reads to an existing genome. "
            "Note that the reads net to be trimmed beforehand."
        )
    )
    parser.add_argument(
        'reference_genome',
        type=pathlib.Path,
        help=(
            "The path to a FASTA file containing the reference genome used to "
            "align the RNASeq reads to."
        ),
    )
    parser.add_argument(
        'rnaseq_reads',
        nargs='+',
        type=pathlib.Path,
        help=(
            "The path(s) to (potentially multiple) FASTQ files containing the "
            "RNASeq reads. Note that no special handling for paired-end reads "
            "is done."
        ),
    )
    parser.add_argument(
        '-s',
        '--species',
        help=(
            "Use the predetermined spliced leader sequence for this species. "
            "If neither this nor the spliced leader sequence is given, "
            "slapquant will attempt to auto-detect the spliced leader "
            "sequence. This is not recommended.\n"
            "Spliced leader sequences are available for the following "
            " species:\n"
            f"{SPECIES_HELP}"
        ),
        default=None,
    )
    parser.add_argument(
        '-S',
        '--spliced-leader-sequence',
        help=(
            "The spliced leader sequence to look for. If neither this nor the "
            "species are given, slapquant will attempt to auto-detect the "
            "spliced leader sequence. This is not recommended."
        ),
        default=None
    )
    parser.add_argument(
        "--sl-length",
        help=(
            "Minimum length of the spliced leader sequence to look for. "
            "(default 9bp)"
        ),
        type=int,
        default=9
    )
    parser.add_argument(
        "--pa-length",
        help=(
            "Minimum length of the polyadenylation sequence to look for. "
            "(default 6bp)"
        ),
        type=int,
        default=6
    )
    parser.add_argument(
        '-v',
        '--verbose',
        help="Give more info about the process.",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
        default=logging.WARNING
    )
    parser.add_argument(
        '-d',
        '--debug',
        help="Debugging info (very verbose)",
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG
    )

    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    sl_sequence = _decide_sl_sequence(args, parser, False)

    gff = slapquant_process_reads(
        args.reference_genome,
        args.rnaseq_reads,
        sl_sequence,
        args.sl_length,
        args.pa_length
    )
    gff.save('/dev/stdout')


def slapidentify_main():
    parser = argparse.ArgumentParser(
        description=(
            "Identify spliced leader sequence from high-quality RNAseq reads. "
            "Note that the reads need to be trimmed beforehand."
        ),
    )
    parser.add_argument(
        'reference_genome',
        type=pathlib.Path,
        help=(
            "The path to a FASTA file containing the reference genome used to "
            "align the RNASeq reads to."""
        ),
    )
    parser.add_argument(
        'rnaseq_reads',
        nargs='+',
        type=pathlib.Path,
        help=(
            "The path(s) to (potentially multiple) FASTQ files containing the "
            "RNASeq reads. Note that no special handling for paired-end reads "
            "is done."
        ),
    )
    parser.add_argument(
        "--sl-length",
        help=(
            "Length of the spliced leader sequence to be identified. "
            "(default 9bp)"
        ),
        type=int,
        default=9
    )
    parser.add_argument(
        '-v',
        '--verbose',
        help="Give more info about the process.",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
        default=logging.WARNING,
    )
    parser.add_argument(
        '-d',
        '--debug',
        help="Debugging info (very verbose)",
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG,
    )

    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    sl_sequence = slapidentify_process_reads(
        args.reference_genome,
        args.rnaseq_reads,
    )
    print(sl_sequence, file=open('/dev/stdout', 'w'))


def slapassign_main():
    parser = argparse.ArgumentParser(
        description=(
            "Assign spliced leader accaptor and polyadenylation sites to gene "
            "models. The resulting GFF file contains the sites assigned to "
            "the given gene models and is printed to standar output."
        ),
    )
    parser.add_argument(
        'gene_models_gff',
        type=pathlib.Path,
        help=(
            "The path to a GFF file containing the gene models the sites "
            "should be matched to. Please note that this must correspond "
            "to the original genome slapquant was run on to identify the "
            "sites."
        ),
    )
    parser.add_argument(
        'slas_pas_gff',
        type=pathlib.Path,
        help=(
            "The path to a GFF file containing the spliced leader accaptor "
            "and polyadenylation sites identified from a prior slapquant run."
        ),
    )
    parser.add_argument(
        '--strip-existing',
        help=(
            "Strip existing SLAS / PAS sites from the gene models GFF file. "
            "Existing SLAS / PAS sites will cause an error if this option is "
            "not given."
        ),
        action="store_const",
        dest="strip_existing",
        const=True,
        default=False,
    )
    parser.add_argument(
        '--min-slas-usage', '-m',
        help=(
            "When assigning sites to genes, if SLAS sites are between a "
            "CDS and a site to be assigned, only consider those with a "
            "minimum usage count. (default 4)"
        ),
        type=int,
        default=4,
    )
    parser.add_argument(
        '--min-pas-usage', '-n',
        help=(
            "When assigning sites to genes, if PAS sites are between a "
            "CDS and a site to be assigned, only consider those with a "
            "minimum usage count. (default 2)"
        ),
        type=int,
        default=4,
    )
    parser.add_argument(
        '-v',
        '--verbose',
        help="Give more info about the process.",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
        default=logging.WARNING,
    )

    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    gff = assign_sites(
        args.gene_models_gff,
        args.slas_pas_gff,
        args.strip_existing,
        min_slas_usage=args.min_slas_usage,
        min_pas_usage=args.min_pas_usage,
    )
    gff.save('/dev/stdout')


def slaputrs_main():
    parser = argparse.ArgumentParser(
        description=(
            "Call untranslated regions (UTRs) of trans-spliced genes based on "
            "spliced leader acceptor sites and polyadenylation sites assigned "
            "to genes. This takes the GFF file created by slapassign and "
            "identifies UTRs. The resulting GFF file contains has the 3' and "
            "5' UTRs added to all genes for which they could be identified."
        ),
    )
    parser.add_argument(
        'gene_models_slas_pas_gff',
        type=pathlib.Path,
        default='/dev/stdin',
        help=(
            "The path to a GFF file containing the gene models and the "
            "spliced leader accaptor and polyadenylation sites identified "
            "and assigned via a prior slapassign run."
        ),
    )
    parser.add_argument(
        '--strip-existing',
        help=(
            "Strip existing UTRs from the gene models GFF file. "
            "Existing UTRs sites will cause an error if this option is "
            "not given."
        ),
        action="store_const",
        dest="strip_existing",
        const=True,
        default=False,
    )
    parser.add_argument(
        '--max-5UTR-length',
        help=(
            "Maximum permissible length of 5'UTRs (default 5000bp)"
        ),
        type=int,
        default=5000,
    )
    parser.add_argument(
        '--max-3UTR-length',
        help=(
            "Maximum permissible length of 3'UTRs (default 5000bp)"
        ),
        type=int,
        default=5000,
    )
    parser.add_argument(
        '-v',
        '--verbose',
        help="Give more info about the process.",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
        default=logging.WARNING,
    )

    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    gff = identify_UTRs(
        args.gene_models_slas_pas_gff,
        args.strip_existing,
        max_5utr_length=args.max_5UTR_length,
        max_3utr_length=args.max_3UTR_length,
    )
    gff.save('/dev/stdout')


def slapspan_main():
    parser = argparse.ArgumentParser(
        description=(
            "Count number of aligned RNASeq reads that span spliced "
            "leader acceptor and polyadenylation sites (i.e. nascent "
            "transcripts). Note that the reads net to be trimmed "
            "beforehand."
        ),
    )
    parser.add_argument(
        'reference_genome',
        type=pathlib.Path,
        help=(
            "The path to a FASTA file containing the reference genome used to "
            "align the RNASeq reads to."
        ),
    )
    parser.add_argument(
        'slas_pas_gff',
        type=pathlib.Path,
        help=(
            "The path to a GFF file containing the SLAS and PAS sites (they "
            "need to be assigned to genes, so probably should have been run "
            "through slapassign)."
        ),
    )
    parser.add_argument(
        'rnaseq_reads',
        nargs='+',
        type=pathlib.Path,
        help=(
            "The path(s) to (potentially multiple) FASTQ files containing the "
            "RNASeq reads. Note that no special handling for paired-end reads "
            "is done."
        ),
    )
    parser.add_argument(
        '-s',
        '--species',
        help=(
            "Use the predetermined spliced leader sequence for this species. "
            "If neither this nor the spliced leader sequence is given, "
            "slapquant will attempt to auto-detect the spliced leader "
            "sequence. This is not recommended.\n"
            "Spliced leader sequences are available for the following "
            " species:\n"
            f"{SPECIES_HELP}"
        ),
        default=None,
    )
    parser.add_argument(
        '-S',
        '--spliced-leader-sequence',
        help=(
            "The spliced leader sequence to look for. If neither this nor the "
            "species are given, slapquant will attempt to auto-detect the "
            "spliced leader sequence. This is not recommended."
        ),
        default=None
    )
    parser.add_argument(
        "--sl-length",
        help=(
            "Minimum length of the spliced leader sequence to look for. "
            "(default 9bp)"
        ),
        type=int,
        default=9
    )
    parser.add_argument(
        "--pa-length",
        help=(
            "Minimum length of the polyadenylation sequence to look for. "
            "(default 6bp)"
        ),
        type=int,
        default=6
    )
    parser.add_argument(
        '-v',
        '--verbose',
        help="Give more info about the process.",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
        default=logging.WARNING,
    )
    parser.add_argument(
        '-d',
        '--debug',
        help="Debugging info (very verbose)",
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG,
    )

    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    sl_sequence: Seq = _decide_sl_sequence(args, parser, True)[:args.sl_length]

    df = slapspan_process_reads(
        args.slas_pas_gff,
        args.reference_genome,
        args.rnaseq_reads,
        sl_sequence,
        args.pa_length,
    )
    df.to_csv('/dev/stdout')
