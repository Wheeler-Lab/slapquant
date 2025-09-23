import geffa
from geffa.geffa import (
    Node,
    SLASNode,
    PASNode,
    GeneNode,
    MRNANode,
    FivePrimeUTRNode,
    ThreePrimeUTRNode
)

import logging
import __main__
import pathlib

import numpy as np

logger = logging.getLogger(pathlib.Path(__main__.__file__).stem)


def check_or_strip_nodes(
    gff: geffa.GffFile,
    types: list[str],
    strip: bool
):
    if strip:
        def node_func(node: Node) -> None:
            node.delete()
        logger.info("Stripping existing SLAS / PAS sites.")
    else:
        def node_func(node: Node) -> None:
            raise ValueError(
                "Existing sites found, cannot proceed. Please use "
                "'--strip-existing' to remove."
            )
    for node in [
        node
        for seqreg in gff.sequence_regions.values()
        for node in seqreg.node_registry.values()
    ]:
        if node.type in types:
            node_func(node)


def assign_sites(
    gene_models_gff: pathlib.Path,
    slas_pas_sites_gff: pathlib.Path,
    strip_existing: bool,
    min_slas_usage: int = 2,
    min_pas_usage: int = 2,
    min_usage_factor: float = 1.0
):
    gene_models = geffa.GffFile(
        gene_models_gff, ignore_unknown_feature_types=True)

    check_or_strip_nodes(gene_models, ["SLAS", "PAS"], strip_existing)

    slas_pas = geffa.GffFile(slas_pas_sites_gff)

    processed_regions: set[str] = set()
    for seqreg in gene_models.sequence_regions.values():
        try:
            sitesreg = slas_pas.sequence_regions[seqreg.name]
        except KeyError:
            continue
        for min_usage, node_type in [
            (min_slas_usage, SLASNode),
            (min_pas_usage, PASNode)
        ]:
            node: SLASNode | PASNode
            for node in sitesreg.nodes_of_type(node_type):
                if int(node.attributes["usage"]) >= min_usage:
                    node_type(
                        node.line_nr,
                        seqreg,
                        node.source,
                        node.type,
                        node.start,
                        node.end,
                        '.',
                        node.strand,
                        '.',
                        f"ID={node.attributes['ID']};"
                        f"usage={node.attributes['usage']}",
                    )
        processed_regions.add(seqreg.name)
    sites_only_regions = set(
        slas_pas.sequence_regions).difference(processed_regions)
    if len(sites_only_regions) > 0:
        logger.warning(
            "Gene models GFF does not contain sequence regions present in the "
            "SLAS/PAS sites GFF - possible mismatch between slapquant run and "
            "given gene models?"
        )

    for seqreg in gene_models.sequence_regions.values():
        SLAS_to_search = seqreg.nodes_of_type(SLASNode)
        PAS_to_search = seqreg.nodes_of_type(PASNode)
        CDS_to_search_SLAS = [
            mRNA.CDS_children()[0]
            for mRNA in seqreg.nodes_of_type(MRNANode)
            if len(mRNA.CDS_children()) > 0
        ]
        CDS_to_search_PAS = [
            mRNA.CDS_children()[-1]
            for mRNA in seqreg.nodes_of_type(MRNANode)
            if len(mRNA.CDS_children()) > 0
        ]
        nodes_SLAS = PAS_to_search + CDS_to_search_SLAS
        nodes_PAS = SLAS_to_search + CDS_to_search_PAS
        for slas in seqreg.nodes_of_type(SLASNode):
            nodes_to_search = sorted(
                (
                    feature for feature in nodes_SLAS
                    if (
                        (
                            (feature.type == 'PAS') and
                            (slas.strand == '+') and
                            (feature.end > slas.end) and
                            (int(feature.attributes['usage']) >
                                int(slas.attributes['usage']) * min_usage_factor)
                        ) or (
                            (feature.type == 'PAS') and
                            (slas.strand == '-') and
                            (feature.start < slas.start) and
                            (int(feature.attributes['usage']) >
                                int(slas.attributes['usage']) * min_usage_factor)
                        ) or (
                            (feature.type == 'CDS') and
                            (slas.strand == '+') and
                            (feature.strand == '+') and
                            (feature.end > slas.end)
                        ) or (
                            (feature.type == 'CDS') and
                            (slas.strand == '-') and
                            (feature.strand == '-') and
                            (feature.start < slas.start)
                        )
                    )
                ),
                key=lambda x: x.start if slas.strand == '+' else -x.end,
            )
            if len(nodes_to_search) == 0:
                logger.info(
                    f"No CDS found we could match {slas.attributes['ID']} to. "
                    "Skipping."
                )
                continue

            closest_node = nodes_to_search[0]

            if closest_node.type == 'PAS':
                logger.info(
                    f"Closest node to {slas.attributes['ID']} is a PAS, no "
                    "CDS could be assigned.")
            else:
                slas.add_parent(closest_node.parents[0].parents[0])

        for pas in seqreg.nodes_of_type(PASNode):
            nodes_to_search = sorted(
                (
                    feature for feature in nodes_PAS
                    if (
                        (
                            (feature.type == 'SLAS') and
                            (pas.strand == '+') and
                            (feature.start < pas.start) and
                            (int(feature.attributes['usage']) >
                                int(pas.attributes['usage']) * min_usage_factor)
                        ) or (
                            (feature.type == 'SLAS') and
                            (pas.strand == '-') and
                            (feature.end > pas.end) and
                            (int(feature.attributes['usage']) >
                                int(pas.attributes['usage']) * min_usage_factor)
                        ) or (
                            (feature.type == "CDS") and
                            (pas.strand == '+') and
                            (feature.strand == '+') and
                            (feature.start < pas.start)
                        ) or (
                            (feature.type == "CDS") and
                            (pas.strand == '-') and
                            (feature.strand == '-') and
                            (feature.end > pas.end)
                        )
                    )
                ),
                key=lambda x: -x.end if pas.strand == '+' else x.start
            )
            if len(nodes_to_search) == 0:
                logger.info(
                    f"No CDS found we could match {pas.attributes['ID']} to. "
                    "Skipping."
                )
                continue

            closest_node = nodes_to_search[0]

            if closest_node.type == 'SLAS':
                logger.info(
                    f"Closest node to {pas.attributes['ID']} is a SLAS, no "
                    "CDS could be assigned."
                )
            else:
                pas.add_parent(closest_node.parents[0].parents[0])

    return gene_models


def identify_UTRs(
    annotations_gff: pathlib.Path,
    strip_existing: bool,
    max_5utr_length: int = 5000,
    max_3utr_length: int = 5000,
    genome_fasta: pathlib.Path | None = None,
    fix_shorten_orfs: bool = False,
    fix_lengthen_orfs: bool = False,
):
    if fix_shorten_orfs and genome_fasta is None:
        raise ValueError(
            "If fix_shorten_orfs is True, you must provide the path of a genome FASTA file."
        )
    if fix_lengthen_orfs and genome_fasta is None:
        raise ValueError(
            "If fix_lengthen_orfs is True, you must provide the path of a genome FASTA file."
        )
    
    gff = geffa.GffFile(annotations_gff, 
        ignore_unknown_feature_types=True,
        fasta_file=genome_fasta)
    check_or_strip_nodes(
        gff, ["five_prime_UTR", "three_prime_UTR"], strip_existing)
    
    cds_shortened=0
    cds_lengthened=0 

    for seqreg in gff.sequence_regions.values():
        for gene in seqreg.nodes_of_type(GeneNode):
            mRNAs = gene.children_of_type(MRNANode)
            if len(mRNAs) == 0:
                logger.info(
                    f"{gene.attributes['ID']} is not a protein coding gene, "
                    "skipping UTR assignment."
                )
                continue
            elif len(mRNAs) > 1:
                logger.warning(
                    f"{gene.attributes['ID']} has multiple mRNAs assigned, "
                    "UTR assignment isn't implemented yet."
                )
                continue
            mRNA: MRNANode = mRNAs[0]
            CDSs = mRNA.CDS_children()
            if len(CDSs) == 0:
                logger.info(
                    f"{mRNA.attributes['ID']} has no CDS assigned, skipping UTR "
                    "assignment."
                )
                continue

            slas_sites: list[SLASNode] = []
            for slas in gene.children_of_type(SLASNode):
                    slas_sites.append(slas)

            pas_sites: list[PASNode] = []
            for pas in gene.children_of_type(PASNode):
                    pas_sites.append(pas)

            if slas_sites:
                slas = sorted(slas_sites, key=lambda x: -
                              int(x.attributes['usage']))[0]
                cds = CDSs[0]

                if mRNA.strand == '+':
                    slas_offset = (cds.start - 1) - (slas.end + 2)
                else:
                    slas_offset = (slas.start - 1) - (cds.end + 1)
                if slas_offset < 0 and fix_shorten_orfs:
                    best_start = 0
                    for i in range(3, len(cds.sequence), 3):
                        if i > slas_offset and cds.sequence[i:i+3] == "ATG":
                            best_start = i
                            break
                    if best_start > 0:
                        cds_shortened += 1
                        logger.debug(
                            f"Start codon for {mRNA.attributes['ID']} updated "
                            f"to codon number {best_start // 3}")
                        if mRNA.strand == '+':
                            cds.start = cds.start + best_start
                        else:
                            cds.end = cds.end - best_start
                    else:
                        logger.warning(
                            f"SLAS for {mRNA.attributes['ID']} is within the CDS, "
                            "and an suitable alternative start codon could not "
                            "be found.")
                
                if mRNA.strand == '+':
                    start = slas.end + 2
                    end = cds.start - 1
                else:
                    end = slas.start - 1
                    start = cds.end + 1
                utr_length = end - start + 1
                if utr_length < 1:
                    logger.warning(
                        f"Could not assign 5'UTR to {mRNA.attributes['ID']}, "
                        "the UTR would have been zero length.")
                elif utr_length > max_5utr_length:
                    logger.info(
                        f"Could not assign 5'UTR to {mRNA.attributes['ID']}, "
                        "the UTR would be longer than the allowed maximum.")
                else:
                    utr = FivePrimeUTRNode(
                        -1,
                        seqreg,
                        'slaputrs',
                        'five_prime_UTR',
                        start,
                        end,
                        '.',
                        mRNA.strand,
                        '.',
                        (
                            f'ID={mRNA.attributes["ID"]}_UTR5;'
                            f'Parent={mRNA.attributes["ID"]}'
                        ),
                    )
                    if fix_lengthen_orfs:
                        best_offset = -1
                        for i in range(len(utr.sequence), 0, -3):
                            if utr.sequence[i-3:i] in ["TAA", "TGA", "TAG"]:
                                best_offset = -1
                                break
                            if utr.sequence[i-3:i] == "ATG":
                                best_offset = i - 3
                        if best_offset > -1:
                            best_offset = len(utr.sequence) - best_offset
                            cds_lengthened += 1
                            logger.debug(
                                f"Start codon for {mRNA.attributes['ID']} updated to "
                                f"codon number {-best_offset // 3}")
                            if mRNA.strand == '+':
                                cds.start = cds.start - best_offset
                                utr.end = utr.end - best_offset
                            else:
                                cds.end = cds.end + best_offset
                                utr.start = utr.start + best_offset
            
            if pas_sites:
                positions = [
                    pas.start
                    for pas in sorted(pas_sites, key=lambda p: p.start)
                    for _ in range(int(pas.attributes["usage"]))
                ]
                if len(positions) % 2 == 0:
                    if mRNA.strand == '+':
                        del positions[0]
                    else:
                        del positions[-1]
                position = int(np.median(positions))
                cds = CDSs[-1]
                if mRNA.strand == '+':
                    start = cds.end + 1
                    end = position
                else:
                    end = cds.start - 1
                    start = position
                utr_length = end - start + 1
                if utr_length < 1:
                    logger.warning(
                        f"Could not assign 3'UTR to {mRNA.attributes['ID']}, "
                        "the UTR would have been zero length.")
                elif utr_length > max_3utr_length:
                    logger.info(
                        f"Could not assign 3'UTR to {mRNA.attributes['ID']}, "
                        "the UTR would be longer than the allowed maximum.")
                else:
                    ThreePrimeUTRNode(
                        -1,
                        seqreg,
                        'slaputrs',
                        'three_prime_UTR',
                        start,
                        end,
                        '.',
                        mRNA.strand,
                        '.',
                        (
                            f'ID={mRNA.attributes["ID"]}_UTR3;'
                            f'Parent={mRNA.attributes["ID"]}'
                        ),
                    )

    if fix_shorten_orfs and fix_lengthen_orfs:
        logger.warning(f"{cds_shortened} CDSs shortened and {cds_lengthened} CDSs lengthened") 
    elif fix_shorten_orfs:
        logger.warning("f{cds_shortened} CDSs shortened")
    elif fix_lengthen_orfs:
        logger.warning(f"{cds_lengthened} CDSs lengthened")

    return gff
