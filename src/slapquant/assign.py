import geffa

import logging
import __main__
import pathlib

logger = logging.getLogger(pathlib.Path(__main__.__file__).stem)


def check_or_strip_nodes(
    gff: geffa.GffFile,
    types: list[str],
    strip: bool
):
    if strip:
        def node_func(node: geffa.geffa.Node) -> None:
            node.delete()
        logger.info("Stripping existing SLAS / PAS sites.")
    else:
        def node_func(node: geffa.geffa.Node) -> None:
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
    min_usage: int = 4,
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
        for SLAS in sitesreg.nodes_of_type(geffa.geffa.SLASNode):
            if int(SLAS.attributes["usage"]) < min_usage:
                continue
            geffa.geffa.SLASNode(
                SLAS.line_nr,
                seqreg,
                SLAS.source,
                'SLAS',
                SLAS.start,
                SLAS.end,
                '.',
                SLAS.strand,
                '.',
                f"ID={SLAS.attributes['ID']};usage={SLAS.attributes['usage']}",
            )
        for PAS in sitesreg.nodes_of_type(geffa.geffa.PASNode):
            if int(SLAS.attributes["usage"]) < min_usage:
                continue
            geffa.geffa.PASNode(
                PAS.line_nr,
                seqreg,
                PAS.source,
                'PAS',
                PAS.start,
                PAS.end,
                '.',
                PAS.strand,
                '.',
                f"ID={PAS.attributes['ID']};usage={PAS.attributes['usage']}"
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
        SLAS_to_search = [
            SLAS for SLAS in seqreg.nodes_of_type(geffa.geffa.SLASNode)
            if int(SLAS.attributes["usage"]) >= min_usage
        ]
        PAS_to_search = [
            PAS for PAS in seqreg.nodes_of_type(geffa.geffa.PASNode)
            if int(PAS.attributes["usage"]) >= min_usage
        ]
        CDS_to_search_SLAS = [
            mRNA.CDS_children()[0]
            for mRNA in seqreg.nodes_of_type(geffa.geffa.MRNANode)
        ]
        CDS_to_search_PAS = [
            mRNA.CDS_children()[-1]
            for mRNA in seqreg.nodes_of_type(geffa.geffa.MRNANode)
        ]
        nodes_SLAS = PAS_to_search + CDS_to_search_SLAS
        nodes_PAS = SLAS_to_search + CDS_to_search_PAS
        for slas in (
            node
            for node in seqreg.node_registry.values()
            if node.type == 'SLAS'
        ):
            nodes_to_search = sorted(
                (
                    feature for feature in nodes_SLAS
                    if (
                        feature.type in ['CDS', 'PAS'] and
                        (
                            (
                                (slas.strand == '+') and
                                (feature.strand == '+') and
                                (feature.start > slas.end)
                            ) or
                            (
                                (slas.strand == '-') and
                                (feature.strand == '-') and
                                (feature.end < slas.start)
                            )
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

        for pas in (
            node
            for node in seqreg.node_registry.values()
            if node.type == 'PAS'
        ):
            nodes_to_search = sorted(
                (
                    feature for feature in nodes_PAS
                    if (
                        feature.type in ['CDS', 'SLAS'] and
                        (
                            (
                                (pas.strand == '+') and
                                (feature.strand == '+') and
                                (feature.end < pas.start)
                            ) or (
                                (pas.strand == '-') and
                                (feature.strand == '-') and
                                (feature.start > pas.end)
                            )
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


def identify_UTRs(annotations_gff: pathlib.Path, strip_existing: bool):
    gff = geffa.GffFile(annotations_gff, ignore_unknown_feature_types=True)
    check_or_strip_nodes(
        gff, ["five_prime_UTR", "three_prime_UTR"], strip_existing)

    for seqreg in gff.sequence_regions.values():
        for gene in [
            feature
            for feature in seqreg.node_registry.values()
            if feature.type == 'gene'
        ]:
            mRNAs = gene.children_of_type(geffa.geffa.MRNANode)
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
            mRNA: geffa.geffa.MRNANode = mRNAs[0]
            CDSs = mRNA.CDS_children()

            slas_sites: list[geffa.geffa.SLASNode] = []
            for slas in gene.children_of_type(geffa.geffa.SLASNode):
                if (
                    (mRNA.strand == '+' and CDSs[0].start < slas.end) or
                    (mRNA.strand == '-' and CDSs[0].end > slas.start)
                ):
                    logger.warning(
                        f"SLAS {slas.attributes['ID']} was assigned wrongly, "
                        "it is behind the start of the CDS. Skipping.")
                else:
                    slas_sites.append(slas)

            pas_sites: list[geffa.geffa.PASNode] = []
            for pas in gene.children_of_type(geffa.geffa.PASNode):
                if (
                    (mRNA.strand == '+' and CDSs[-1].end > pas.start) or
                    (mRNA.strand == '-' and CDSs[-1].start < pas.end)
                ):
                    logger.warning(
                        f"PAS {pas.attributes['ID']} was assigned wrongly, "
                        "it comes before the start of the CDS. Skipping.")
                else:
                    pas_sites.append(pas)

            if slas_sites:
                slas = sorted(slas_sites, key=lambda x: -
                              int(x.attributes['usage']))[0]
                cds = CDSs[0]
                if mRNA.strand == '+':
                    start = slas.end
                    end = cds.start - 1
                else:
                    end = slas.start
                    start = cds.end + 1
                geffa.geffa.FivePrimeUTRNode(
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
            if pas_sites:
                pas = sorted(pas_sites, key=lambda x: -
                             int(x.attributes['usage']))[0]
                cds = CDSs[-1]
                if mRNA.strand == '+':
                    start = cds.end + 1
                    end = pas.start
                else:
                    end = cds.start - 1
                    start = pas.end
                geffa.geffa.ThreePrimeUTRNode(
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

    return gff
