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
        def node_func(node):
            node.delete()
        logger.info("Stripping existing SLAS / PAS sites.")
    else:
        def node_func(_):
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
):
    gene_models = geffa.GffFile(
        gene_models_gff, ignore_unknown_feature_types=True)

    check_or_strip_nodes(gene_models, ["SLAS", "PAS"], strip_existing)

    slas_pas = geffa.GffFile(slas_pas_sites_gff)

    processed_regions = set()
    for seqreg in gene_models.sequence_regions.values():
        try:
            sitesreg = slas_pas.sequence_regions[seqreg.name]
        except KeyError:
            logger.warning(
                "SLAS/PAS sites GFF does not contain sequence region "
                f"{seqreg.name} - possible mismatch between slapquant "
                "run and given gene models?"
            )
            continue
        for SLAS in (
            feature
            for feature in sitesreg.node_registry.values()
            if feature.type == 'SLAS'
        ):
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
        for PAS in (
            feature
            for feature in sitesreg.node_registry.values()
            if feature.type == 'PAS'
        ):
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
        for slas in (
            node
            for node in seqreg.node_registry.values()
            if node.type == 'SLAS'
        ):
            nodes_to_search = sorted(
                (
                    feature for feature in seqreg.node_registry.values()
                    if (
                        feature.type in ['CDS', 'PAS'] and
                        (
                            (
                                (slas.strand == '+') and
                                (feature.strand == '+') and
                                (feature.start >= slas.end) or
                                (slas.strand == '-') and
                                (feature.strand == '-') and
                                (feature.end <= slas.start)
                            )
                        )
                    )
                ),
                key=lambda x: x.start if slas.strand == '+' else -x.end
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
                    feature for feature in seqreg.node_registry.values()
                    if (
                        feature.type in ['CDS', 'SLAS'] and
                        (
                            (
                                (pas.strand == '+') and
                                (feature.strand == '+') and
                                (feature.end <= pas.start) or
                                (pas.strand == '-') and
                                (feature.strand == '-') and
                                (feature.start >= pas.end)
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
            mRNAs = [
                feature
                for feature in gene.children
                if feature.type == 'mRNA'
            ]
            if len(mRNAs) == 0:
                logger.warning(
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

            slas_sites = [
                feature for feature in gene.children if feature.type == 'SLAS']
            pas_sites = [
                feature for feature in gene.children if feature.type == 'PAS']

            if slas_sites:
                slas = sorted(slas_sites, key=lambda x: -
                              int(x.attributes['usage']))[0]
                CDS = CDSs[0]
                if mRNA.strand == '+':
                    start = slas.end
                    end = CDS.start
                else:
                    end = slas.start
                    start = CDS.end
                geffa.geffa.FivePrimeUTRNode(
                    -1,
                    seqreg,
                    'RNASeq',
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
                CDS = CDSs[-1]
                if mRNA.strand == '+':
                    start = CDS.end
                    end = pas.start
                else:
                    end = CDS.start
                    start = pas.end
                geffa.geffa.ThreePrimeUTRNode(
                    -1,
                    seqreg,
                    'RNASeq',
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
        break

    return gff
