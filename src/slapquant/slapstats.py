import logging
import pathlib
from geffa.geffa import GffFile, GeneNode, SLASNode, PASNode
import pandas as pd

logger = logging.getLogger('slapstats')


def get_usage(node: SLASNode | PASNode):
    try:
        return int(node.attributes["usage"])
    except KeyError:
        logger.warning(
            "SLAS site has no attribute 'usage', substituting zero"
            "(suppressing further warnings of this type)."
        )
    return 0


def gene_stats(gene: GeneNode):
    SLAS_nodes = gene.children_of_type(SLASNode)
    PAS_nodes = gene.children_of_type(PASNode)
    return {
        "SLAS_site_count": len(SLAS_nodes),
        "SLAS_site_total_usage": sum(get_usage(s) for s in SLAS_nodes),
        "PAS_site_count": len(PAS_nodes),
        "PAS_site_total_usage": sum(get_usage(p) for p in PAS_nodes),
    }


def gather_stats(gff: GffFile | pathlib.Path):
    if not isinstance(gff, GffFile):
        gff = GffFile(gff, ignore_unknown_feature_types=True)

    data = []
    index = []
    for seqreg in gff.sequence_regions.values():
        for gene in seqreg.nodes_of_type(GeneNode):
            index.append((gene.sequence_region.name, gene.attributes["ID"]))
            data.append(gene_stats(gene))

    return pd.DataFrame(
        data,
        index=pd.MultiIndex.from_tuples(index, names=("contig", "gene"))
    ).sort_index()
