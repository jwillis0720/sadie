"""
J Gene Reference Data
=====================

Contains CDR3 end positions and reading frames for J gene alleles.
Data sourced from IMGT reference and G3 aux files.
"""

# Human J gene reference data
# Format: {allele: (reading_frame, cdr3_end, is_functional)}
HUMAN_J_GENE_DATA = {
    # Heavy chain (JH)
    "IGHJ1*01": (0, 17, 1),
    "IGHJ2*01": (1, 18, 1),
    "IGHJ3*01": (1, 15, 1),
    "IGHJ3*02": (1, 15, 1),
    "IGHJ4*01": (2, 13, 1),
    "IGHJ4*02": (2, 13, 1),
    "IGHJ4*03": (2, 13, 1),
    "IGHJ5*01": (2, 16, 1),
    "IGHJ5*02": (2, 16, 1),
    "IGHJ5*03": (2, 16, 1),
    "IGHJ5*04": (2, 16, 1),
    "IGHJ6*01": (2, 28, 1),
    "IGHJ6*02": (2, 28, 0),
    "IGHJ6*03": (2, 28, 0),
    "IGHJ6*04": (2, 28, 1),
    # Kappa chain (JK)
    "IGKJ1*01": (1, 6, 1),
    "IGKJ2*01": (2, 7, 1),
    "IGKJ2*02": (1, 6, 1),
    "IGKJ2*03": (2, 7, 1),
    "IGKJ2*04": (2, 7, 1),
    "IGKJ3*01": (1, 6, 1),
    "IGKJ4*01": (1, 6, 1),
    "IGKJ4*02": (1, 6, 1),
    "IGKJ4*03": (1, 6, 1),
    "IGKJ5*01": (1, 6, 1),
    # Lambda chain (JL)
    "IGLJ1*01": (1, 6, 1),
    "IGLJ2*01": (1, 6, 1),
    "IGLJ2A*01": (1, 6, 1),
    "IGLJ3*01": (1, 6, 1),
    "IGLJ3*02": (1, 6, 1),
    "IGLJ4*01": (1, 6, 1),
    "IGLJ5*01": (1, 6, 1),
    "IGLJ5*02": (1, 6, 1),
    "IGLJ6*01": (1, 6, 1),
    "IGLJ7*01": (1, 6, 1),
    "IGLJ7*02": (1, 6, 1),
}

# Chain type mapping
CHAIN_TYPE_MAP = {
    "H": "JH",
    "K": "JK",
    "L": "JL",
}


# J gene sequence lengths (nucleotides) for complete_vdj calculation
# Sourced from SADIE germline database FASTA files
J_GENE_LENGTHS = {
    # Heavy chain
    "IGHJ1*01": 52,
    "IGHJ2*01": 53,
    "IGHJ3*01": 50,
    "IGHJ3*02": 50,
    "IGHJ4*01": 48,
    "IGHJ4*02": 48,
    "IGHJ4*03": 48,
    "IGHJ5*01": 51,
    "IGHJ5*02": 51,
    "IGHJ5*03": 51,
    "IGHJ5*04": 51,
    "IGHJ6*01": 63,
    "IGHJ6*02": 63,
    "IGHJ6*03": 63,
    "IGHJ6*04": 63,
    # Kappa chain
    "IGKJ1*01": 38,
    "IGKJ2*01": 39,
    "IGKJ2*02": 38,
    "IGKJ2*03": 39,
    "IGKJ2*04": 39,
    "IGKJ3*01": 38,
    "IGKJ4*01": 38,
    "IGKJ4*02": 38,
    "IGKJ4*03": 38,
    "IGKJ5*01": 38,
    # Lambda chain
    "IGLJ1*01": 38,
    "IGLJ2*01": 38,
    "IGLJ3*02": 38,
    "IGLJ4*01": 38,
    "IGLJ5*01": 38,
    "IGLJ5*02": 38,
    "IGLJ6*01": 38,
    "IGLJ7*01": 38,
    "IGLJ7*02": 38,
}


def get_j_gene_length(allele_name: str) -> int | None:
    """Get expected J gene length for complete_vdj calculation.

    Parameters
    ----------
    allele_name : str
        Full allele name (e.g., "IGHJ1*01")

    Returns
    -------
    int | None
        Expected J gene length in nucleotides, or None if not found
    """
    if allele_name in J_GENE_LENGTHS:
        return J_GENE_LENGTHS[allele_name]
    # Try gene family (remove allele suffix, use *01)
    gene = allele_name.split("*")[0] + "*01"
    return J_GENE_LENGTHS.get(gene)


def get_j_gene_data(allele_name: str, chain: str) -> tuple:
    """
    Get J gene reference data for an allele.

    Parameters
    ----------
    allele_name : str
        Full allele name (e.g., "IGHJ1*01")
    chain : str
        Chain type (H, K, or L)

    Returns
    -------
    tuple
        (reading_frame, chain_type, cdr3_end, is_functional)
        Returns default values if allele not found.
    """
    chain_type = CHAIN_TYPE_MAP.get(chain, f"J{chain}")

    # Check known reference data
    if allele_name in HUMAN_J_GENE_DATA:
        rf, cdr3_end, is_func = HUMAN_J_GENE_DATA[allele_name]
        return (rf, chain_type, cdr3_end, is_func)

    # Default fallback values based on chain type
    # These defaults are based on most common values
    defaults = {
        "H": (1, "JH", 15, 1),  # Most IGHJ are RF1, CDR3 ~15
        "K": (1, "JK", 6, 1),  # Most IGKJ are RF1, CDR3 6
        "L": (1, "JL", 6, 1),  # Most IGLJ are RF1, CDR3 6
    }

    return defaults.get(chain, (1, chain_type, 10, 1))
