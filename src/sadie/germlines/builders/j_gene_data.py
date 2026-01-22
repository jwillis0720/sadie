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
        "K": (1, "JK", 6, 1),   # Most IGKJ are RF1, CDR3 6
        "L": (1, "JL", 6, 1),   # Most IGLJ are RF1, CDR3 6
    }

    return defaults.get(chain, (1, chain_type, 10, 1))
