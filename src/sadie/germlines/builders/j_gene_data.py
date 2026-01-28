"""
J Gene Reference Data
=====================

Parses J gene sequences to extract CDR3 end positions, reading frames, and
FWR4 boundaries using species-specific amino acid motif patterns.

This module provides full parity with the G3 database's J gene parsing logic,
using the same motif lookup patterns to identify the conserved FWR4 start.

The FWR4 region in J genes starts with a conserved motif:
- Heavy chain: typically W-G-x-G (Trp-Gly-any-Gly)
- Kappa chain: typically F-G-x-G (Phe-Gly-any-Gly)
- Lambda chain: typically F-G-x-G (Phe-Gly-any-Gly)

These patterns vary by species and are defined in j_gene_motif.json.
"""

import json
import logging
import re
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

logger = logging.getLogger(__name__)

# Load motif lookup from JSON file
_MOTIF_FILE = Path(__file__).parent / "data" / "j_gene_motif.json"
try:
    MOTIF_LOOKUP: Dict = json.loads(_MOTIF_FILE.read_text())
except FileNotFoundError:
    logger.warning(f"Motif file not found: {_MOTIF_FILE}")
    MOTIF_LOOKUP = {}

# Human J gene reference data (validated fallback)
# Format: {allele: (reading_frame, cdr3_end, extra_bps)}
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

# Species name normalization (handles variations in naming)
SPECIES_ALIASES = {
    "rhesus_monkey": "macaque",
    "rhesus": "macaque",
    "cynomolgus_monkey": "cynomolgus",
    "crab-eating_macaque": "crab_eating_macaque",
    "sus_scrofa": "pig",
    "bos_taurus": "cow",
    "canis_familiaris": "dog",
    "felis_catus": "cat",
    "oryctolagus_cuniculus": "rabbit",
    "mus_musculus": "mouse",
    "rattus_norvegicus": "rat",
    "homo_sapiens": "human",
    "danio_rerio": "zebrafish",
    "oncorhynchus_mykiss": "rainbow_trout",
    "salmo_salar": "atlantic_salmon",
    "gallus_gallus": "chicken",
}


def _normalize_species(species: str) -> str:
    """Normalize species name to match motif lookup keys."""
    normalized = species.lower().replace(" ", "_").replace("-", "_")
    return SPECIES_ALIASES.get(normalized, normalized)


def _get_gene_type(gene_name: str) -> str:
    """
    Extract gene type from gene name (e.g., 'IGHJ1*01' -> 'IGHJ').

    Parameters
    ----------
    gene_name : str
        Full gene name with allele

    Returns
    -------
    str
        Gene type (first 4 characters, e.g., 'IGHJ', 'IGKJ', 'IGLJ')
    """
    return gene_name[:4].upper()


def _get_short_name(gene_name: str) -> str:
    """
    Extract short gene name without allele (e.g., 'IGHJ1*01' -> 'IGHJ1').

    Parameters
    ----------
    gene_name : str
        Full gene name with allele

    Returns
    -------
    str
        Short name without allele
    """
    return gene_name.split("*")[0]


def _translate_sequence(sequence: str, reading_frame: int = 0) -> str:
    """
    Translate nucleotide sequence to amino acids.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence
    reading_frame : int
        Reading frame offset (0, 1, or 2)

    Returns
    -------
    str
        Amino acid sequence
    """
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }

    seq = sequence.upper()[reading_frame:]
    aa_seq = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = codon_table.get(codon, "X")
        aa_seq.append(aa)
    return "".join(aa_seq)


def _get_reading_frames(sequence: str) -> Dict[int, str]:
    """
    Translate sequence in all three reading frames.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence

    Returns
    -------
    Dict[int, str]
        Mapping of reading frame (0, 1, 2) to amino acid sequence
    """
    return {rf: _translate_sequence(sequence, rf) for rf in [0, 1, 2]}


def _get_longest_frame(reading_frames: Dict[int, str]) -> int:
    """
    Find reading frame with longest sequence before stop codon.

    Parameters
    ----------
    reading_frames : Dict[int, str]
        Mapping of reading frame to amino acid sequence

    Returns
    -------
    int
        Best reading frame (0, 1, or 2)
    """
    longest = 0
    longest_rf = 0
    for rf, aa_seq in reading_frames.items():
        if "*" in aa_seq:
            stop_index = aa_seq.index("*")
        else:
            stop_index = len(aa_seq)
        if stop_index > longest:
            longest_rf = rf
            longest = stop_index
    return longest_rf


def _infer_reading_frame(sequence: str, gene_name: str) -> int:
    """
    Infer the correct reading frame for a J gene.

    Uses the last amino acid as a hint:
    - Heavy chain (H): typically ends with S
    - Kappa chain (K): typically ends with K
    - Lambda chain (L): typically ends with L

    Parameters
    ----------
    sequence : str
        Nucleotide sequence
    gene_name : str
        Gene name to determine chain type

    Returns
    -------
    int
        Inferred reading frame (0, 1, or 2)
    """
    potential_frames = _get_reading_frames(sequence)

    # Determine expected last amino acid based on chain
    chain_char = gene_name[2].upper() if len(gene_name) > 2 else ""
    look_for = {"H": "S", "K": "K", "L": "L"}.get(chain_char, "")

    if look_for:
        # Find frames ending with the expected amino acid
        matching_frames = {
            rf: aa_seq for rf, aa_seq in potential_frames.items()
            if aa_seq and aa_seq[-1] == look_for
        }
        if len(matching_frames) == 1:
            return list(matching_frames.keys())[0]

    # Fallback: use longest frame before stop codon
    return _get_longest_frame(potential_frames)


def parse_j_gene(
    species: str,
    gene_name: str,
    sequence: str,
    reading_frame: Optional[int] = None
) -> Dict:
    """
    Parse J gene sequence to extract CDR3/FWR4 boundaries.

    Uses species-specific motif patterns from the motif lookup to identify
    the conserved FWR4 start position.

    Parameters
    ----------
    species : str
        Species name (e.g., "human", "macaque", "mouse")
    gene_name : str
        Full gene name (e.g., "IGHJ1*01")
    sequence : str
        J gene nucleotide sequence (ungapped)
    reading_frame : int, optional
        Reading frame (0, 1, or 2). If None, will be inferred.

    Returns
    -------
    Dict
        Dictionary containing:
        - sequence_gapped_aa: Amino acid sequence
        - cdr3: CDR3 nucleotide sequence
        - cdr3_aa: CDR3 amino acid sequence
        - cdr3_start: CDR3 start position (always 0 for J genes)
        - cdr3_end: CDR3 end position (nucleotide, 0-based)
        - fwr4: FWR4 nucleotide sequence
        - fwr4_aa: FWR4 amino acid sequence
        - fwr4_start: FWR4 start position
        - fwr4_end: FWR4 end position
        - reading_frame: Reading frame used
        - remainder: Extra nucleotides beyond last complete codon
        - ignored: Whether gene is in ignore list
        - not_implemented: Whether species/gene type not in motif lookup
        - expression_match: Whether motif was found
        - expression: List of motif matches found
    """
    sequence = sequence.upper().replace(".", "").replace("-", "")
    species_normalized = _normalize_species(species)

    # Initialize return values
    cdr3_aa = ""
    cdr3_nt = ""
    fwr4_aa = ""
    fwr4_nt = ""
    ignored = False
    not_implemented = False
    cdr3_end_nt_index = None
    fwr4_start = None
    fwr4_end = None
    expression_match = False
    remainder = None
    expression: List[str] = []

    # Determine reading frame
    if reading_frame is None:
        reading_frame = _infer_reading_frame(sequence, gene_name)

    # Translate to amino acids
    j_gene_aa = _translate_sequence(sequence, reading_frame)

    # Check if species is in motif lookup
    if species_normalized in MOTIF_LOOKUP:
        short_name = _get_short_name(gene_name)
        gene_type = _get_gene_type(gene_name)
        species_data = MOTIF_LOOKUP[species_normalized]

        if gene_type in species_data:
            motif_pattern = species_data[gene_type]
            ignore_list = species_data.get("ignore", [])

            # Check if gene should be ignored
            if isinstance(ignore_list, list) and short_name in ignore_list:
                ignored = True
            elif isinstance(ignore_list, str) and ignore_list and short_name == ignore_list:
                ignored = True

            if not ignored:
                # Search for motif in amino acid sequence
                matches = re.findall(motif_pattern, j_gene_aa)
                expression = matches

                if matches:
                    end_motif = matches[0]
                    end_cdr3_index = j_gene_aa.index(end_motif)
                    end_index_codon = end_cdr3_index * 3 + reading_frame

                    cdr3_nt = sequence[:end_index_codon]
                    remainder = len(sequence[reading_frame:]) % 3
                    if remainder:
                        fwr4_nt = sequence[end_index_codon:-remainder]
                    else:
                        fwr4_nt = sequence[end_index_codon:]

                    cdr3_aa = j_gene_aa[:end_cdr3_index]
                    fwr4_aa = j_gene_aa[end_cdr3_index:]
                    cdr3_end_nt_index = end_index_codon
                    expression_match = True
                    fwr4_start = cdr3_end_nt_index
                    fwr4_end = len(sequence) - (remainder or 0)
        else:
            not_implemented = True
    else:
        not_implemented = True

    return {
        "sequence_gapped_aa": j_gene_aa,
        "cdr3": cdr3_nt,
        "cdr3_aa": cdr3_aa,
        "cdr3_start": 0,
        "cdr3_end": cdr3_end_nt_index,
        "fwr4": fwr4_nt,
        "fwr4_aa": fwr4_aa,
        "fwr4_start": fwr4_start,
        "fwr4_end": fwr4_end,
        "reading_frame": reading_frame,
        "ignored": ignored,
        "remainder": remainder,
        "not_implemented": not_implemented,
        "expression_match": expression_match,
        "expression": expression,
    }


def get_j_gene_data(
    allele_name: str,
    chain: str,
    sequence: Optional[str] = None,
    species: Optional[str] = None
) -> Tuple[int, str, int, int]:
    """
    Get J gene reference data for an allele.

    Uses species-specific motif patterns to calculate CDR3 end position
    from the sequence. Falls back to human reference data or defaults
    if sequence is not provided or parsing fails.

    Parameters
    ----------
    allele_name : str
        Full allele name (e.g., "IGHJ1*01")
    chain : str
        Chain type (H, K, or L)
    sequence : str, optional
        J gene nucleotide sequence for dynamic calculation
    species : str, optional
        Species name for motif lookup (defaults to "human")

    Returns
    -------
    tuple
        (reading_frame, chain_type, cdr3_end, extra_bps)
        - reading_frame: First coding frame start position (0, 1, or 2)
        - chain_type: Chain type string (JH, JK, JL)
        - cdr3_end: CDR3 stop position (0-based)
        - extra_bps: Extra base pairs beyond J coding end (typically 1)
    """
    chain_type = CHAIN_TYPE_MAP.get(chain, f"J{chain}")

    # Check known human reference data first (for validated human genes)
    if allele_name in HUMAN_J_GENE_DATA and (species is None or species.lower() == "human"):
        rf, cdr3_end, extra_bps = HUMAN_J_GENE_DATA[allele_name]
        return (rf, chain_type, cdr3_end, extra_bps)

    # Calculate from sequence if provided
    if sequence:
        # Use species if provided, otherwise try to infer or default to human
        parse_species = species or "human"
        result = parse_j_gene(parse_species, allele_name, sequence)

        if result["expression_match"] and result["cdr3_end"] is not None:
            cdr3_end = result["cdr3_end"]
            reading_frame = result["reading_frame"]
            remainder = result["remainder"] or 1
            return (reading_frame, chain_type, cdr3_end, remainder)

        # If motif not found but we have the sequence, use reading frame calculation
        if not result["not_implemented"]:
            reading_frame = result["reading_frame"]
            # Default CDR3 end based on chain type if motif not found
            defaults = {"H": 15, "K": 6, "L": 6}
            cdr3_end = defaults.get(chain, 10)
            return (reading_frame, chain_type, cdr3_end, 1)

    # Default fallback values based on chain type
    defaults = {
        "H": (1, "JH", 15, 1),
        "K": (1, "JK", 6, 1),
        "L": (1, "JL", 6, 1),
    }

    return defaults.get(chain, (1, chain_type, 10, 1))
