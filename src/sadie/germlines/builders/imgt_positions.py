"""
IMGT Region Position Utilities
==============================

Helpers for deriving ungapped IMGT region positions from gapped V-gene sequences.
Positions are 0-based (inclusive) by default to match G3 API output.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

# IMGT V gene regions in gapped coordinates (1-based, inclusive)
IMGT_V_REGIONS = {
    "FWR1": (1, 78),
    "CDR1": (79, 114),
    "FWR2": (115, 165),
    "CDR2": (166, 195),
    "FWR3": (196, 312),
}

GAP_CHARACTERS = {".", "-"}


def calculate_imgt_v_positions(gapped_sequence: str, *, zero_based: bool = True) -> Dict[str, Tuple[int, int]]:
    """
    Calculate ungapped IMGT region positions for V genes.

    Parameters
    ----------
    gapped_sequence : str
        IMGT-gapped nucleotide sequence.
    zero_based : bool, optional
        If True, return 0-based positions (default). If False, return 1-based positions.

    Returns
    -------
    Dict[str, Tuple[int, int]]
        Mapping of region name to (start, end) positions in ungapped sequence.
    """
    if not gapped_sequence:
        return {}

    pos_map: Dict[int, int] = {}
    ungapped_pos = 0

    for gapped_pos, char in enumerate(gapped_sequence, start=1):
        if char not in GAP_CHARACTERS:
            pos_map[gapped_pos] = ungapped_pos
            ungapped_pos += 1

    regions: Dict[str, Tuple[int, int]] = {}
    for region_name, (gapped_start, gapped_end) in IMGT_V_REGIONS.items():
        start_pos = None
        end_pos = None

        for gapped_pos in range(gapped_start, min(gapped_end, len(gapped_sequence)) + 1):
            if gapped_pos in pos_map:
                if start_pos is None:
                    start_pos = pos_map[gapped_pos]
                end_pos = pos_map[gapped_pos]

        if start_pos is not None and end_pos is not None:
            if zero_based:
                regions[region_name.lower()] = (start_pos, end_pos)
            else:
                regions[region_name.lower()] = (start_pos + 1, end_pos + 1)

    return regions


def derive_imgt_v_regions(
    gapped_sequence: str,
    ungapped_sequence: Optional[str] = None,
    *,
    zero_based: bool = True,
) -> Tuple[Dict[str, Tuple[int, int]], Dict[str, str]]:
    """
    Derive IMGT V-region positions and sequences from a gapped sequence.

    Parameters
    ----------
    gapped_sequence : str
        IMGT-gapped nucleotide sequence.
    ungapped_sequence : str, optional
        Ungapped nucleotide sequence. Derived from gapped_sequence if not provided.
    zero_based : bool, optional
        If True, return 0-based positions (default). If False, return 1-based positions.

    Returns
    -------
    Tuple[Dict[str, Tuple[int, int]], Dict[str, str]]
        Region positions and region nucleotide sequences.
    """
    positions = calculate_imgt_v_positions(gapped_sequence, zero_based=zero_based)
    if not positions:
        return {}, {}

    if ungapped_sequence is None:
        ungapped_sequence = gapped_sequence.replace(".", "").replace("-", "")
    ungapped_sequence = ungapped_sequence.upper()

    regions: Dict[str, str] = {}
    for region_name, (start, end) in positions.items():
        if zero_based:
            start_index = start
            end_index = end + 1
        else:
            start_index = max(start - 1, 0)
            end_index = end

        if start_index < len(ungapped_sequence):
            regions[region_name] = ungapped_sequence[start_index:end_index]

    return positions, regions
