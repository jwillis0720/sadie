#!/usr/bin/env python
"""
Build IgBLAST internal_data directories for species.

Creates Ig/internal_data/{species}/ with:
1. Symlinks to BLAST databases from database/{species}/
2. Generated .ndm.imgt files from gapped sequences

Usage:
    python build_internal_data.py mouse rhesus_macaque chicken
"""

import argparse
import logging
import os
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# IMGT V gene regions (nucleotide positions in gapped sequence)
# These are the standard IMGT positions for a complete V gene
IMGT_V_REGIONS = {
    "FR1": (1, 78),      # 26 codons
    "CDR1": (79, 114),   # 12 codons (can have insertions)
    "FR2": (115, 165),   # 17 codons
    "CDR2": (166, 195),  # 10 codons (can have insertions)
    "FR3": (196, 312),   # 39 codons
}

# Chain type mapping
CHAIN_TYPES = {
    "H": "VH",
    "K": "VK", 
    "L": "VL",
}


def get_germlines_root() -> Path:
    """Get the germlines module root directory."""
    return Path(__file__).parent.parent


def calculate_ungapped_positions(gapped_seq: str) -> Dict[str, tuple]:
    """
    Calculate ungapped FR/CDR positions from an IMGT-gapped sequence.
    
    The IMGT gapping uses dots to maintain alignment. We need to calculate
    where each region falls in the ungapped sequence.
    
    Parameters
    ----------
    gapped_seq : str
        IMGT-gapped nucleotide sequence
        
    Returns
    -------
    Dict[str, tuple]
        Mapping of region name to (start, end) positions in ungapped sequence
    """
    # Remove gaps to get ungapped sequence
    ungapped_seq = gapped_seq.replace(".", "").replace("-", "")
    
    # Track position mapping: gapped_pos -> ungapped_pos
    ungapped_pos = 0
    pos_map = {}
    
    for gapped_pos, char in enumerate(gapped_seq, 1):
        if char not in ".-":
            ungapped_pos += 1
            pos_map[gapped_pos] = ungapped_pos
    
    # Calculate ungapped positions for each region
    regions = {}
    for region_name, (gapped_start, gapped_end) in IMGT_V_REGIONS.items():
        # Find the ungapped positions
        # Start: first non-gap position >= gapped_start
        # End: last non-gap position <= gapped_end
        
        start_pos = None
        end_pos = None
        
        for g_pos in range(gapped_start, min(gapped_end + 1, len(gapped_seq) + 1)):
            if g_pos in pos_map:
                if start_pos is None:
                    start_pos = pos_map[g_pos]
                end_pos = pos_map[g_pos]
        
        if start_pos is not None and end_pos is not None:
            regions[region_name] = (start_pos, end_pos)
    
    return regions, len(ungapped_seq)


def generate_ndm_entry(gene_name: str, gapped_seq: str, chain: str) -> Optional[str]:
    """
    Generate an NDM file entry for a V gene.
    
    NDM format:
    GENE_NAME  FR1_START  FR1_END  CDR1_START  CDR1_END  FR2_START  FR2_END  CDR2_START  CDR2_END  FR3_START  SEQ_LEN  CHAIN_TYPE  FLAG
    
    Parameters
    ----------
    gene_name : str
        Gene name (e.g., IGHV1-2*01)
    gapped_seq : str
        IMGT-gapped nucleotide sequence
    chain : str
        Chain type (H, K, L)
        
    Returns
    -------
    str or None
        NDM entry line or None if cannot be calculated
    """
    if not gapped_seq or "." not in gapped_seq:
        return None
    
    regions, seq_len = calculate_ungapped_positions(gapped_seq)
    
    # Need at least FR1 through FR3 start
    required = ["FR1", "CDR1", "FR2", "CDR2", "FR3"]
    if not all(r in regions for r in required):
        return None
    
    chain_type = CHAIN_TYPES.get(chain, "VH")
    
    # Format: gene  fr1_start  fr1_end  cdr1_start  cdr1_end  fr2_start  fr2_end  cdr2_start  cdr2_end  fr3_start  seq_len  chain_type  flag
    entry = (
        f"{gene_name}\t"
        f"{regions['FR1'][0]}\t{regions['FR1'][1]}\t"
        f"{regions['CDR1'][0]}\t{regions['CDR1'][1]}\t"
        f"{regions['FR2'][0]}\t{regions['FR2'][1]}\t"
        f"{regions['CDR2'][0]}\t{regions['CDR2'][1]}\t"
        f"{regions['FR3'][0]}\t{seq_len}\t"
        f"{chain_type}\t0"
    )
    
    return entry


def build_ndm_file(species: str, germlines_root: Path) -> List[str]:
    """
    Build NDM file content for a species from IMGT gapped sequences.
    
    Parameters
    ----------
    species : str
        Species name
    germlines_root : Path
        Germlines module root directory
        
    Returns
    -------
    List[str]
        NDM file lines
    """
    imgt_dir = germlines_root / "sources" / "imgt" / species
    entries = []
    
    # Process V genes for each chain
    for chain in ["H", "K", "L"]:
        gapped_fasta = imgt_dir / f"IG{chain}V_gapped.fasta"
        
        if not gapped_fasta.exists():
            logger.debug(f"No gapped V file: {gapped_fasta}")
            continue
        
        for record in SeqIO.parse(gapped_fasta, "fasta"):
            # Parse IMGT header: >ACCESSION|GENE_NAME|SPECIES|...
            parts = record.id.split("|")
            gene_name = parts[1] if len(parts) > 1 else parts[0]
            
            gapped_seq = str(record.seq).upper()
            entry = generate_ndm_entry(gene_name, gapped_seq, chain)
            
            if entry:
                entries.append(entry)
    
    return entries


def build_internal_data(species: str, germlines_root: Path) -> bool:
    """
    Build internal_data directory for a species.
    
    Creates:
    - Ig/internal_data/{species}/ directory
    - Symlinks to BLAST databases
    - {species}.ndm.imgt file
    
    Parameters
    ----------
    species : str
        Species name
    germlines_root : Path
        Germlines module root directory
        
    Returns
    -------
    bool
        True if successful
    """
    database_dir = germlines_root / "igblast" / "database" / species
    internal_data_dir = germlines_root / "igblast" / "Ig" / "internal_data" / species
    
    # Check database exists
    if not database_dir.exists():
        logger.error(f"Database directory not found: {database_dir}")
        return False
    
    # Create internal_data directory
    internal_data_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created: {internal_data_dir}")
    
    # Symlink BLAST database files
    for db_file in database_dir.glob(f"{species}_*"):
        link_path = internal_data_dir / db_file.name
        if link_path.exists() or link_path.is_symlink():
            link_path.unlink()
        
        # Use relative symlink
        rel_path = os.path.relpath(db_file, internal_data_dir)
        link_path.symlink_to(rel_path)
    
    logger.info(f"Symlinked BLAST databases from {database_dir}")
    
    # Generate NDM file
    ndm_entries = build_ndm_file(species, germlines_root)
    
    if not ndm_entries:
        logger.warning(f"No NDM entries generated for {species}")
        return False
    
    ndm_path = internal_data_dir / f"{species}.ndm.imgt"
    with open(ndm_path, "w") as f:
        f.write("\n".join(ndm_entries) + "\n")
    
    logger.info(f"Generated {ndm_path} with {len(ndm_entries)} entries")
    
    return True


def main():
    parser = argparse.ArgumentParser(description="Build IgBLAST internal_data for species")
    parser.add_argument("species", nargs="+", help="Species names to build")
    args = parser.parse_args()
    
    germlines_root = get_germlines_root()
    
    for species in args.species:
        logger.info(f"\n=== Building internal_data for {species} ===")
        success = build_internal_data(species, germlines_root)
        if success:
            logger.info(f"SUCCESS: {species} internal_data ready")
        else:
            logger.error(f"FAILED: {species}")


if __name__ == "__main__":
    main()
