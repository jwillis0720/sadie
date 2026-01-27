#!/usr/bin/env python
"""
Build IgBLAST internal_data directories for species.

Creates Ig/internal_data/{species}/ with:
1. Combined VDJC FASTA file ({species}_V.fasta) containing all segments
2. BLAST database built from combined FASTA
3. Generated .ndm.imgt files from gapped sequences

The combined VDJC file is required for IgBLAST's complete_vdj calculation.
IgBLAST expects internal_data/{species}/{species}_V to contain ALL gene segments
(V, D, J, C) for proper annotation.

Usage:
    python build_internal_data.py mouse rhesus_macaque chicken
"""

import argparse
import logging
import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO

from sadie.germlines.builders.imgt_positions import calculate_imgt_v_positions

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# Chain type mapping
CHAIN_TYPES = {
    "H": "VH",
    "K": "VK",
    "L": "VL",
}


def build_combined_fasta(species: str, database_dir: Path, internal_data_dir: Path) -> Path:
    """
    Concatenate V+D+J+C FASTAs into single combined file with deduplication.

    IgBLAST requires a single file named {species}_V.fasta in internal_data/
    that contains ALL gene segments (V, D, J, C) for proper complete_vdj calculation.

    Sequences are deduplicated by ID to avoid BLAST database errors when the
    same gene appears in multiple segment files.

    Parameters
    ----------
    species : str
        Species name
    database_dir : Path
        Path to database/{species}/ containing separate segment FASTAs
    internal_data_dir : Path
        Path to internal_data/{species}/ for output

    Returns
    -------
    Path
        Path to combined FASTA file
    """
    combined_fasta = internal_data_dir / f"{species}_V.fasta"

    # Ordered segments: V, D, J, C (D may not exist for all species)
    segments = ["V", "D", "J", "C"]

    # Track seen IDs to deduplicate
    seen_ids: set = set()
    total_sequences = 0
    duplicates_skipped = 0

    with open(combined_fasta, "w") as out:
        for segment in segments:
            fasta_file = database_dir / f"{species}_{segment}.fasta"
            if fasta_file.exists():
                segment_count = 0
                for record in SeqIO.parse(fasta_file, "fasta"):
                    seq_id = record.id
                    if seq_id in seen_ids:
                        duplicates_skipped += 1
                        continue
                    seen_ids.add(seq_id)
                    out.write(f">{seq_id}\n{str(record.seq)}\n")
                    segment_count += 1
                    total_sequences += 1
                logger.info(f"  Added {segment_count} {segment} genes from {fasta_file.name}")
            else:
                logger.debug(f"  No {segment} genes found at {fasta_file}")

    if duplicates_skipped > 0:
        logger.info(f"  Skipped {duplicates_skipped} duplicate sequences")
    logger.info(f"  Combined FASTA: {total_sequences} unique sequences")
    return combined_fasta


def build_blast_db(fasta_path: Path, db_prefix: Path) -> bool:
    """
    Build BLAST database from FASTA file.

    Parameters
    ----------
    fasta_path : Path
        Path to input FASTA file
    db_prefix : Path
        Output database prefix (e.g., internal_data/human/human_V)

    Returns
    -------
    bool
        True if successful
    """
    # Find makeblastdb binary
    system = platform.system().lower()
    reference_bin_dir = Path(__file__).parent.parent.parent / "reference" / "bin" / system
    make_blast_db_bin = reference_bin_dir / "makeblastdb"

    if not shutil.which(str(make_blast_db_bin)):
        # Try system makeblastdb
        make_blast_db_bin = shutil.which("makeblastdb")  # type: ignore
        if not make_blast_db_bin:
            logger.error("makeblastdb not found - ensure BLAST+ is installed")
            return False

    cmd = [
        str(make_blast_db_bin),
        "-in",
        str(fasta_path),
        "-dbtype",
        "nucl",
        "-hash_index",
        "-parse_seqids",
        "-out",
        str(db_prefix),
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  Built BLAST database at {db_prefix}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"makeblastdb failed: {e.stderr}")
        return False
    except FileNotFoundError:
        logger.error("makeblastdb not found - ensure BLAST+ is installed")
        return False


def get_germlines_root() -> Path:
    """Get the germlines module root directory."""
    return Path(__file__).parent.parent


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

    regions = calculate_imgt_v_positions(gapped_seq, zero_based=False)

    # Need at least FR1 through FR3 start
    required = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"]
    if not all(r in regions for r in required):
        return None

    chain_type = CHAIN_TYPES.get(chain, "VH")

    # Format: gene  fr1_start  fr1_end  cdr1_start  cdr1_end  fr2_start  fr2_end  cdr2_start  cdr2_end  fr3_start  seq_len  chain_type  flag
    entry = (
        f"{gene_name}\t"
        f"{regions['fwr1'][0]}\t{regions['fwr1'][1]}\t"
        f"{regions['cdr1'][0]}\t{regions['cdr1'][1]}\t"
        f"{regions['fwr2'][0]}\t{regions['fwr2'][1]}\t"
        f"{regions['cdr2'][0]}\t{regions['cdr2'][1]}\t"
        f"{regions['fwr3'][0]}\t{regions['fwr3'][1]}\t"
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
    - Combined VDJC FASTA ({species}_V.fasta) containing all segments
    - BLAST database built from combined FASTA
    - {species}.ndm.imgt file (V genes only, for region annotations)

    The combined VDJC file is required for IgBLAST's complete_vdj calculation.
    IgBLAST reads {species}_V from internal_data/ for its internal annotation
    database, which must contain ALL gene segments.

    Note: The separate V/D/J/C BLAST databases in database/{species}/ are used
    by IgBLAST's -germline_db_V, -germline_db_D, etc. options for alignment.

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

    # Clean internal_data directory - remove old symlinks and files
    if internal_data_dir.exists():
        for item in internal_data_dir.iterdir():
            # Remove symlinks and old segment files (but keep .ndm.imgt for now)
            if item.is_symlink() or item.name.startswith(f"{species}_"):
                item.unlink()
                logger.debug(f"  Removed: {item.name}")

    # Create internal_data directory
    internal_data_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created: {internal_data_dir}")

    # Create combined VDJC FASTA
    combined_fasta = build_combined_fasta(species, database_dir, internal_data_dir)
    if not combined_fasta.exists() or combined_fasta.stat().st_size == 0:
        logger.error(f"Failed to create combined FASTA for {species}")
        return False

    # Build BLAST database from combined file
    db_prefix = internal_data_dir / f"{species}_V"
    if not build_blast_db(combined_fasta, db_prefix):
        logger.error(f"Failed to build BLAST database for {species}")
        return False

    # Generate NDM file (V genes only - for framework/CDR region annotations)
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
