#!/usr/bin/env python3
"""
Build Auxiliary Files for IgBLAST
=================================

Script to generate *.aux files for IgBLAST from normalized germline data.

This is a CLI wrapper around the existing AuxFileBuilder in builders/aux.py.

Usage:
    python -m sadie.germlines.scripts.build_aux_files --species human
    python -m sadie.germlines.scripts.build_aux_files --all-species
    python -m sadie.germlines.scripts.build_aux_files --help

Output:
    src/sadie/germlines/igblast/aux_db/{species}_gl.aux
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger(__name__)

# Module paths
GERMLINES_ROOT = Path(__file__).parent.parent
NORMALIZED_DIR = GERMLINES_ROOT / "normalized"
AUX_OUTPUT_DIR = GERMLINES_ROOT / "igblast" / "aux_db"
DATABASE_DIR = GERMLINES_ROOT / "igblast" / "database"


def get_available_species() -> List[str]:
    """
    Get list of species with built BLAST databases.
    
    Returns
    -------
    List[str]
        Species names with databases
    """
    species = []
    if DATABASE_DIR.exists():
        for d in DATABASE_DIR.iterdir():
            if d.is_dir() and not d.name.startswith('.'):
                # Check if it has actual database files
                v_db = d / f"{d.name}_V.nsq"
                if v_db.exists():
                    species.append(d.name)
    return sorted(species)


def get_normalized_species() -> List[str]:
    """
    Get list of species with normalized data.
    
    Returns
    -------
    List[str]
        Species names with normalized data
    """
    species = []
    if NORMALIZED_DIR.exists():
        for d in NORMALIZED_DIR.iterdir():
            if d.is_dir() and not d.name.startswith('.'):
                gapped_dir = d / "gapped"
                if gapped_dir.exists() and list(gapped_dir.glob("*.fasta")):
                    species.append(d.name)
    return sorted(species)


def build_aux_file(species: str, force: bool = False) -> Path:
    """
    Build aux file for a specific species using the AuxFileBuilder.
    
    Parameters
    ----------
    species : str
        Species name
    force : bool
        Force regeneration even if file exists
        
    Returns
    -------
    Path
        Path to generated aux file
    """
    from sadie.germlines.builders.aux import AuxFileBuilder
    
    output_file = AUX_OUTPUT_DIR / f"{species}_gl.aux"
    source_dir = NORMALIZED_DIR / species / "gapped"
    
    if output_file.exists() and not force:
        logger.info(f"Aux file already exists: {output_file}")
        return output_file
    
    if not source_dir.exists():
        logger.warning(f"No normalized gapped data for {species}")
        return output_file
    
    builder = AuxFileBuilder()
    builder.build_for_species(
        species,
        source_dir=source_dir,
        output_file=output_file
    )
    
    return output_file


def build_all(force: bool = False) -> Dict[str, Path]:
    """
    Build aux files for all species with normalized data.
    
    Parameters
    ----------
    force : bool
        Force regeneration
        
    Returns
    -------
    Dict[str, Path]
        Mapping of species to aux file paths
    """
    results = {}
    species_list = get_normalized_species()
    
    if not species_list:
        logger.warning("No species with normalized data found")
        return results
    
    logger.info(f"Building aux files for {len(species_list)} species")
    
    for species in species_list:
        try:
            path = build_aux_file(species, force=force)
            results[species] = path
        except Exception as e:
            logger.error(f"Failed to build aux file for {species}: {e}")
    
    return results


def main():
    """Main entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    parser = argparse.ArgumentParser(
        description="Build IgBLAST auxiliary files from germline data"
    )
    parser.add_argument(
        "--species",
        nargs="+",
        help="Species to process (e.g., human mouse)"
    )
    parser.add_argument(
        "--all-species",
        action="store_true",
        help="Process all species with normalized data"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force regeneration even if files exist"
    )
    parser.add_argument(
        "--list-species",
        action="store_true",
        help="List species with available normalized data"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    if args.list_species:
        print("Species with normalized data:")
        for species in get_normalized_species():
            aux_exists = (AUX_OUTPUT_DIR / f"{species}_gl.aux").exists()
            status = "✓" if aux_exists else "✗"
            print(f"  {status} {species}")
        return
    
    if args.all_species:
        results = build_all(force=args.force)
        print(f"\nBuilt aux files for {len(results)} species")
        for species, path in results.items():
            print(f"  {species}: {path}")
    elif args.species:
        for species in args.species:
            try:
                path = build_aux_file(species, force=args.force)
                print(f"{species}: {path}")
            except Exception as e:
                logger.error(f"Failed for {species}: {e}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
