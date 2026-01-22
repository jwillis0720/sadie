#!/usr/bin/env python3
"""
Update organism.yaml Configuration
===================================

Script to update organism.yaml with species configurations for IgBLAST.

Scans the BLAST database directory and aux_db directory to find available
species and generates proper organism.yaml entries.

Usage:
    python -m sadie.germlines.scripts.update_organism_yaml --all-species
    python -m sadie.germlines.scripts.update_organism_yaml --species human mouse
    python -m sadie.germlines.scripts.update_organism_yaml --list-species
    python -m sadie.germlines.scripts.update_organism_yaml --help

Output:
    Updates src/sadie/germlines/igblast/internal_data/organism.yaml
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional
import yaml

logger = logging.getLogger(__name__)

# Module paths
GERMLINES_ROOT = Path(__file__).parent.parent
IGBLAST_DIR = GERMLINES_ROOT / "igblast"
DATABASE_DIR = IGBLAST_DIR / "database"
AUX_DB_DIR = IGBLAST_DIR / "aux_db"
ORGANISM_YAML = IGBLAST_DIR / "internal_data" / "organism.yaml"

# Standard segments
SEGMENTS = ["V", "D", "J"]


class OrganismYamlUpdater:
    """Update organism.yaml with species configurations."""
    
    def __init__(
        self,
        database_dir: Optional[Path] = None,
        aux_db_dir: Optional[Path] = None,
        organism_yaml: Optional[Path] = None
    ):
        """
        Initialize updater.
        
        Parameters
        ----------
        database_dir : Path, optional
            Directory containing BLAST databases
        aux_db_dir : Path, optional
            Directory containing aux files
        organism_yaml : Path, optional
            Path to organism.yaml file
        """
        self.database_dir = database_dir or DATABASE_DIR
        self.aux_db_dir = aux_db_dir or AUX_DB_DIR
        self.organism_yaml = organism_yaml or ORGANISM_YAML
    
    def get_available_species(self) -> Dict[str, Dict]:
        """
        Get available species with their database status.
        
        Returns
        -------
        Dict[str, Dict]
            Species info with database and aux file status
        """
        species_info = {}
        
        if not self.database_dir.exists():
            return species_info
        
        for species_dir in self.database_dir.iterdir():
            if not species_dir.is_dir() or species_dir.name.startswith('.'):
                continue
            
            species = species_dir.name
            
            # Check which segments exist
            segments = []
            for segment in SEGMENTS:
                db_path = species_dir / f"{species}_{segment}.nsq"
                if db_path.exists():
                    segments.append(segment)
            
            if not segments:
                continue
            
            # Check if aux file exists
            aux_file = self.aux_db_dir / f"{species}_gl.aux"
            
            species_info[species] = {
                "segments": segments,
                "has_aux": aux_file.exists(),
                "database_path": f"../database/{species}",
                "aux_file": f"../aux_db/{species}_gl.aux" if aux_file.exists() else None
            }
        
        return species_info
    
    def load_existing_config(self) -> Dict:
        """
        Load existing organism.yaml configuration.
        
        Returns
        -------
        Dict
            Existing configuration
        """
        if self.organism_yaml.exists():
            with open(self.organism_yaml, "r") as f:
                return yaml.safe_load(f) or {"organisms": {}}
        return {"organisms": {}}
    
    def update(
        self,
        species: Optional[List[str]] = None,
        force: bool = False
    ) -> Dict:
        """
        Update organism.yaml with species configurations.
        
        Parameters
        ----------
        species : List[str], optional
            Specific species to update. If None, updates all available.
        force : bool
            Force overwrite existing entries
            
        Returns
        -------
        Dict
            Updated configuration
        """
        available = self.get_available_species()
        config = self.load_existing_config()
        
        if "organisms" not in config:
            config["organisms"] = {}
        
        # Filter to requested species if specified
        if species:
            species_to_update = {s: available[s] for s in species if s in available}
            missing = set(species) - set(available.keys())
            if missing:
                logger.warning(f"Species not found in databases: {missing}")
        else:
            species_to_update = available
        
        updated_count = 0
        for sp, info in species_to_update.items():
            if sp in config["organisms"] and not force:
                logger.debug(f"Skipping existing: {sp}")
                continue
            
            # Create species entry
            entry = {
                "database_path": info["database_path"],
                "segments": info["segments"]
            }
            
            if info["has_aux"]:
                entry["aux_file"] = info["aux_file"]
            
            config["organisms"][sp] = entry
            updated_count += 1
            logger.info(f"Added/updated: {sp}")
        
        # Sort organisms alphabetically
        config["organisms"] = dict(sorted(config["organisms"].items()))
        
        # Write updated config
        self.organism_yaml.parent.mkdir(parents=True, exist_ok=True)
        with open(self.organism_yaml, "w") as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        
        logger.info(f"Updated {updated_count} species in {self.organism_yaml}")
        return config
    
    def validate(self) -> Dict[str, List[str]]:
        """
        Validate organism.yaml configuration.
        
        Returns
        -------
        Dict[str, List[str]]
            Validation results with 'errors' and 'warnings' lists
        """
        results = {"errors": [], "warnings": []}
        
        config = self.load_existing_config()
        
        for species, entry in config.get("organisms", {}).items():
            # Check database path
            db_path = self.database_dir / species
            if not db_path.exists():
                results["errors"].append(f"{species}: database directory not found")
                continue
            
            # Check segments
            for segment in entry.get("segments", []):
                db_file = db_path / f"{species}_{segment}.nsq"
                if not db_file.exists():
                    results["errors"].append(f"{species}: {segment} database not found")
            
            # Check aux file
            if "aux_file" in entry:
                aux_path = self.aux_db_dir / f"{species}_gl.aux"
                if not aux_path.exists():
                    results["warnings"].append(f"{species}: aux file configured but not found")
        
        return results


def main():
    """Main entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    parser = argparse.ArgumentParser(
        description="Update organism.yaml with species configurations"
    )
    parser.add_argument(
        "--species",
        nargs="+",
        help="Species to add/update"
    )
    parser.add_argument(
        "--all-species",
        action="store_true",
        help="Update all species with available databases"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force overwrite existing entries"
    )
    parser.add_argument(
        "--list-species",
        action="store_true",
        help="List available species and their status"
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate existing configuration"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    updater = OrganismYamlUpdater()
    
    if args.list_species:
        print("Available species:")
        print("-" * 60)
        available = updater.get_available_species()
        config = updater.load_existing_config()
        configured = set(config.get("organisms", {}).keys())
        
        for species, info in sorted(available.items()):
            in_yaml = "✓" if species in configured else "✗"
            has_aux = "aux" if info["has_aux"] else "   "
            segments = ",".join(info["segments"])
            print(f"  {in_yaml} {species:25} [{segments:5}] {has_aux}")
        
        print(f"\nTotal: {len(available)} species available")
        return
    
    if args.validate:
        results = updater.validate()
        if results["errors"]:
            print("Errors:")
            for err in results["errors"]:
                print(f"  ✗ {err}")
        if results["warnings"]:
            print("Warnings:")
            for warn in results["warnings"]:
                print(f"  ! {warn}")
        if not results["errors"] and not results["warnings"]:
            print("✓ Configuration is valid")
        return
    
    if args.all_species:
        config = updater.update(force=args.force)
        print(f"\nUpdated organism.yaml with {len(config['organisms'])} species")
    elif args.species:
        config = updater.update(species=args.species, force=args.force)
        print(f"\nProcessed species: {args.species}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
