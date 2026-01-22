#!/usr/bin/env python3
"""
Convert G3 Custom JSON to FASTA Files
=====================================

Script to convert G3-format custom germline JSON files to FASTA files
organized by species and segment.

Data Source: G3.custom.json (exported from G3 database)

The G3 JSON format contains:
- source: "custom"
- common: common species name (e.g., "cat", "macaque")
- latin: latin species name
- gene: gene name (e.g., "IGHV1-1*01")
- gene_segment: segment type (V, D, J)
- sequence: ungapped nucleotide sequence
- imgt.sequence_gapped: IMGT-gapped sequence (for validation)

Usage:
    python download_custom.py --input G3.custom.json
    python download_custom.py --input G3.custom.json --output-dir ./custom
    python download_custom.py --input G3.custom.json --list-species

Output Directory Structure:
    sources/custom/
    ├── cat/
    │   ├── IGHV.fasta
    │   └── ...
    ├── crab_eating_macaque/
    │   ├── IGHV.fasta
    │   └── ...
    └── macaque/
        └── ...

Note: This script outputs ONLY ungapped sequences. Gapped sequences are
programmatically created by the GapperService using IMGT-gapped templates.
"""

import argparse
import json
import logging
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Any

logger = logging.getLogger(__name__)


class CustomConverter:
    """Convert G3 custom JSON to FASTA files."""
    
    def __init__(
        self,
        output_dir: Optional[Path] = None,
    ):
        """
        Initialize custom converter.
        
        Parameters
        ----------
        output_dir : Path, optional
            Output directory for FASTA files.
            Defaults to sources/custom/
        """
        if output_dir is None:
            output_dir = Path(__file__).parent.parent / "sources" / "custom"
        
        self.output_dir = Path(output_dir)
    
    def _normalize_species_name(self, common: str) -> str:
        """
        Normalize species common name for directory naming.
        
        Parameters
        ----------
        common : str
            Common species name from JSON (e.g., "crab_eating_macaque")
            
        Returns
        -------
        str
            Normalized name suitable for directory
        """
        # Replace spaces with underscores, lowercase
        normalized = common.lower().replace(" ", "_").replace("-", "_")
        # Remove any non-alphanumeric chars except underscore
        normalized = re.sub(r'[^a-z0-9_]', '', normalized)
        return normalized
    
    def _parse_gene_info(self, gene_name: str) -> Dict[str, str]:
        """
        Parse gene name to extract chain and segment.
        
        Parameters
        ----------
        gene_name : str
            Gene name (e.g., "IGHV1-69*01")
            
        Returns
        -------
        Dict[str, str]
            Dict with 'chain' and 'segment' keys
        """
        upper_name = gene_name.upper()
        
        # Try standard IMGT format: IG{chain}{segment}...
        match = re.match(r'IG([HKL])([VDJ])', upper_name)
        if match:
            return {
                'chain': match.group(1),
                'segment': match.group(2)
            }
        
        # Handle custom kappa/lambda naming patterns like:
        # IGK-lib4kappa_1 -> IGKV
        # IGL-F124_lambda_10 -> IGLV
        if upper_name.startswith('IGK') and 'KAPPA' in upper_name:
            return {'chain': 'K', 'segment': 'V'}
        if upper_name.startswith('IGL') and 'LAMBDA' in upper_name:
            return {'chain': 'L', 'segment': 'V'}
        
        # Handle TR (T-cell receptor) - skip or default
        if upper_name.startswith('TR'):
            logger.debug(f"Skipping T-cell receptor gene: {gene_name}")
            return {'chain': 'H', 'segment': 'V'}  # Will be filtered if TR support not needed
        
        # Default to heavy chain V if can't parse
        logger.warning(f"Could not parse gene name: {gene_name}, defaulting to IGHV")
        return {'chain': 'H', 'segment': 'V'}
    
    def load_json(self, json_path: Path) -> List[Dict[str, Any]]:
        """
        Load G3 custom JSON file.
        
        Parameters
        ----------
        json_path : Path
            Path to JSON file
            
        Returns
        -------
        List[Dict[str, Any]]
            List of sequence entries
        """
        logger.info(f"Loading JSON from {json_path}")
        
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        if isinstance(data, list):
            return data
        elif isinstance(data, dict) and 'sequences' in data:
            return data['sequences']
        else:
            raise ValueError("Unexpected JSON structure")
    
    def list_species(self, entries: List[Dict[str, Any]]) -> Dict[str, int]:
        """
        List all species in the JSON data.
        
        Parameters
        ----------
        entries : List[Dict[str, Any]]
            Sequence entries
            
        Returns
        -------
        Dict[str, int]
            Species to count mapping
        """
        species_count = defaultdict(int)
        
        for entry in entries:
            common = entry.get('common', 'unknown')
            species_count[common] += 1
        
        return dict(sorted(species_count.items()))
    
    def convert(
        self,
        json_path: Path,
        species_filter: Optional[List[str]] = None,
        include_gapped: bool = False,
        force: bool = False
    ) -> Dict[str, int]:
        """
        Convert G3 JSON to FASTA files.
        
        Parameters
        ----------
        json_path : Path
            Path to input JSON file
        species_filter : List[str], optional
            Only process these species
        include_gapped : bool
            Also write gapped FASTA files (for validation)
        force : bool
            Overwrite existing files
            
        Returns
        -------
        Dict[str, int]
            Dictionary mapping species to sequence count
        """
        start_time = time.time()
        
        # Load JSON
        entries = self.load_json(json_path)
        logger.info(f"Loaded {len(entries)} entries from JSON")
        
        # Organize sequences by species/chain/segment
        organized: Dict[str, Dict[str, List[Dict]]] = defaultdict(
            lambda: defaultdict(list)
        )
        
        for entry in entries:
            # Extract species
            common = entry.get('common', 'unknown')
            species = self._normalize_species_name(common)
            
            # Filter species if requested
            if species_filter:
                if species not in species_filter and common not in species_filter:
                    continue
            
            # Extract gene info
            gene_name = entry.get('gene', 'UNKNOWN')
            gene_info = self._parse_gene_info(gene_name)
            
            # Use gene_segment from JSON if available
            segment = entry.get('gene_segment', gene_info['segment'])
            chain = gene_info['chain']
            
            # Get sequences
            ungapped = entry.get('sequence', '')
            gapped = entry.get('imgt', {}).get('sequence_gapped', '')
            
            if not ungapped:
                logger.warning(f"No sequence for {gene_name}, skipping")
                continue
            
            # Key for file grouping
            file_key = f"IG{chain}{segment}"
            
            organized[species][file_key].append({
                'gene_name': gene_name,
                'ungapped': ungapped.upper(),
                'gapped': gapped.upper() if gapped else None,
                'common': common,
                'latin': entry.get('latin', ''),
            })
        
        # Write FASTA files
        results = {}
        
        for species, files in organized.items():
            species_dir = self.output_dir / species
            species_dir.mkdir(parents=True, exist_ok=True)
            
            species_count = 0
            
            for file_key, sequences in files.items():
                # Write ungapped FASTA
                ungapped_path = species_dir / f"{file_key}.fasta"
                
                if ungapped_path.exists() and not force:
                    logger.info(f"Skipping {ungapped_path} (exists)")
                    species_count += len(sequences)
                    continue
                
                with open(ungapped_path, 'w') as f:
                    for seq in sequences:
                        f.write(f">{seq['gene_name']}\n")
                        f.write(f"{seq['ungapped']}\n")
                
                logger.info(f"  {file_key}: {len(sequences)} sequences")
                species_count += len(sequences)
                
                # Optionally write gapped FASTA (for validation/comparison)
                if include_gapped:
                    gapped_path = species_dir / f"{file_key}_gapped.fasta"
                    with open(gapped_path, 'w') as f:
                        for seq in sequences:
                            if seq['gapped']:
                                f.write(f">{seq['gene_name']}\n")
                                f.write(f"{seq['gapped']}\n")
                    logger.info(f"  {file_key}_gapped: wrote validation file")
            
            results[species] = species_count
            logger.info(f"Converted {species_count} sequences for {species}")
        
        duration_ms = int((time.time() - start_time) * 1000)
        total_seqs = sum(results.values())
        logger.info(
            f"operation=convert provider=custom "
            f"species={len(results)} sequences={total_seqs} "
            f"duration_ms={duration_ms} status=success"
        )
        
        return results


def main():
    """Main entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    parser = argparse.ArgumentParser(
        description="Convert G3 custom JSON to FASTA files"
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input JSON file (G3.custom.json)"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=None,
        help="Output directory for FASTA files"
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=None,
        help="Only process these species"
    )
    parser.add_argument(
        "--list-species",
        action="store_true",
        help="List all species in JSON and exit"
    )
    parser.add_argument(
        "--include-gapped",
        action="store_true",
        help="Also write gapped FASTA files (for validation)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    converter = CustomConverter(output_dir=args.output_dir)
    
    try:
        if args.list_species:
            entries = converter.load_json(args.input)
            species_count = converter.list_species(entries)
            
            print(f"\nSpecies in {args.input}:")
            print("-" * 50)
            for species, count in species_count.items():
                normalized = converter._normalize_species_name(species)
                print(f"  {species:30} ({normalized}): {count}")
            print("-" * 50)
            print(f"Total: {len(species_count)} species, {sum(species_count.values())} sequences")
            return
        
        results = converter.convert(
            args.input,
            species_filter=args.species,
            include_gapped=args.include_gapped,
            force=args.force
        )
        
        print(f"\nCustom data converted successfully to {converter.output_dir}")
        print("\nSummary:")
        for species, count in results.items():
            print(f"  {species}: {count} sequences")
        
    except Exception as e:
        logger.error(f"Conversion failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
