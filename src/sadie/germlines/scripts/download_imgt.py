#!/usr/bin/env python3
"""
Download IMGT Germline Data
============================

Script to download IMGT germline sequences from the V-QUEST reference directory.

Data Source: https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/

The IMGT reference directory contains:
- FASTA files with IMGT-gapped sequences (dots indicate gaps)
- Files organized by species and chain type (IG/TR)
- Header format: >accession|gene_name|species|functionality|region|positions|...

Usage:
    python download_imgt.py --species human
    python download_imgt.py --species human mouse --output-dir ./data
    python download_imgt.py --list-species
    python download_imgt.py --help

Output Directory Structure:
    sources/imgt/
    ├── human/
    │   ├── IGHV.fasta          # Ungapped sequences (dots removed)
    │   ├── IGHV_gapped.fasta   # IMGT-gapped sequences (original)
    │   ├── IGHD.fasta
    │   ├── IGHD_gapped.fasta
    │   ├── IGHJ.fasta
    │   ├── IGHJ_gapped.fasta
    │   └── ...
    └── mouse/
        └── ...
"""

import argparse
import json
import logging
import re
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

logger = logging.getLogger(__name__)

CHECKPOINT_FILE = ".download_progress.json"


def _load_checkpoint(checkpoint_path: Path) -> dict:
    if checkpoint_path.exists():
        try:
            return json.loads(checkpoint_path.read_text())
        except Exception:
            return {}
    return {}


def _save_checkpoint(checkpoint_path: Path, data: dict):
    checkpoint_path.write_text(json.dumps(data, indent=2))

# IMGT V-QUEST reference directory base URL
IMGT_BASE_URL = "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory"

# Species name mapping (internal names to IMGT directory names)
SPECIES_MAP = {
    "human": "Homo_sapiens",
    "mouse": "Mus_musculus",
    "mouse_c57bl6j": "Mus_musculus_C57BL6J",
    "rat": "Rattus_norvegicus",
    "rabbit": "Oryctolagus_cuniculus",
    "rhesus_macaque": "Macaca_mulatta",
    "cynomolgus": "Macaca_fascicularis",
    "dog": "Canis_lupus_familiaris",
    "cat": "Felis_catus",
    "pig": "Sus_scrofa",
    "cow": "Bos_taurus",
    "sheep": "Ovis_aries",
    "goat": "Capra_hircus",
    "horse": "Equus_caballus",
    "chicken": "Gallus_gallus",
    "alpaca": "Vicugna_pacos",
    "camel": "Camelus_dromedarius",
    "zebrafish": "Danio_rerio",
    "gorilla": "Gorilla_gorilla_gorilla",
    "chimpanzee": "Pan_troglodytes",
    "orangutan_sumatran": "Pongo_abelii",
    "orangutan_bornean": "Pongo_pygmaeus",
    "ferret": "Mustela_putorius_furo",
    "mink": "Neogale_vison",
    "dolphin": "Tursiops_truncatus",
    "platypus": "Ornithorhynchus_anatinus",
    "atlantic_salmon": "Salmo_salar",
    "rainbow_trout": "Oncorhynchus_mykiss",
    "atlantic_cod": "Gadus_morhua",
    "channel_catfish": "Ictalurus_punctatus",
    "lemur": "Lemur_catta",
    "owl_monkey": "Aotus_nancymaae",
    "naked_mole_rat": "Heterocephalus_glaber",
}

SPECIES_MAP_REVERSE = {v: k for k, v in SPECIES_MAP.items()}

# Standard IG segment files
IG_SEGMENTS = ["IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ"]

# TR segment files (T-cell receptors)
TR_SEGMENTS = ["TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ", "TRDV", "TRDD", "TRDJ", "TRGV", "TRGJ"]


class IMGTDownloader:
    """Download and process IMGT reference sequences."""
    
    def __init__(
        self,
        output_dir: Optional[Path] = None,
        include_tr: bool = False,
        timeout: int = 30
    ):
        """
        Initialize IMGT downloader.
        
        Parameters
        ----------
        output_dir : Path, optional
            Output directory for FASTA files.
            Defaults to sources/imgt/
        include_tr : bool
            Include T-cell receptor (TR) sequences
        timeout : int
            HTTP request timeout in seconds
        """
        if output_dir is None:
            output_dir = Path(__file__).parent.parent / "sources" / "imgt"
        
        self.output_dir = Path(output_dir)
        self.include_tr = include_tr
        self.timeout = timeout
    
    def list_available_species(self) -> List[str]:
        """
        List all available species from IMGT.
        
        Returns
        -------
        List[str]
            List of species directory names
        """
        url = f"{IMGT_BASE_URL}/"
        logger.info(f"Fetching species list from {url}")
        
        try:
            req = Request(url, headers={"User-Agent": "SADIE-Germlines/1.0"})
            with urlopen(req, timeout=self.timeout) as response:
                html = response.read().decode("utf-8")
        except (URLError, HTTPError) as e:
            raise RuntimeError(f"Failed to fetch species list: {e}")
        
        # Parse directory listing for species folders
        # Looking for: href="Species_name/"
        species_pattern = re.compile(r'href="([A-Z][a-z_]+(?:_[a-z]+)*)/?"')
        species = []
        
        for match in species_pattern.finditer(html):
            name = match.group(1)
            # Filter out non-species directories
            if name not in ("icons", "images", "css"):
                species.append(name)
        
        return sorted(species)
    
    def download(
        self,
        species: List[str],
        segments: Optional[List[str]] = None,
        force: bool = False
    ) -> Dict[str, int]:
        """
        Download IMGT data for specified species.
        
        Parameters
        ----------
        species : List[str]
            Species to download (internal names or IMGT names)
        segments : List[str], optional
            Specific segments to download (e.g., ["IGHV", "IGHJ"])
            Defaults to all IG segments (and TR if include_tr=True)
        force : bool
            Force re-download even if files exist
            
        Returns
        -------
        Dict[str, int]
            Dictionary mapping species to sequence count
        """
        start_time = time.time()
        results = {}
        
        # Determine segments to download
        if segments is None:
            segments = IG_SEGMENTS.copy()
            if self.include_tr:
                segments.extend(TR_SEGMENTS)
        
        for sp in species:
            # Map internal name to IMGT name if needed
            imgt_name = SPECIES_MAP.get(sp.lower(), sp)
            internal_name = SPECIES_MAP_REVERSE.get(imgt_name, sp.lower().replace(" ", "_"))
            
            logger.info(f"Downloading IMGT data for {internal_name} ({imgt_name})...")
            
            try:
                count = self._download_species(imgt_name, internal_name, segments, force)
                results[internal_name] = count
                logger.info(f"Downloaded {count} sequences for {internal_name}")
            except Exception as e:
                logger.error(f"Failed to download {internal_name}: {e}")
                results[internal_name] = 0
        
        duration_ms = int((time.time() - start_time) * 1000)
        total_seqs = sum(results.values())
        logger.info(
            f"operation=download provider=imgt "
            f"species={','.join(species)} sequences={total_seqs} "
            f"duration_ms={duration_ms} status=success"
        )
        
        return results
    
    def _download_species(
        self,
        imgt_name: str,
        internal_name: str,
        segments: List[str],
        force: bool
    ) -> int:
        """
        Download all segments for a species.

        Parameters
        ----------
        imgt_name : str
            IMGT species directory name (e.g., "Homo_sapiens")
        internal_name : str
            Internal species name (e.g., "human")
        segments : List[str]
            Segments to download
        force : bool
            Force re-download

        Returns
        -------
        int
            Total number of sequences downloaded
        """
        species_dir = self.output_dir / internal_name
        species_dir.mkdir(parents=True, exist_ok=True)

        checkpoint_path = species_dir / CHECKPOINT_FILE
        checkpoint = _load_checkpoint(checkpoint_path)
        completed_segments = set(checkpoint.get("completed_segments", []))

        total_count = 0
        total_segments = len(segments)
        completed_count = 0

        for segment in segments:
            receptor_type = "IG" if segment.startswith("IG") else "TR"
            url = f"{IMGT_BASE_URL}/{imgt_name}/{receptor_type}/{segment}.fasta"

            gapped_path = species_dir / f"{segment}_gapped.fasta"
            ungapped_path = species_dir / f"{segment}.fasta"

            if segment in completed_segments and not force:
                if gapped_path.exists():
                    count = self._count_sequences(gapped_path)
                    logger.debug(f"Skipping {segment} (checkpoint, {count} sequences)")
                    total_count += count
                    completed_count += 1
                    continue

            if gapped_path.exists() and ungapped_path.exists() and not force:
                count = self._count_sequences(gapped_path)
                logger.debug(f"Skipping {segment} (exists with {count} sequences)")
                total_count += count
                completed_segments.add(segment)
                completed_count += 1
                continue

            try:
                count = self._download_segment(url, gapped_path, ungapped_path)
                total_count += count
                completed_segments.add(segment)
                completed_count += 1

                if count > 0:
                    logger.info(f"  {segment}: {count} sequences")

                if completed_count % 10 == 0 or completed_count == total_segments:
                    pct = int(100 * completed_count / total_segments)
                    logger.info(f"Downloaded {completed_count}/{total_segments} files ({pct}%)")

                _save_checkpoint(checkpoint_path, {
                    "completed_segments": list(completed_segments),
                    "total_segments": total_segments,
                    "timestamp": datetime.now().isoformat(),
                    "last_segment": segment
                })

            except HTTPError as e:
                if e.code == 404:
                    logger.debug(f"  {segment}: not available for this species")
                    completed_segments.add(segment)
                else:
                    logger.warning(f"  {segment}: HTTP error {e.code}")
            except Exception as e:
                logger.warning(f"  {segment}: failed - {e}")

        if checkpoint_path.exists():
            checkpoint_path.unlink()

        return total_count
    
    def _download_segment(
        self,
        url: str,
        gapped_path: Path,
        ungapped_path: Path
    ) -> int:
        """
        Download a single segment file.
        
        Parameters
        ----------
        url : str
            URL to download from
        gapped_path : Path
            Output path for gapped FASTA
        ungapped_path : Path
            Output path for ungapped FASTA
            
        Returns
        -------
        int
            Number of sequences downloaded
        """
        # Download file
        req = Request(url, headers={"User-Agent": "SADIE-Germlines/1.0"})
        with urlopen(req, timeout=self.timeout) as response:
            content = response.read().decode("utf-8")
        
        if not content.strip():
            return 0
        
        # Parse and process sequences
        sequences = self._parse_imgt_fasta(content)
        
        if not sequences:
            return 0
        
        # Write gapped FASTA (original IMGT format)
        with open(gapped_path, "w") as f:
            for header, seq in sequences:
                f.write(f">{header}\n")
                # Write sequence in original multiline format
                f.write(f"{seq}\n")
        
        # Write ungapped FASTA (dots removed)
        with open(ungapped_path, "w") as f:
            for header, seq in sequences:
                # Remove gaps (dots) and join lines
                ungapped = seq.replace(".", "").replace("\n", "")
                f.write(f">{header}\n")
                f.write(f"{ungapped}\n")
        
        return len(sequences)
    
    def _parse_imgt_fasta(self, content: str) -> List[Tuple[str, str]]:
        """
        Parse IMGT FASTA content.
        
        IMGT FASTA format:
        >accession|gene_name|species|functionality|region|positions|length|...
        cag.gtgcagctggtgcag...tctggggctgag...gtgaag...
        
        Parameters
        ----------
        content : str
            Raw FASTA content
            
        Returns
        -------
        List[Tuple[str, str]]
            List of (header, sequence) tuples
        """
        sequences = []
        current_header = None
        current_seq = []
        
        for line in content.split("\n"):
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                # Save previous sequence
                if current_header and current_seq:
                    sequences.append((current_header, "".join(current_seq)))
                
                # Start new sequence
                current_header = line[1:]  # Remove ">"
                current_seq = []
            else:
                # Sequence line - keep as-is (includes dots for gaps)
                current_seq.append(line)
        
        # Save last sequence
        if current_header and current_seq:
            sequences.append((current_header, "".join(current_seq)))
        
        return sequences
    
    def _count_sequences(self, fasta_path: Path) -> int:
        """
        Count sequences in a FASTA file.
        
        Parameters
        ----------
        fasta_path : Path
            Path to FASTA file
            
        Returns
        -------
        int
            Number of sequences
        """
        count = 0
        with open(fasta_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
        return count


def main():
    """Main entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    parser = argparse.ArgumentParser(
        description="Download IMGT germline data from V-QUEST reference directory"
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=["human"],
        help="Species to download (e.g., human mouse rabbit)"
    )
    parser.add_argument(
        "--segments",
        nargs="+",
        default=None,
        help="Specific segments to download (e.g., IGHV IGHJ)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for FASTA files"
    )
    parser.add_argument(
        "--include-tr",
        action="store_true",
        help="Include T-cell receptor (TR) sequences"
    )
    parser.add_argument(
        "--list-species",
        action="store_true",
        help="List all available species and exit"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if files exist"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    downloader = IMGTDownloader(
        output_dir=args.output_dir,
        include_tr=args.include_tr
    )
    
    try:
        if args.list_species:
            print("Available IMGT species:")
            print("-" * 40)
            for species in downloader.list_available_species():
                internal = SPECIES_MAP_REVERSE.get(species, species.lower())
                print(f"  {internal:20} ({species})")
            return
        
        results = downloader.download(
            args.species,
            segments=args.segments,
            force=args.force
        )
        
        print(f"\nIMGT data downloaded successfully to {downloader.output_dir}")
        print("\nSummary:")
        for species, count in results.items():
            print(f"  {species}: {count} sequences")
        
    except Exception as e:
        logger.error(f"Download failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
