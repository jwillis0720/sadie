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
from typing import Dict, List, Optional, Set, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

logger = logging.getLogger(__name__)

CHECKPOINT_FILE = ".download_progress.json"
SPECIES_MAPPING_FILE = "species_mapping.json"


def _load_checkpoint(checkpoint_path: Path) -> dict:
    if checkpoint_path.exists():
        try:
            return json.loads(checkpoint_path.read_text())
        except Exception:
            return {}
    return {}


def _save_checkpoint(checkpoint_path: Path, data: dict):
    checkpoint_path.write_text(json.dumps(data, indent=2))


def _load_species_mapping(mapping_path: Path) -> Set[str]:
    """Load species mapping from JSON file."""
    if mapping_path.exists():
        try:
            data = json.loads(mapping_path.read_text())
            return set(data.get("species_variants", []))
        except Exception:
            return set()
    return set()


def _save_species_mapping(mapping_path: Path, species_variants: Set[str], internal_name: str):
    """Save species mapping to JSON file."""
    data = {
        "internal_name": internal_name,
        "species_variants": sorted(species_variants),
        "count": len(species_variants),
        "updated_at": datetime.now().isoformat(),
    }
    mapping_path.write_text(json.dumps(data, indent=2))


# IMGT V-QUEST reference directory base URL
IMGT_BASE_URL = "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory"

# IMGT GENE-DB URL for C genes (not available in V-QUEST)
IMGT_GENEDB_URL = "https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP"

# Species name mapping (internal names to IMGT directory names)
SPECIES_MAP = {
    "human": "Homo_sapiens",
    "mouse": "Mus_musculus",
    "mouse_c57bl6j": "Mus_musculus_C57BL6J",
    "rat": "Rattus_norvegicus",
    "rabbit": "Oryctolagus_cuniculus",
    "rhesus_macaque": "Macaca_mulatta",
    "macaque": "Macaca_mulatta",  # Alias for canonical name
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

    def __init__(self, output_dir: Optional[Path] = None, include_tr: bool = False, timeout: int = 30):
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

    def download(self, species: List[str], segments: Optional[List[str]] = None, force: bool = False) -> Dict[str, int]:
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

    def _download_species(self, imgt_name: str, internal_name: str, segments: List[str], force: bool) -> int:
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
        # Collect all species/strain names found in V-QUEST FASTA headers
        all_species_variants: Set[str] = set()

        for segment in segments:
            receptor_type = "IG" if segment.startswith("IG") else "TR"
            url = f"{IMGT_BASE_URL}/{imgt_name}/{receptor_type}/{segment}.fasta"

            gapped_path = species_dir / f"{segment}_gapped.fasta"
            ungapped_path = species_dir / f"{segment}.fasta"

            if segment in completed_segments and not force:
                if gapped_path.exists():
                    count = self._count_sequences(gapped_path)
                    # Extract species names from existing file
                    species_from_file = self._extract_species_from_fasta(gapped_path)
                    all_species_variants.update(species_from_file)
                    logger.debug(f"Skipping {segment} (checkpoint, {count} sequences)")
                    total_count += count
                    completed_count += 1
                    continue

            if gapped_path.exists() and ungapped_path.exists() and not force:
                count = self._count_sequences(gapped_path)
                # Extract species names from existing file
                species_from_file = self._extract_species_from_fasta(gapped_path)
                all_species_variants.update(species_from_file)
                logger.debug(f"Skipping {segment} (exists with {count} sequences)")
                total_count += count
                completed_segments.add(segment)
                completed_count += 1
                continue

            try:
                count, species_names = self._download_segment(url, gapped_path, ungapped_path)
                total_count += count
                all_species_variants.update(species_names)
                completed_segments.add(segment)
                completed_count += 1

                if count > 0:
                    logger.info(f"  {segment}: {count} sequences")

                if completed_count % 10 == 0 or completed_count == total_segments:
                    pct = int(100 * completed_count / total_segments)
                    logger.info(f"Downloaded {completed_count}/{total_segments} files ({pct}%)")

                _save_checkpoint(
                    checkpoint_path,
                    {
                        "completed_segments": list(completed_segments),
                        "total_segments": total_segments,
                        "timestamp": datetime.now().isoformat(),
                        "last_segment": segment,
                    },
                )

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

        # Save species mapping to JSON for future C gene downloads
        mapping_path = species_dir / SPECIES_MAPPING_FILE
        if all_species_variants:
            logger.info(f"Found {len(all_species_variants)} species/strain variants in V-QUEST data")
            _save_species_mapping(mapping_path, all_species_variants, internal_name)
            logger.info(f"Saved species mapping to {mapping_path}")

        # Download C genes from GENE-DB (not available in V-QUEST)
        # Use exact species names collected from V-QUEST headers (or load from JSON)
        c_count = self._download_c_genes_from_genedb(internal_name, all_species_variants, force)
        total_count += c_count

        return total_count

    def _extract_species_from_fasta(self, fasta_path: Path) -> Set[str]:
        """
        Extract unique species names from FASTA headers (field 3).

        Parameters
        ----------
        fasta_path : Path
            Path to FASTA file

        Returns
        -------
        Set[str]
            Set of species names found in headers
        """
        species_names: Set[str] = set()
        with open(fasta_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    parts = line[1:].split("|")
                    if len(parts) >= 3:
                        species_names.add(parts[2])
        return species_names

    def _download_c_genes_from_genedb(self, internal_name: str, species_variants: Set[str], force: bool = False) -> int:
        """
        Download C genes from IMGT GENE-DB.

        C genes are not available in V-QUEST reference directory, they must
        be fetched from the GENE-DB bulk FASTA file.

        Parameters
        ----------
        internal_name : str
            Internal species name (e.g., "human")
        species_variants : Set[str]
            Set of exact species/strain names to match (collected from V-QUEST headers)
        force : bool
            Force re-download

        Returns
        -------
        int
            Number of C gene sequences downloaded
        """
        species_dir = self.output_dir / internal_name

        # Check if C gene files already exist
        c_files_exist = all((species_dir / f"IG{chain}C.fasta").exists() for chain in ["H", "K", "L"])
        if c_files_exist and not force:
            # Count existing sequences
            count = sum(
                self._count_sequences(species_dir / f"IG{chain}C.fasta")
                for chain in ["H", "K", "L"]
                if (species_dir / f"IG{chain}C.fasta").exists()
            )
            logger.debug(f"Skipping C genes (exists with {count} sequences)")
            return count

        # If no species variants provided, try to load from saved mapping
        if not species_variants:
            mapping_path = species_dir / SPECIES_MAPPING_FILE
            species_variants = _load_species_mapping(mapping_path)
            if species_variants:
                logger.info(f"Loaded {len(species_variants)} species variants from {mapping_path}")

        if not species_variants:
            logger.debug(f"No species variants found for {internal_name}, skipping C genes")
            return 0

        logger.info(f"Downloading C genes from GENE-DB for {internal_name} ({len(species_variants)} variants)...")

        # Download and cache GENE-DB file
        genedb_content = self._get_genedb_content(force)
        if not genedb_content:
            logger.warning("Failed to get GENE-DB content")
            return 0

        # Parse and filter C genes for this species using exact matching
        c_genes = self._parse_genedb_c_genes(genedb_content, species_variants)

        if not c_genes or all(len(v) == 0 for v in c_genes.values()):
            logger.info(f"No C genes found in GENE-DB for {internal_name}")
            return 0

        # Group by chain and write FASTA files
        total_count = 0
        for chain in ["H", "K", "L"]:
            chain_genes = c_genes.get(chain, [])
            if chain_genes:
                fasta_path = species_dir / f"IG{chain}C.fasta"
                with open(fasta_path, "w") as f:
                    for gene_name, sequence in chain_genes:
                        f.write(f">{gene_name}\n{sequence}\n")
                logger.info(f"  IG{chain}C: {len(chain_genes)} sequences")
                total_count += len(chain_genes)

        return total_count

    def _get_genedb_content(self, force: bool = False) -> Optional[str]:
        """
        Get GENE-DB FASTA content, with caching.

        Parameters
        ----------
        force : bool
            Force re-download

        Returns
        -------
        str or None
            GENE-DB FASTA content
        """
        cache_dir = Path.home() / ".cache" / "sadie" / "imgt"
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / "GENEDB_C_genes.fasta"

        if cache_path.exists() and not force:
            logger.debug("Using cached GENE-DB file")
            return cache_path.read_text()

        try:
            logger.info(f"Downloading GENE-DB from {IMGT_GENEDB_URL}...")
            req = Request(IMGT_GENEDB_URL, headers={"User-Agent": "SADIE-Germlines/1.0"})
            with urlopen(req, timeout=self.timeout) as response:
                content = response.read().decode("utf-8")

            # Cache the content
            cache_path.write_text(content)
            logger.info(f"Cached GENE-DB to {cache_path}")
            return content

        except (URLError, HTTPError) as e:
            logger.error(f"Failed to download GENE-DB: {e}")
            return None

    def _parse_genedb_c_genes(self, content: str, species_variants: Set[str]) -> Dict[str, List[Tuple[str, str]]]:
        """
        Parse GENE-DB FASTA to extract C genes for a species.

        GENE-DB header format:
        >ACCESSION|GENE*ALLELE|SPECIES|FUNCTIONALITY|REGION|...

        Parameters
        ----------
        content : str
            GENE-DB FASTA content
        species_variants : Set[str]
            Set of exact species/strain names to match (from V-QUEST headers)

        Returns
        -------
        Dict[str, List[Tuple[str, str]]]
            C genes grouped by chain (H, K, L) as (gene_name, sequence) tuples
        """
        c_genes: Dict[str, List[Tuple[str, str]]] = {"H": [], "K": [], "L": []}

        # Extract base species prefixes for prefix matching
        # e.g., "Macaca mulatta_17573" -> "Macaca mulatta"
        # This allows matching GENE-DB strains not in V-QUEST (e.g., "Macaca mulatta_RUp15")
        species_prefixes: Set[str] = set()
        for variant in species_variants:
            # Split on underscore and take genus + species (first two words)
            # Handle cases like "Macaca mulatta" (no underscore) and "Macaca mulatta_17573"
            if "_" in variant:
                # Check if it's genus_species_strain format
                parts = variant.split("_")
                # If first part has space, it's already "Genus species"
                if " " in parts[0]:
                    species_prefixes.add(parts[0])
                else:
                    # It's "Genus_species_strain" format - reconstruct base
                    species_prefixes.add(f"{parts[0]} {parts[1]}" if len(parts) > 1 else parts[0])
            else:
                # No underscore - use as-is (e.g., "Homo sapiens")
                species_prefixes.add(variant)
        
        logger.debug(f"Using {len(species_prefixes)} species prefixes for C gene matching: {species_prefixes}")

        # C gene name patterns (including isotype genes)
        # IGHC: IGHA1, IGHA2, IGHD, IGHE, IGHG1-4, IGHM, IGHGP
        # IGKC: IGKC
        # IGLC: IGLC1-7
        c_gene_patterns = {
            "H": re.compile(r"^IGH[ADEGM]", re.IGNORECASE),  # IGHA, IGHD, IGHE, IGHG, IGHM
            "K": re.compile(r"^IGKC", re.IGNORECASE),  # IGKC
            "L": re.compile(r"^IGLC", re.IGNORECASE),  # IGLC1-7
        }

        current_header = None
        current_seq = []

        for line in content.split("\n"):
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save previous sequence
                if current_header:
                    gene_info = self._parse_genedb_header(current_header, species_prefixes, c_gene_patterns)
                    if gene_info:
                        chain, gene_name = gene_info
                        sequence = "".join(current_seq)
                        c_genes[chain].append((gene_name, sequence))

                current_header = line[1:]  # Remove ">"
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_header:
            gene_info = self._parse_genedb_header(current_header, species_prefixes, c_gene_patterns)
            if gene_info:
                chain, gene_name = gene_info
                sequence = "".join(current_seq)
                c_genes[chain].append((gene_name, sequence))

        return c_genes

    def _parse_genedb_header(
        self, header: str, species_prefixes: Set[str], c_gene_patterns: Dict[str, re.Pattern]
    ) -> Optional[Tuple[str, str]]:
        """
        Parse GENE-DB header to extract C gene info.

        Parameters
        ----------
        header : str
            FASTA header (without ">")
        species_prefixes : Set[str]
            Set of species name prefixes to match (e.g., "Macaca mulatta")
        c_gene_patterns : Dict
            Regex patterns for C gene identification

        Returns
        -------
        Tuple[str, str] or None
            (chain, gene_name) or None if not a C gene for this species
        """
        parts = header.split("|")
        if len(parts) < 4:
            return None

        # GENE-DB format: ACCESSION|GENE*ALLELE|SPECIES|FUNCTIONALITY|...
        gene_name = parts[1]
        species = parts[2]
        functionality = parts[3] if len(parts) > 3 else "F"
        # Strip parentheses and brackets from functionality (e.g., "(F)" -> "F", "[ORF]" -> "ORF")
        functionality = functionality.strip("()[]")

        # Filter by species using prefix matching
        # This captures all strains like "Macaca mulatta_RUp15" when prefix is "Macaca mulatta"
        species_match = any(species.startswith(prefix) for prefix in species_prefixes)
        if not species_match:
            return None

        # Filter by functionality (F = functional, ORF = open reading frame)
        # Include F and ORF, exclude P (pseudogene)
        if functionality not in ["F", "ORF"]:
            return None

        # Check if this is a C gene
        for chain, pattern in c_gene_patterns.items():
            if pattern.match(gene_name):
                return (chain, gene_name)

        return None

    def _download_segment(self, url: str, gapped_path: Path, ungapped_path: Path) -> Tuple[int, Set[str]]:
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
        Tuple[int, Set[str]]
            Number of sequences downloaded and set of species names found in headers
        """
        # Download file
        req = Request(url, headers={"User-Agent": "SADIE-Germlines/1.0"})
        with urlopen(req, timeout=self.timeout) as response:
            content = response.read().decode("utf-8")

        if not content.strip():
            return 0, set()

        # Parse and process sequences
        sequences = self._parse_imgt_fasta(content)

        if not sequences:
            return 0, set()

        # Extract species names from headers (field 3, pipe-delimited)
        species_names: Set[str] = set()
        for header, _ in sequences:
            parts = header.split("|")
            if len(parts) >= 3:
                species_names.add(parts[2])

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

        return len(sequences), species_names

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
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    parser = argparse.ArgumentParser(description="Download IMGT germline data from V-QUEST reference directory")
    parser.add_argument(
        "--species", nargs="+", default=["human"], help="Species to download (e.g., human mouse rabbit)"
    )
    parser.add_argument("--segments", nargs="+", default=None, help="Specific segments to download (e.g., IGHV IGHJ)")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output directory for FASTA files")
    parser.add_argument("--include-tr", action="store_true", help="Include T-cell receptor (TR) sequences")
    parser.add_argument("--list-species", action="store_true", help="List all available species and exit")
    parser.add_argument("--force", action="store_true", help="Force re-download even if files exist")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    downloader = IMGTDownloader(output_dir=args.output_dir, include_tr=args.include_tr)

    try:
        if args.list_species:
            print("Available IMGT species:")
            print("-" * 40)
            for species in downloader.list_available_species():
                internal = SPECIES_MAP_REVERSE.get(species, species.lower())
                print(f"  {internal:20} ({species})")
            return

        results = downloader.download(args.species, segments=args.segments, force=args.force)

        print(f"\nIMGT data downloaded successfully to {downloader.output_dir}")
        print("\nSummary:")
        for species, count in results.items():
            print(f"  {species}: {count} sequences")

    except Exception as e:
        logger.error(f"Download failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
