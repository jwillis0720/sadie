"""
Germline Pipeline - Staged Data Processing
===========================================

Simple 3-stage pipeline:
1. sources/ → detect changes
2. sources/ → normalize → normalized/
3. normalized/ → build → igblast/

Design Principles:
- Simple file timestamp comparison for change detection
- Pre-build everything (no on-demand transformations)
- Each stage is independent and testable
"""

import hashlib
import logging
import time
from pathlib import Path
from typing import List
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .models import GermlineGene


logger = logging.getLogger(__name__)


def _log_timing(operation: str, start_time: float, **kwargs):
    duration_ms = (time.time() - start_time) * 1000
    extra = " ".join(f"{k}={v}" for k, v in kwargs.items())
    logger.info(f"operation={operation} duration_ms={duration_ms:.1f} {extra}")


def _file_hash(path: Path) -> str:
    h = hashlib.md5()
    h.update(path.read_bytes())
    return h.hexdigest()[:8]


def _check_offline_ready(sources_dir: Path, species: str) -> tuple[bool, str]:
    species_found = False
    for provider in ["imgt", "ogrdb", "vdjbase", "custom"]:
        provider_species_dir = sources_dir / provider / species
        if provider_species_dir.exists():
            fasta_files = list(provider_species_dir.glob("*.fasta"))
            if fasta_files:
                species_found = True
                break

    if not species_found:
        return False, (
            f"No germline data found for '{species}'. "
            f"Run download scripts first:\n"
            f"  python src/sadie/germlines/scripts/download_imgt.py --species {species}\n"
            f"  python src/sadie/germlines/scripts/download_ogrdb.py --species {species}\n"
            f"Or add custom sequences to: sources/custom/{species}/"
        )
    return True, ""


# Constants (extracted magic numbers per Zen)
SEGMENTS = ["V", "D", "J"]
CHAINS = ["H", "K", "L"]


class GermlinePipeline:
    """
    Simple pipeline using file timestamps for change detection.

    Pipeline Stages:
    1. Detect changes in sources/
    2. If changed → rebuild normalized/ (gap, merge, deduplicate)
    3. If normalized/ changed → rebuild igblast/ (BLAST DBs, aux files)

    Examples
    --------
    >>> from pathlib import Path
    >>> pipeline = GermlinePipeline(Path("src/sadie/germlines"))
    >>> pipeline.update("human")  # Auto-detects and rebuilds as needed
    """

    def __init__(self, base_dir: Path):
        """
        Initialize pipeline.

        Parameters
        ----------
        base_dir : Path
            Base directory containing sources/, normalized/, igblast/
        """
        self.base_dir = base_dir
        self.sources_dir = base_dir / "sources"
        self.normalized_dir = base_dir / "normalized"
        self.igblast_dir = base_dir / "igblast"

    def update(self, species: str) -> None:
        """
        Update pipeline for species. Auto-detects what needs rebuilding.

        Parameters
        ----------
        species : str
            Species name (e.g., "human", "mouse")
        """
        logger.info(f"Checking {species} germlines (offline mode)...")

        ready, error_msg = _check_offline_ready(self.sources_dir, species)
        if not ready:
            logger.error(f"offline_check_failed species={species}")
            raise RuntimeError(error_msg)

        if self._sources_changed(species):
            logger.info("Sources changed, rebuilding normalized data...")
            self._rebuild_normalized(species)

        if self._normalized_changed(species):
            logger.info("Normalized data changed, rebuilding IgBLAST...")
            self._rebuild_igblast(species)

        logger.info(f"{species} germlines up to date!")

    def force_rebuild(self, species: str) -> None:
        """
        Force complete rebuild regardless of timestamps.

        Parameters
        ----------
        species : str
            Species name
        """
        start = time.time()
        logger.info(f"Force rebuilding {species}...")
        self._rebuild_normalized(species)
        self._rebuild_igblast(species)
        _log_timing("force_rebuild", start, species=species, status="complete")

    def _sources_changed(self, species: str) -> bool:
        """
        Check if any source FASTA is newer than normalized.

        Parameters
        ----------
        species : str
            Species name

        Returns
        -------
        bool
            True if sources changed
        """
        normalized_dir = self.normalized_dir / species / "gapped"

        # Guard: normalized never built
        if not normalized_dir.exists():
            return True

        normalized_files = list(normalized_dir.glob("*.fasta"))
        if not normalized_files:
            return True

        latest_normalized = self._get_latest_mtime(normalized_files)

        # Check all source directories
        for source_name in ["custom", "ogrdb", "vdjbase", "imgt"]:
            if self._source_newer_than(source_name, species, latest_normalized):
                return True

        return False

    def _source_newer_than(
        self,
        source_name: str,
        species: str,
        threshold: float
    ) -> bool:
        """
        Check if source directory has files newer than threshold.

        Parameters
        ----------
        source_name : str
            Source name (imgt, ogrdb, custom)
        species : str
            Species name
        threshold : float
            Timestamp threshold

        Returns
        -------
        bool
            True if source has newer files
        """
        source_dir = self.sources_dir / source_name / species

        # Guard: source doesn't exist
        if not source_dir.exists():
            return False

        for fasta_file in source_dir.glob("**/*.fasta"):
            if fasta_file.stat().st_mtime > threshold:
                file_hash = _file_hash(fasta_file)
                logger.info(
                    f"change_detected path={fasta_file} hash={file_hash} "
                    f"change_type=modified provider={source_name}"
                )
                return True

        return False

    def _normalized_changed(self, species: str) -> bool:
        """
        Check if normalized files newer than IgBLAST files.

        Parameters
        ----------
        species : str
            Species name

        Returns
        -------
        bool
            True if normalized changed
        """
        igblast_db_dir = self.igblast_dir / "database" / species

        # Guard: igblast never built
        if not igblast_db_dir.exists():
            return True

        blast_files = list(igblast_db_dir.glob("*.nhr"))
        if not blast_files:
            return True

        latest_blast = self._get_latest_mtime(blast_files)

        # Check normalized
        normalized_dir = self.normalized_dir / species / "gapped"
        normalized_files = list(normalized_dir.glob("*.fasta"))

        if not normalized_files:
            return False

        latest_normalized = self._get_latest_mtime(normalized_files)

        return latest_normalized > latest_blast

    def _get_latest_mtime(self, files: List[Path]) -> float:
        """
        Get latest modification time from list of files.

        Parameters
        ----------
        files : List[Path]
            Files to check

        Returns
        -------
        float
            Latest modification timestamp
        """
        return max(f.stat().st_mtime for f in files)

    def _rebuild_normalized(self, species: str) -> None:
        """
        Rebuild normalized/ from sources/.

        Steps:
        1. Collect genes from all sources (priority-based)
        2. Write gapped FASTA files
        3. Write ungapped FASTA files

        Parameters
        ----------
        species : str
            Species name
        """
        start = time.time()
        from .manager import GermlineManager

        manager = GermlineManager()  # Uses default priority

        # Process each segment/chain combination
        for chain in CHAINS:
            for segment in SEGMENTS:
                self._normalize_segment(manager, species, chain, segment)
        _log_timing("rebuild_normalized", start, species=species)

    def _normalize_segment(
        self,
        manager,
        species: str,
        chain: str,
        segment: str
    ) -> None:
        """
        Normalize single segment/chain combination.

        Parameters
        ----------
        manager : GermlineManager
            Manager instance
        species : str
            Species name
        chain : str
            Chain type
        segment : str
            Segment type
        """
        genes = manager.get_genes(species, segment, chain)

        # Guard: no genes found
        if not genes:
            return

        logger.info(
            f"Processing {species} IG{chain}{segment}: {len(genes)} genes"
        )

        self._write_gapped_fasta(genes, species, chain, segment)
        self._write_ungapped_fasta(genes, species, chain, segment)

    def _write_gapped_fasta(
        self,
        genes: List[GermlineGene],
        species: str,
        chain: str,
        segment: str
    ) -> None:
        """
        Write gapped FASTA file.

        Parameters
        ----------
        genes : List[GermlineGene]
            Genes to write
        species : str
            Species name
        chain : str
            Chain type
        segment : str
            Segment type
        """
        gapped_records = [
            SeqRecord(
                seq=Seq(str(gene.sequence_gapped or gene.sequence)),
                id=gene.name,
                description=f"source={gene.source}"
            )
            for gene in genes
            if gene.sequence_gapped or segment == "D"
        ]

        # Guard: no gapped sequences
        if not gapped_records:
            return

        gapped_dir = self.normalized_dir / species / "gapped"
        gapped_dir.mkdir(parents=True, exist_ok=True)
        gapped_path = gapped_dir / f"IG{chain}{segment}.fasta"

        SeqIO.write(gapped_records, gapped_path, "fasta")
        logger.info(f"Wrote {len(gapped_records)} gapped to {gapped_path}")

    def _write_ungapped_fasta(
        self,
        genes: List[GermlineGene],
        species: str,
        chain: str,
        segment: str
    ) -> None:
        """
        Write ungapped FASTA file.

        Parameters
        ----------
        genes : List[GermlineGene]
            Genes to write
        species : str
            Species name
        chain : str
            Chain type
        segment : str
            Segment type
        """
        ungapped_records = [
            SeqRecord(
                seq=Seq(str(gene.sequence)),
                id=gene.name,
                description=f"source={gene.source}"
            )
            for gene in genes
        ]

        ungapped_dir = self.normalized_dir / species / "ungapped"
        ungapped_dir.mkdir(parents=True, exist_ok=True)
        ungapped_path = ungapped_dir / f"IG{chain}{segment}.fasta"

        SeqIO.write(ungapped_records, ungapped_path, "fasta")
        logger.info(f"Wrote {len(ungapped_records)} ungapped to {ungapped_path}")

    def _rebuild_igblast(self, species: str) -> None:
        """
        Rebuild igblast/ from normalized/.

        Steps:
        1. Build BLAST databases from ungapped
        2. Build aux file from gapped
        3. Generate organism.yaml

        Parameters
        ----------
        species : str
            Species name
        """
        start = time.time()
        from .builders.blast import BlastDBBuilder
        from .builders.aux import AuxFileBuilder

        blast_builder = BlastDBBuilder()
        blast_builder.build_for_species(
            species,
            source_dir=self.normalized_dir / species / "ungapped",
            output_dir=self.igblast_dir / "database" / species
        )

        aux_builder = AuxFileBuilder()
        aux_builder.build_for_species(
            species,
            source_dir=self.normalized_dir / species / "gapped",
            output_file=self.igblast_dir / "aux_db" / f"{species}_gl.aux"
        )

        self._generate_organism_yaml(species)
        _log_timing("rebuild_igblast", start, species=species)

    def _generate_organism_yaml(self, species: str) -> None:
        """
        Generate organism.yaml for IgBLAST internal_data.

        Parameters
        ----------
        species : str
            Species name
        """
        import yaml

        internal_data_dir = self.igblast_dir / "internal_data"
        internal_data_dir.mkdir(parents=True, exist_ok=True)

        yaml_path = internal_data_dir / "organism.yaml"

        if yaml_path.exists():
            with open(yaml_path) as f:
                data = yaml.safe_load(f) or {}
        else:
            data = {}

        if "organisms" not in data:
            data["organisms"] = {}

        data["organisms"][species] = {
            "database_path": f"../database/{species}",
            "aux_file": f"../aux_db/{species}_gl.aux",
            "segments": ["V", "D", "J"],
        }

        with open(yaml_path, "w") as f:
            yaml.dump(data, f, default_flow_style=False)

        logger.info(f"Updated organism.yaml for {species}")
