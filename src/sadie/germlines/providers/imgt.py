"""
IMGT Provider
=============

Provider for IMGT (International ImMunoGeneTics information system) database.

IMGT is the reference database for immunoglobulins, providing:
- Curated germline sequences
- IMGT-gapped sequences with standardized numbering
- Functional annotations
- Multiple species support

Data Source: https://www.imgt.org/
Reference Sequences: https://www.imgt.org/vquest/refseqh.html

TODO: Implement IMGT data download and parsing
Current Status: Stub implementation - reads from pre-downloaded FASTA files
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
from Bio import SeqIO

from .base import GermlineProvider
from ..models import GermlineGene, ProviderMetadata


logger = logging.getLogger(__name__)


class IMGTProvider(GermlineProvider):
    """
    Provider for IMGT germline database.

    IMGT provides curated, validated germline sequences with
    IMGT-standard gapped numbering.

    Expected Directory Structure:
        sources/imgt/
        ├── human/
        │   ├── IGHV.fasta    # Heavy chain V (IMGT-gapped)
        │   ├── IGHD.fasta    # Heavy chain D (ungapped)
        │   ├── IGHJ.fasta    # Heavy chain J (IMGT-gapped)
        │   ├── IGKV.fasta    # Kappa V (IMGT-gapped)
        │   ├── IGKJ.fasta    # Kappa J (IMGT-gapped)
        │   ├── IGLV.fasta    # Lambda V (IMGT-gapped)
        │   └── IGLJ.fasta    # Lambda J (IMGT-gapped)
        └── mouse/
            └── ...

    FASTA Format (from IMGT):
        >IGHV1-2*01|Homo sapiens|F|...
        cag.gtgcagctggtgcag...tctggggctgag...gtgaag...
        (dots indicate IMGT gaps)

    Examples
    --------
    >>> provider = IMGTProvider()
    >>> genes = provider.fetch_genes("human", "V", "H")
    >>> print(f"Found {len(genes)} IMGT IGHV genes")
    """

    def fetch_genes(
        self,
        species: str,
        segment: str,
        chain: str
    ) -> List[GermlineGene]:
        """
        Fetch IMGT genes from pre-downloaded FASTA files.

        TODO: Implement automatic download from IMGT if files missing

        Parameters
        ----------
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type

        Returns
        -------
        List[GermlineGene]
            IMGT genes
        """
        fasta_path = self.get_fasta_path(species, segment, chain)
        gapped_path = self._get_gapped_fasta_path(species, segment, chain)

        # Guard: file doesn't exist
        if not fasta_path.exists():
            logger.debug(f"No IMGT file: {fasta_path}")
            logger.info(
                "Run download script or add FASTA manually. "
                "See sources/imgt/README.md"
            )
            return []

        # Load gapped sequences if available
        gapped_sequences: Dict[str, str] = {}
        if gapped_path.exists():
            gapped_sequences = self._load_gapped_sequences(gapped_path)
            logger.debug(f"Loaded {len(gapped_sequences)} gapped sequences from {gapped_path}")

        logger.info(f"Loading IMGT FASTA: {fasta_path}")

        genes = self._parse_imgt_fasta(fasta_path, species, segment, chain, gapped_sequences)

        logger.info(f"Loaded {len(genes)} IMGT genes")

        return genes

    def _parse_imgt_fasta(
        self,
        fasta_path: Path,
        species: str,
        segment: str,
        chain: str,
        gapped_sequences: Optional[Dict[str, str]] = None
    ) -> List[GermlineGene]:
        """
        Parse IMGT FASTA file.

        IMGT FASTA headers contain metadata separated by pipes.
        IMGT sequences are typically pre-gapped with dots.

        TODO: Parse IMGT header metadata (functionality, accession, etc.)

        Parameters
        ----------
        fasta_path : Path
            Path to FASTA file
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type
        gapped_sequences : Dict[str, str], optional
            Pre-loaded gapped sequences mapping gene name to gapped sequence

        Returns
        -------
        List[GermlineGene]
            Parsed genes
        """
        genes = []
        gapped_sequences = gapped_sequences or {}

        try:
            records = list(SeqIO.parse(fasta_path, "fasta"))
        except Exception as e:
            logger.error(f"Failed to parse {fasta_path}: {e}")
            return []

        for record in records:
            gene = self._create_imgt_gene(record, species, segment, chain, gapped_sequences)
            if gene:
                genes.append(gene)

        return genes

    def _create_imgt_gene(
        self,
        record,
        species: str,
        segment: str,
        chain: str,
        gapped_sequences: Optional[Dict[str, str]] = None
    ) -> Optional[GermlineGene]:
        """
        Create GermlineGene from IMGT SeqRecord.

        IMGT sequences are typically gapped with dots.
        Header format: >GENENAME|Species|Functionality|...

        TODO: Parse all IMGT metadata fields

        Parameters
        ----------
        record : SeqRecord
            BioPython sequence record
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type
        gapped_sequences : Dict[str, str], optional
            Pre-loaded gapped sequences from *_gapped.fasta file

        Returns
        -------
        GermlineGene or None
            Gene object if successful
        """
        gapped_sequences = gapped_sequences or {}

        # Parse header
        # IMGT format: >ACCESSION|GENE_NAME|SPECIES|FUNCTIONALITY|REGION|...
        header_parts = record.id.split("|")
        gene_name = header_parts[1] if len(header_parts) > 1 else header_parts[0]

        # Extract functionality if present (index 3 in IMGT format)
        functionality = "F"  # Default to functional
        if len(header_parts) > 3:
            functionality = header_parts[3]

        # Get sequence from main file
        sequence_raw = str(record.seq).upper()
        is_gapped = "." in sequence_raw

        if is_gapped:
            # Main file has gapped sequences (e.g., human)
            sequence_gapped = sequence_raw
            sequence_ungapped = sequence_raw.replace(".", "").replace("-", "")
        else:
            # Main file is ungapped, look up from gapped file
            sequence_ungapped = sequence_raw
            sequence_gapped = gapped_sequences.get(gene_name)

        # Determine if functional
        is_functional = functionality == "F"

        try:
            gene = GermlineGene(
                name=gene_name,
                species=species,
                segment=segment,
                chain=chain,
                sequence=sequence_ungapped,
                sequence_gapped=sequence_gapped,
                is_functional=is_functional,
                functionality=functionality,
                source="imgt",
            )
            return gene

        except Exception as e:
            logger.error(f"Failed to create IMGT gene {gene_name}: {e}")
            return None

    def _get_gapped_fasta_path(
        self,
        species: str,
        segment: str,
        chain: str
    ) -> Path:
        """
        Get path to gapped FASTA file.

        Parameters
        ----------
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type

        Returns
        -------
        Path
            Path to gapped FASTA file
        """
        return self.data_dir / species / f"IG{chain}{segment}_gapped.fasta"

    def _load_gapped_sequences(self, fasta_path: Path) -> Dict[str, str]:
        """
        Load gapped sequences from FASTA file.

        Parameters
        ----------
        fasta_path : Path
            Path to gapped FASTA file

        Returns
        -------
        Dict[str, str]
            Mapping of gene name to gapped sequence
        """
        gapped = {}
        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                # IMGT format: >ACCESSION|GENE_NAME|...
                header_parts = record.id.split("|")
                gene_name = header_parts[1] if len(header_parts) > 1 else header_parts[0]
                gapped[gene_name] = str(record.seq).upper()
        except Exception as e:
            logger.warning(f"Failed to load gapped sequences from {fasta_path}: {e}")
        return gapped

    def fetch_gene_by_name(
        self,
        name: str,
        species: str
    ) -> Optional[GermlineGene]:
        """
        Fetch specific IMGT gene by name.

        Parameters
        ----------
        name : str
            Gene name
        species : str
            Species name

        Returns
        -------
        GermlineGene or None
            Gene if found
        """
        for segment in ["V", "D", "J"]:
            for chain in ["H", "K", "L"]:
                genes = self.fetch_genes(species, segment, chain)
                for gene in genes:
                    if gene.name == name:
                        return gene

        return None

    def is_available(self, species: str) -> bool:
        """
        Check if IMGT data is available for species.

        Parameters
        ----------
        species : str
            Species name

        Returns
        -------
        bool
            True if data available
        """
        species_dir = self.data_dir / species

        if not species_dir.exists():
            return False

        fasta_files = list(species_dir.glob("IG*.fasta"))
        return len(fasta_files) > 0

    def get_metadata(self) -> ProviderMetadata:
        """
        Get IMGT provider metadata.

        Returns
        -------
        ProviderMetadata
            Provider metadata
        """
        species = []
        if self.data_dir.exists():
            for species_dir in self.data_dir.iterdir():
                if species_dir.is_dir():
                    if list(species_dir.glob("IG*.fasta")):
                        species.append(species_dir.name)

        return ProviderMetadata(
            name="imgt",
            version="release-202501",  # TODO: Extract from data
            last_updated=datetime.now(),
            species_available=species,
            url="https://www.imgt.org/",
        )

    def download(self, species: List[str], force: bool = False) -> None:
        """
        Download IMGT data for species.

        Uses the IMGTDownloader to fetch germline sequences from
        the IMGT V-QUEST reference directory.

        Parameters
        ----------
        species : List[str]
            Species to download. If empty, downloads all available species.
        force : bool
            Force re-download even if files exist
        """
        from ..scripts.download_imgt import IMGTDownloader, SPECIES_MAP

        downloader = IMGTDownloader(output_dir=self.data_dir)

        if not species:
            species = list(SPECIES_MAP.keys())

        logger.info(f"Downloading IMGT data for {len(species)} species...")

        results = downloader.download(species, force=force)

        total_seqs = sum(results.values())
        logger.info(f"Downloaded {total_seqs} sequences for {len(results)} species")
