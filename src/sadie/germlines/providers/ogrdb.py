"""
OGRDB Provider
==============

Provider for OGRDB (Open Germline Receptor Database).

OGRDB is a community-curated database providing:
- Novel alleles discovered through repertoire sequencing
- Inferred germline sequences with supporting evidence
- AIRR Community governance and review

Data Source: https://ogrdb.airr-community.org/
Zenodo Archive: https://zenodo.org/records/18145568

The OGRDB Zenodo archive contains:
- SQL dump with gene_description table containing:
  - sequence: ungapped nucleotide sequence
  - coding_seq_imgt: IMGT-gapped nucleotide sequence

Download data using:
    python -m sadie.germlines.scripts.download_ogrdb --species human
"""

import logging
import time
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
from Bio import SeqIO

from .base import GermlineProvider
from ..models import GermlineGene, ProviderMetadata


logger = logging.getLogger(__name__)


class OGRDBProvider(GermlineProvider):
    """
    Provider for OGRDB germline database.

    OGRDB provides novel alleles and inferred germline sequences
    discovered through repertoire sequencing.

    Expected Directory Structure:
        sources/ogrdb/
        ├── human/
        │   ├── IGHV.fasta          # Ungapped sequences
        │   ├── IGHV_gapped.fasta   # IMGT-gapped sequences (optional)
        │   ├── IGHD.fasta
        │   ├── IGHJ.fasta
        │   ├── IGHJ_gapped.fasta
        │   ├── IGKV.fasta
        │   ├── IGKJ.fasta
        │   ├── IGLV.fasta
        │   └── IGLJ.fasta
        └── mouse/
            └── ...

    Examples
    --------
    >>> provider = OGRDBProvider()
    >>> genes = provider.fetch_genes("human", "V", "H")
    >>> print(f"Found {len(genes)} OGRDB IGHV genes")
    
    >>> # Download data from Zenodo archive
    >>> provider.download(["human", "mouse"])
    """

    def fetch_genes(
        self,
        species: str,
        segment: str,
        chain: str
    ) -> List[GermlineGene]:
        """
        Fetch OGRDB genes from downloaded FASTA files.

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
            OGRDB genes
        """
        fasta_path = self.get_fasta_path(species, segment, chain)
        gapped_path = self._get_gapped_fasta_path(species, segment, chain)

        # Guard: file doesn't exist
        if not fasta_path.exists():
            logger.debug(f"No OGRDB file: {fasta_path}")
            logger.info(
                "Run 'python -m sadie.germlines.scripts.download_ogrdb --species {species}' "
                "or add FASTA manually. See sources/ogrdb/OGRDB_DATA.md"
            )
            return []

        # Load gapped sequences if available
        gapped_sequences: Dict[str, str] = {}
        if gapped_path.exists():
            gapped_sequences = self._load_gapped_sequences(gapped_path)
            logger.debug(f"Loaded {len(gapped_sequences)} gapped sequences from {gapped_path}")

        logger.info(f"Loading OGRDB FASTA: {fasta_path}")

        genes = self._parse_ogrdb_fasta(fasta_path, species, segment, chain, gapped_sequences)

        logger.info(
            f"operation=load_ogrdb provider=ogrdb "
            f"species={species} segment={segment} chain={chain} "
            f"gene_count={len(genes)} status=success"
        )

        return genes

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
                gene_name = record.id.split("|")[0]
                gapped[gene_name] = str(record.seq).upper()
        except Exception as e:
            logger.warning(f"Failed to load gapped sequences from {fasta_path}: {e}")
        return gapped

    def _parse_ogrdb_fasta(
        self,
        fasta_path: Path,
        species: str,
        segment: str,
        chain: str,
        gapped_sequences: Optional[Dict[str, str]] = None
    ) -> List[GermlineGene]:
        """
        Parse OGRDB FASTA file.

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
            gene = self._create_ogrdb_gene(
                record, species, segment, chain, gapped_sequences
            )
            if gene:
                genes.append(gene)

        return genes

    def _create_ogrdb_gene(
        self,
        record,
        species: str,
        segment: str,
        chain: str,
        gapped_sequences: Optional[Dict[str, str]] = None
    ) -> Optional[GermlineGene]:
        """
        Create GermlineGene from OGRDB SeqRecord.

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
            Pre-loaded gapped sequences

        Returns
        -------
        GermlineGene or None
            Gene object if successful
        """
        gene_name = record.id.split("|")[0]
        gapped_sequences = gapped_sequences or {}

        sequence = str(record.seq).upper()
        is_gapped = "." in sequence or "-" in sequence

        if is_gapped:
            sequence_gapped = sequence
            sequence_ungapped = sequence.replace(".", "").replace("-", "")
        else:
            sequence_ungapped = sequence
            # Look up gapped sequence from pre-loaded file
            sequence_gapped = gapped_sequences.get(gene_name)

        try:
            gene = GermlineGene(
                name=gene_name,
                species=species,
                segment=segment,
                chain=chain,
                sequence=sequence_ungapped,
                sequence_gapped=sequence_gapped,
                is_functional=True,  # OGRDB sequences are typically functional
                functionality="F",
                source="ogrdb",
            )
            return gene

        except Exception as e:
            logger.error(f"Failed to create OGRDB gene {gene_name}: {e}")
            return None

    def fetch_gene_by_name(
        self,
        name: str,
        species: str
    ) -> Optional[GermlineGene]:
        """
        Fetch specific OGRDB gene by name.

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
        Check if OGRDB data is available for species.

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
        Get OGRDB provider metadata.

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
            name="ogrdb",
            version="latest",  # TODO: Extract from OGRDB API
            last_updated=datetime.now(),
            species_available=species,
            url="https://ogrdb.airr-community.org/",
        )

    def download(self, species: List[str]) -> None:
        """
        Download OGRDB data from Zenodo archive.

        Downloads the OGRDB archive from Zenodo and extracts
        gapped and ungapped FASTA sequences for the specified species.

        Parameters
        ----------
        species : List[str]
            Species to download (e.g., ["human", "mouse"])

        Examples
        --------
        >>> provider = OGRDBProvider()
        >>> provider.download(["human"])
        """
        from ..scripts.download_ogrdb import OGRDBDownloader
        
        start_time = time.time()
        
        downloader = OGRDBDownloader(output_dir=self.data_dir)
        downloader.download(species)
        
        duration_ms = int((time.time() - start_time) * 1000)
        logger.info(
            f"operation=download provider=ogrdb "
            f"species={','.join(species)} duration_ms={duration_ms} status=success"
        )
