"""
VDJbase Provider
================

Provider for VDJbase (Adaptive Immune Receptor Gene Database).

VDJbase provides population-specific germline alleles from:
- Repertoire sequencing data
- Genotype inference pipelines
- Multiple population studies

Data Source: https://vdjbase.org/
API Documentation: https://vdjbase.org/admin/api/

Supported Species:
- Human (IGH, IGK, IGL)
- Rhesus Macaque (IGH, IGK, IGL)

Note: VDJbase sequences are UNGAPPED - require gapper module for IMGT formatting.
"""

import logging
import json
import time
from pathlib import Path
from typing import List, Optional, Dict, Any
from datetime import datetime
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
from urllib.parse import quote
from Bio import SeqIO
from io import StringIO

from .base import GermlineProvider
from ..models import GermlineGene, ProviderMetadata


logger = logging.getLogger(__name__)


# VDJbase API configuration
VDJBASE_API_BASE = "https://vdjbase.org/api"
VDJBASE_ADMIN_API_BASE = "https://vdjbase.org/admin/api"

# Species name mapping (VDJbase API names to internal names)
SPECIES_MAP = {
    "Human": "human",
    "Rhesus Macaque": "rhesus_macaque",
    "Mouse": "mouse",
}

SPECIES_MAP_REVERSE = {v: k for k, v in SPECIES_MAP.items()}

# Chain to dataset mapping
CHAIN_TO_DATASET = {
    "H": "IGH",
    "K": "IGK", 
    "L": "IGL",
}


class VDJbaseProvider(GermlineProvider):
    """
    Provider for VDJbase germline database.

    VDJbase provides population-specific germline alleles discovered
    through repertoire sequencing and genotype inference.

    Expected Directory Structure:
        sources/vdjbase/
        ├── human/
        │   ├── IGHV.fasta
        │   ├── IGHD.fasta
        │   ├── IGHJ.fasta
        │   ├── IGKV.fasta
        │   ├── IGKJ.fasta
        │   ├── IGLV.fasta
        │   └── IGLJ.fasta
        ├── rhesus_macaque/
        │   └── ...
        └── mouse/
            └── ...

    Examples
    --------
    >>> provider = VDJbaseProvider()
    >>> genes = provider.fetch_genes("human", "V", "H")
    >>> print(f"Found {len(genes)} VDJbase IGHV genes")
    
    >>> # Download data from VDJbase API
    >>> provider.download(["human", "rhesus_macaque"])
    """

    def __init__(self, data_dir: Optional[Path] = None):
        """
        Initialize VDJbase provider.

        Parameters
        ----------
        data_dir : Path, optional
            Base directory for VDJbase data.
            Defaults to sources/vdjbase/
        """
        if data_dir is None:
            data_dir = Path(__file__).parent.parent / "sources" / "vdjbase"
        
        super().__init__(data_dir)
        self.name = "vdjbase"
        self._api_cache: Dict[str, Any] = {}

    def fetch_genes(
        self,
        species: str,
        segment: str,
        chain: str
    ) -> List[GermlineGene]:
        """
        Fetch VDJbase genes from local FASTA files.

        Parameters
        ----------
        species : str
            Species name (e.g., "human", "rhesus_macaque")
        segment : str
            Segment type: "V", "D", or "J"
        chain : str
            Chain type: "H", "K", or "L"

        Returns
        -------
        List[GermlineGene]
            VDJbase genes
        """
        fasta_path = self.get_fasta_path(species, segment, chain)

        if not fasta_path.exists():
            logger.warning(
                f"VDJbase data not found at {fasta_path}. "
                "Skipping VDJbase provider."
            )
            return []

        logger.info(f"Loading VDJbase FASTA: {fasta_path}")
        
        genes = self._parse_vdjbase_fasta(fasta_path, species, segment, chain)
        
        logger.info(
            f"operation=load_vdjbase provider=vdjbase "
            f"species={species} segment={segment} chain={chain} "
            f"gene_count={len(genes)} status=success"
        )

        return genes

    def _parse_vdjbase_fasta(
        self,
        fasta_path: Path,
        species: str,
        segment: str,
        chain: str
    ) -> List[GermlineGene]:
        """
        Parse VDJbase FASTA file.

        VDJbase FASTA header format:
        >{gene_name}|{species}|{segment}|{chain}[|population={pop}][|genotype={gt}]

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

        Returns
        -------
        List[GermlineGene]
            Parsed genes
        """
        genes = []

        try:
            records = list(SeqIO.parse(fasta_path, "fasta"))
        except Exception as e:
            logger.error(f"Failed to parse {fasta_path}: {e}")
            return []

        for record in records:
            gene = self._create_vdjbase_gene(record, species, segment, chain)
            if gene:
                genes.append(gene)

        return genes

    def _create_vdjbase_gene(
        self,
        record,
        species: str,
        segment: str,
        chain: str
    ) -> Optional[GermlineGene]:
        """
        Create GermlineGene from VDJbase SeqRecord.

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

        Returns
        -------
        GermlineGene or None
            Gene object if successful
        """
        # Parse header - VDJbase format or simple gene name
        header_parts = record.description.split("|")
        gene_name = header_parts[0].strip()
        
        # Extract metadata from header if present
        metadata = {}
        for part in header_parts[1:]:
            if "=" in part:
                key, value = part.split("=", 1)
                metadata[key.strip()] = value.strip()

        sequence = str(record.seq).upper()
        
        # VDJbase sequences are typically ungapped
        is_gapped = "." in sequence or "-" in sequence

        if is_gapped:
            sequence_gapped = sequence
            sequence_ungapped = sequence.replace(".", "").replace("-", "")
        else:
            sequence_ungapped = sequence
            sequence_gapped = None  # Will be gapped by gapper module

        # Extract novel/confidence flags if present in description
        is_novel = "novel" in record.description.lower()
        low_confidence = "low_confidence" in record.description.lower()

        try:
            gene = GermlineGene(
                name=gene_name,
                species=species,
                segment=segment,
                chain=chain,
                sequence=sequence_ungapped,
                sequence_gapped=sequence_gapped,
                is_functional=True,
                functionality="F",
                source="vdjbase",
            )
            return gene

        except Exception as e:
            logger.error(f"Failed to create VDJbase gene {gene_name}: {e}")
            return None

    def fetch_gene_by_name(
        self,
        name: str,
        species: str
    ) -> Optional[GermlineGene]:
        """
        Fetch specific VDJbase gene by name.

        Parameters
        ----------
        name : str
            Gene name (e.g., "IGHV1-69*01")
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
        Check if VDJbase data is available for species.

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
        Get VDJbase provider metadata.

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
            name="vdjbase",
            version="latest",
            last_updated=datetime.now(),
            species_available=species,
            url="https://vdjbase.org/",
        )

    def download(self, species: List[str]) -> None:
        """
        Download VDJbase data for specified species.

        Fetches germline sequences from VDJbase API and saves
        to local FASTA files.

        Parameters
        ----------
        species : List[str]
            Species to download (e.g., ["human", "rhesus_macaque"])
        """
        start_time = time.time()
        
        for sp in species:
            logger.info(f"Downloading VDJbase data for {sp}...")
            self._download_species(sp)
        
        duration_ms = int((time.time() - start_time) * 1000)
        logger.info(
            f"operation=download provider=vdjbase "
            f"species={','.join(species)} duration_ms={duration_ms} status=success"
        )

    def _download_species(self, species: str) -> None:
        """
        Download all chains for a species.

        Parameters
        ----------
        species : str
            Species name (internal format, e.g., "human")
        """
        # Create species directory
        species_dir = self.data_dir / species
        species_dir.mkdir(parents=True, exist_ok=True)

        # Map to VDJbase API species name
        api_species = SPECIES_MAP_REVERSE.get(species, species)
        
        # Download each chain
        for chain, dataset in CHAIN_TO_DATASET.items():
            self._download_chain(species, api_species, chain, dataset)

    def _download_chain(
        self,
        species: str,
        api_species: str,
        chain: str,
        dataset: str
    ) -> None:
        """
        Download sequences for a specific chain.

        Parameters
        ----------
        species : str
            Internal species name
        api_species : str
            VDJbase API species name
        chain : str
            Chain type (H, K, L)
        dataset : str
            Dataset name (IGH, IGK, IGL)
        """
        logger.info(f"Fetching {dataset} for {api_species}...")
        
        # Fetch sequences from API
        sequences = self._fetch_sequences_from_api(api_species, dataset)
        
        if not sequences:
            logger.warning(f"No sequences found for {api_species} {dataset}")
            return

        # Deduplicate by gene name - keep first occurrence (highest confidence)
        unique_sequences = {}
        for seq_data in sequences:
            name = seq_data.get("name", "")
            if name and name not in unique_sequences:
                unique_sequences[name] = seq_data
        
        logger.info(
            f"Deduplicated {len(sequences)} sequences to {len(unique_sequences)} unique alleles"
        )

        # Group by segment type
        segments = {"V": [], "D": [], "J": []}
        
        for seq_data in unique_sequences.values():
            seq_type = seq_data.get("type", "")
            if seq_type.endswith("V"):
                segments["V"].append(seq_data)
            elif seq_type.endswith("D"):
                segments["D"].append(seq_data)
            elif seq_type.endswith("J"):
                segments["J"].append(seq_data)

        # Write FASTA files for each segment
        for segment, seqs in segments.items():
            if seqs:
                self._write_fasta(species, chain, segment, seqs)

    def _fetch_sequences_from_api(
        self,
        api_species: str,
        dataset: str,
        page_size: int = 1000,
        max_pages: int = 100
    ) -> List[Dict[str, Any]]:
        """
        Fetch all sequences from VDJbase API with pagination.

        Parameters
        ----------
        api_species : str
            VDJbase API species name
        dataset : str
            Dataset name
        page_size : int
            Number of items per page
        max_pages : int
            Maximum number of pages to fetch (safety limit)

        Returns
        -------
        List[Dict]
            All sequence records
        """
        all_sequences = []
        page = 1
        
        while page <= max_pages:
            url = (
                f"{VDJBASE_ADMIN_API_BASE}/repseq/sequences/"
                f"{quote(api_species)}/{dataset}?page={page}&per_page={page_size}"
            )
            
            try:
                logger.debug(f"Fetching: {url}")
                
                request = Request(url)
                request.add_header("Accept", "application/json")
                
                with urlopen(request, timeout=30) as response:
                    data = json.loads(response.read().decode("utf-8"))
                
                samples = data.get("samples", [])
                
                if not samples:
                    break
                
                all_sequences.extend(samples)
                
                logger.info(
                    f"Downloaded {len(all_sequences)} sequences "
                    f"(page {page}/{max_pages})"
                )
                
                # Check if we've reached the last page
                if len(samples) < page_size:
                    break
                
                page += 1
                
                # Rate limiting - minimal delay
                time.sleep(0.1)
                
            except HTTPError as e:
                if e.code == 404:
                    logger.warning(f"No data available for {api_species}/{dataset}")
                    break
                else:
                    logger.error(f"HTTP error fetching {url}: {e}")
                    break
            except URLError as e:
                logger.error(f"URL error fetching {url}: {e}")
                break
            except Exception as e:
                logger.error(f"Error fetching {url}: {e}")
                break

        if page > max_pages:
            logger.warning(
                f"Reached max pages limit ({max_pages}). "
                f"Some sequences may be missing."
            )

        return all_sequences

    def _write_fasta(
        self,
        species: str,
        chain: str,
        segment: str,
        sequences: List[Dict[str, Any]]
    ) -> None:
        """
        Write sequences to FASTA file.

        Parameters
        ----------
        species : str
            Species name
        chain : str
            Chain type
        segment : str
            Segment type
        sequences : List[Dict]
            Sequence records from API
        """
        fasta_path = self.get_fasta_path(species, segment, chain)
        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Writing {len(sequences)} sequences to {fasta_path}")
        
        with open(fasta_path, "w") as f:
            for seq_data in sequences:
                name = seq_data.get("name", "unknown")
                sequence = seq_data.get("seq", "")
                
                if not sequence:
                    continue
                
                # Build header with metadata
                # Format: >gene_name|species|segment|chain[|metadata]
                header_parts = [name]
                
                # Add optional metadata
                if seq_data.get("novel"):
                    header_parts.append("novel=true")
                if seq_data.get("low_confidence"):
                    header_parts.append("low_confidence=true")
                if seq_data.get("appears"):
                    header_parts.append(f"appears={seq_data['appears']}")
                
                header = "|".join(header_parts)
                
                f.write(f">{header}\n")
                f.write(f"{sequence.upper()}\n")
        
        logger.info(f"Wrote {fasta_path}")

    def get_available_species(self) -> List[str]:
        """
        Get list of species available from VDJbase API.

        Returns
        -------
        List[str]
            Available species names (internal format)
        """
        try:
            url = f"{VDJBASE_ADMIN_API_BASE}/genomic/species"
            
            request = Request(url)
            request.add_header("Accept", "application/json")
            
            with urlopen(request, timeout=10) as response:
                data = json.loads(response.read().decode("utf-8"))
            
            # Map API names to internal names
            return [
                SPECIES_MAP.get(sp, sp.lower().replace(" ", "_"))
                for sp in data
            ]
            
        except Exception as e:
            logger.error(f"Failed to fetch species list: {e}")
            return []

    def get_available_datasets(self, species: str) -> List[str]:
        """
        Get available datasets for a species.

        Parameters
        ----------
        species : str
            Species name (internal format)

        Returns
        -------
        List[str]
            Available dataset names (IGH, IGK, IGL)
        """
        api_species = SPECIES_MAP_REVERSE.get(species, species)
        
        try:
            url = f"{VDJBASE_ADMIN_API_BASE}/genomic/data_sets/{quote(api_species)}"
            
            request = Request(url)
            request.add_header("Accept", "application/json")
            
            with urlopen(request, timeout=10) as response:
                data = json.loads(response.read().decode("utf-8"))
            
            return [ds.get("dataset", "") for ds in data if ds.get("dataset")]
            
        except Exception as e:
            logger.error(f"Failed to fetch datasets for {species}: {e}")
            return []
