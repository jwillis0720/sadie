"""
Custom Germline Provider
=========================

Allows users to add their own germline sequences by dropping FASTA files
into the custom/ directory.

Features:
- Auto-detects ungapped vs gapped sequences
- Auto-aligns ungapped sequences to IMGT reference using GapperService
- Validates sequences
- Works offline

Directory Structure:
    sources/custom/
    ├── README.md          # Instructions for users
    ├── human/
    │   ├── IGHV.fasta    # User adds sequences here
    │   ├── IGHD.fasta
    │   └── .processed/   # Auto-generated metadata
    └── mouse/
        └── ...

FASTA Format (flexible):
    Headers: >IGHV1-69*01 or >my_novel_seq
    Sequences: Ungapped (ACGT) or IMGT-gapped (with dots)
"""

import logging
from pathlib import Path
from typing import List, Optional
from datetime import datetime
from Bio import SeqIO

from .base import GermlineProvider
from ..models import GermlineGene, ProviderMetadata
from ..builders.gapper import GapperService


logger = logging.getLogger(__name__)

VALID_NUCLEOTIDES = set("ACGTNacgtn.")
IUPAC_AMBIGUOUS = set("RYSWKMBDHVryswkmbdhv-")


def _validate_sequence(sequence: str, gene_name: str) -> tuple[bool, str]:
    if not sequence:
        return False, f"Empty sequence for {gene_name}"
    clean_seq = sequence.replace(".", "").replace("-", "")
    if not clean_seq:
        return False, f"Sequence contains only gap characters for {gene_name}"
    valid_chars = VALID_NUCLEOTIDES | IUPAC_AMBIGUOUS
    invalid = set(sequence) - valid_chars
    if invalid:
        return False, f"Invalid characters in {gene_name}: {invalid}"
    return True, ""


class CustomProvider(GermlineProvider):
    """
    Provider for user-supplied custom germline sequences.

    This is the highest priority provider by default, allowing users
    to override standard databases with their novel sequences.

    Examples
    --------
    >>> # User adds FASTA file to sources/custom/human/IGHV.fasta
    >>> provider = CustomProvider()
    >>> genes = provider.fetch_genes("human", "V", "H")
    >>> print(f"Found {len(genes)} custom genes")
    """

    def __init__(self, data_dir: Optional[Path] = None, template_dir: Optional[Path] = None):
        """
        Initialize CustomProvider with optional gapping template directory.

        Parameters
        ----------
        data_dir : Path, optional
            Directory containing custom FASTA files.
            Defaults to sources/custom/
        template_dir : Path, optional
            Directory containing IMGT-gapped templates for auto-gapping.
            Defaults to sources/imgt/ (sibling of custom/)
        """
        super().__init__(data_dir=data_dir)
        
        # Set up template directory for gapping
        # Template dir should be sibling to custom (both under sources/)
        if template_dir is None:
            # self.data_dir is sources/custom, so parent is sources, then add imgt
            template_dir = self.data_dir.parent / "imgt"
        self.template_dir = template_dir
        
        # Lazy initialization of gapper (per-species)
        self._gappers: dict = {}

    def _get_gapper(self, species: str) -> tuple:
        """
        Get or create a GapperService for the given species.

        Falls back to human IMGT templates if species-specific templates
        aren't available, since IG gene structure is conserved across species.

        Parameters
        ----------
        species : str
            Species name

        Returns
        -------
        tuple
            (GapperService, template_species) - gapper and species to use for templates
        """
        cache_key = species
        
        if cache_key not in self._gappers:
            species_template_dir = self.template_dir / species
            if species_template_dir.exists():
                # Use species-specific templates
                self._gappers[cache_key] = (GapperService(template_dir=self.template_dir), species)
                logger.info(f"Initialized gapper with {species} templates")
            else:
                # Fall back to human templates - IG gene structure is conserved
                human_template_dir = self.template_dir / "human"
                if human_template_dir.exists():
                    self._gappers[cache_key] = (GapperService(template_dir=self.template_dir), "human")
                    logger.info(f"Using human IMGT templates for {species} (no species-specific templates)")
                else:
                    # No templates available at all
                    self._gappers[cache_key] = (GapperService(template_dir=None), species)
                    logger.warning(f"No IMGT templates found, gapping disabled for {species}")
        
        return self._gappers[cache_key]

    def fetch_genes(
        self,
        species: str,
        segment: str,
        chain: str
    ) -> List[GermlineGene]:
        """
        Fetch custom genes from user-supplied FASTA files.

        Automatically processes ungapped sequences if needed.

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
            Custom genes
        """
        fasta_path = self.get_fasta_path(species, segment, chain)

        # Guard: file doesn't exist
        if not fasta_path.exists():
            logger.debug(f"No custom file: {fasta_path}")
            return []

        logger.info(f"Processing custom FASTA: {fasta_path}")

        genes = self._parse_fasta_file(fasta_path, species, segment, chain)

        logger.info(f"Loaded {len(genes)} custom genes")

        return genes

    def _parse_fasta_file(
        self,
        fasta_path: Path,
        species: str,
        segment: str,
        chain: str
    ) -> List[GermlineGene]:
        """
        Parse FASTA file and create GermlineGene objects.

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
            gene = self._create_gene_from_record(
                record,
                species,
                segment,
                chain
            )
            if gene:
                genes.append(gene)

        return genes

    def _create_gene_from_record(
        self,
        record,
        species: str,
        segment: str,
        chain: str
    ) -> Optional[GermlineGene]:
        """
        Create GermlineGene from SeqRecord.

        Handles both gapped and ungapped sequences. For ungapped sequences,
        attempts to auto-gap using the GapperService with IMGT templates.

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
        sequence = str(record.seq).upper()
        gene_name = self._clean_gene_name(record.id, segment, chain)

        valid, error_msg = _validate_sequence(sequence, gene_name)
        if not valid:
            logger.warning(f"Skipping invalid sequence: {error_msg}")
            return None

        is_gapped = "." in sequence or "-" in sequence

        # Determine gapped/ungapped versions
        if is_gapped:
            sequence_gapped = sequence
            sequence_ungapped = sequence.replace(".", "").replace("-", "")
        else:
            sequence_ungapped = sequence
            # Auto-gap using GapperService
            gapper, template_species = self._get_gapper(species)
            sequence_gapped = gapper.gap_sequence(
                sequence=sequence_ungapped,
                segment=segment,
                chain=chain,
                gene_name=gene_name,
                species=template_species  # Use template species for lookup
            )
            if sequence_gapped:
                logger.debug(f"Auto-gapped {gene_name} using {template_species} templates")
            else:
                logger.debug(f"Could not auto-gap {gene_name}, storing ungapped only")

        try:
            gene = GermlineGene(
                name=gene_name,
                species=species,
                segment=segment,
                chain=chain,
                sequence=sequence_ungapped,
                sequence_gapped=sequence_gapped,
                is_functional=True,  # Assume custom sequences are functional
                functionality="F",
                source="custom",
            )
            return gene

        except Exception as e:
            logger.error(f"Failed to create gene {gene_name}: {e}")
            return None

    def _clean_gene_name(
        self,
        record_id: str,
        segment: str,
        chain: str
    ) -> str:
        """
        Clean and standardize gene names from FASTA headers.

        Handles various formats:
        - IGHV1-69*01
        - IGHV1-69*01|description
        - my_custom_sequence

        Parameters
        ----------
        record_id : str
            FASTA record ID
        segment : str
            Segment type
        chain : str
            Chain type

        Returns
        -------
        str
            Standardized gene name
        """
        # Remove description after pipe
        name = record_id.split("|")[0].strip()

        # Check if already properly formatted
        if name.startswith(f"IG{chain}{segment}"):
            return name

        # Check if missing chain
        if name.startswith(f"IG{segment}"):
            return f"IG{chain}" + name[2:]

        # Completely custom name - make IMGT-like
        if not name.startswith("IG"):
            return f"IG{chain}{segment}-CUSTOM-{name}*01"

        return name

    def fetch_gene_by_name(
        self,
        name: str,
        species: str
    ) -> Optional[GermlineGene]:
        """
        Fetch specific custom gene by name.

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
        # Try each segment/chain combination
        for segment in ["V", "D", "J"]:
            for chain in ["H", "K", "L"]:
                genes = self.fetch_genes(species, segment, chain)
                for gene in genes:
                    if gene.name == name:
                        return gene

        return None

    def is_available(self, species: str) -> bool:
        """
        Check if custom sequences exist for species.

        Parameters
        ----------
        species : str
            Species name

        Returns
        -------
        bool
            True if custom data exists
        """
        species_dir = self.data_dir / species

        if not species_dir.exists():
            return False

        # Check if any FASTA files exist
        fasta_files = list(species_dir.glob("IG*.fasta"))
        return len(fasta_files) > 0

    def get_metadata(self) -> ProviderMetadata:
        """
        Get metadata about custom provider.

        Returns
        -------
        ProviderMetadata
            Provider metadata
        """
        # Find all species with custom data
        species = []
        if self.data_dir.exists():
            for species_dir in self.data_dir.iterdir():
                if species_dir.is_dir() and not species_dir.name.startswith("_"):
                    if list(species_dir.glob("IG*.fasta")):
                        species.append(species_dir.name)

        return ProviderMetadata(
            name="custom",
            version="user-supplied",
            last_updated=datetime.now(),
            species_available=species,
            url=None,
        )

    def download(self, species: List[str]) -> None:
        """
        Not applicable for custom provider.

        Raises
        ------
        NotImplementedError
            Custom provider doesn't support downloads
        """
        raise NotImplementedError(
            "Custom provider does not support download. "
            "Add FASTA files manually to sources/custom/<species>/ directory."
        )
