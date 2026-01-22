"""
Data Models for Germline Sequences
===================================

Unified data models for germline genes across all providers (IMGT, OGRDB, custom).

Design Principles:
- Explicit field names (no abbreviations)
- Optional fields clearly marked
- Validation through Pydantic
- Source tracking for provenance
"""

from typing import Optional, Dict, List, Tuple, Any
from pydantic import BaseModel, Field, field_validator
from datetime import datetime


class GermlineGene(BaseModel):
    """
    Unified germline gene model across all providers.

    This model represents a single germline gene sequence with all
    associated metadata. It normalizes data from IMGT, OGRDB, and
    custom sources into a common format.

    Attributes
    ----------
    name : str
        Gene name in IMGT format (e.g., "IGHV1-69*01")
    species : str
        Species name (e.g., "human", "mouse")
    segment : str
        Segment type: "V", "D", or "J"
    chain : str
        Chain type: "H" (Heavy), "K" (Kappa), or "L" (Lambda)
    sequence : str
        Ungapped nucleotide sequence (ACGT only)
    sequence_gapped : str, optional
        IMGT-gapped nucleotide sequence (dots for gaps)
    sequence_aa : str, optional
        Ungapped amino acid sequence
    sequence_aa_gapped : str, optional
        IMGT-gapped amino acid sequence
    is_functional : bool
        Whether gene is functional (default: True)
    functionality : str
        Functionality code: "F" (functional), "ORF" (open reading frame),
        "P" (pseudogene)
    regions : dict, optional
        Sequence regions (CDR1, CDR2, CDR3, FWR1-4) if available
    region_positions : dict, optional
        Start/end positions for each region
    source : str
        Data source: "imgt", "ogrdb", "custom"
    source_version : str, optional
        Version/date of source data
    allele : str, optional
        Allele designation
    gene_family : str, optional
        Gene family classification
    accession : str, optional
        GenBank/EMBL accession number
    """

    # Core identifiers
    name: str = Field(..., description="Gene name (e.g., IGHV1-69*01)")
    species: str = Field(..., description="Species (e.g., human)")
    segment: str = Field(..., description="Segment: V, D, or J")
    chain: str = Field(..., description="Chain: H, K, or L")

    # Sequences
    sequence: str = Field(..., description="Ungapped nucleotide sequence")
    sequence_gapped: Optional[str] = Field(None, description="IMGT-gapped nucleotide")
    sequence_aa: Optional[str] = Field(None, description="Ungapped amino acid")
    sequence_aa_gapped: Optional[str] = Field(None, description="IMGT-gapped amino acid")

    # Functional annotation
    is_functional: bool = Field(True, description="Is gene functional")
    functionality: str = Field("F", description="F, ORF, or P")

    # IMGT regions (if available)
    regions: Optional[Dict[str, str]] = Field(
        None,
        description="Sequence regions (CDR1, CDR2, CDR3, FWR1-4)"
    )
    region_positions: Optional[Dict[str, Tuple[int, int]]] = Field(
        None,
        description="Start/end positions for regions"
    )

    # Source tracking
    source: str = Field(..., description="Data source: imgt, ogrdb, custom")
    source_version: Optional[str] = Field(None, description="Source version/date")

    # Metadata
    allele: Optional[str] = Field(None, description="Allele designation")
    gene_family: Optional[str] = Field(None, description="Gene family")
    accession: Optional[str] = Field(None, description="GenBank/EMBL accession")

    @field_validator("segment")
    @classmethod
    def validate_segment(cls, v: str) -> str:
        """Validate segment is V, D, or J."""
        v = v.upper()
        if v not in ["V", "D", "J"]:
            raise ValueError(f"Segment must be V, D, or J, got: {v}")
        return v

    @field_validator("chain")
    @classmethod
    def validate_chain(cls, v: str) -> str:
        """Validate chain is H, K, or L."""
        v = v.upper()
        if v not in ["H", "K", "L"]:
            raise ValueError(f"Chain must be H, K, or L, got: {v}")
        return v

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        """Validate sequence contains only valid nucleotides."""
        v = v.upper()
        valid_chars = set("ACGTN")
        invalid = set(v) - valid_chars
        if invalid:
            raise ValueError(f"Sequence contains invalid characters: {invalid}")
        return v

    @field_validator("functionality")
    @classmethod
    def validate_functionality(cls, v: str) -> str:
        """Validate functionality is F, ORF, or P."""
        v = v.upper()
        if v not in ["F", "ORF", "P"]:
            raise ValueError(f"Functionality must be F, ORF, or P, got: {v}")
        return v

    def __str__(self) -> str:
        """String representation."""
        return f"{self.name} ({self.source})"

    def __repr__(self) -> str:
        """Detailed representation."""
        return (
            f"GermlineGene(name='{self.name}', "
            f"species='{self.species}', "
            f"segment='{self.segment}', "
            f"chain='{self.chain}', "
            f"source='{self.source}')"
        )


class ProviderMetadata(BaseModel):
    """
    Provider metadata for version tracking.

    Tracks information about each germline data provider
    (IMGT, OGRDB, custom) including version and available species.

    Attributes
    ----------
    name : str
        Provider name (e.g., "imgt", "ogrdb", "custom")
    version : str
        Version identifier or date
    last_updated : datetime
        When data was last updated
    species_available : List[str]
        Species with available data
    url : str, optional
        Source URL if applicable
    """

    name: str = Field(..., description="Provider name")
    version: str = Field(..., description="Version or date")
    last_updated: datetime = Field(..., description="Last update time")
    species_available: List[str] = Field(
        default_factory=list,
        description="Available species"
    )
    url: Optional[str] = Field(None, description="Source URL")

    def __str__(self) -> str:
        """String representation."""
        return f"{self.name} v{self.version} ({len(self.species_available)} species)"


class ProcessingMetadata(BaseModel):
    """
    Metadata about processed germline files.

    Tracks when files were processed and their content
    for change detection.

    Attributes
    ----------
    source_file : str
        Path to source FASTA file
    processed_at : datetime
        When processing occurred
    num_sequences : int
        Number of sequences processed
    file_hash : str
        Hash of source file for change detection
    sequences : List[dict]
        Summary of processed sequences
    """

    source_file: str = Field(..., description="Source file path")
    processed_at: datetime = Field(..., description="Processing timestamp")
    num_sequences: int = Field(..., description="Number of sequences")
    file_hash: str = Field(..., description="File hash for change detection")
    sequences: List[Dict[str, Any]] = Field(
        default_factory=list,
        description="Sequence summaries"
    )

    def __str__(self) -> str:
        """String representation."""
        return f"Processed {self.num_sequences} sequences at {self.processed_at}"
