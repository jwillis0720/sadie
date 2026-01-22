"""
Base Provider Interface
=======================

Abstract base class defining the interface for germline data providers.

Design Principles:
- Single responsibility: fetch germline sequences
- Open/closed: easy to add new providers
- Dependency inversion: depend on abstractions
"""

import logging
from abc import ABC, abstractmethod
from typing import List, Optional
from pathlib import Path

from ..models import GermlineGene, ProviderMetadata

logger = logging.getLogger(__name__)


ERROR_TEMPLATES = {
    "parse_failed": "provider={provider} action=parse_failed path={path} error={error}",
    "gene_creation_failed": "provider={provider} action=gene_creation_failed gene={gene} error={error}",
    "http_error": "provider={provider} action=http_error url={url} status={status}",
    "data_not_found": "provider={provider} action=data_not_found species={species}",
    "validation_failed": "provider={provider} action=validation_failed gene={gene} reason={reason}",
}


class GermlineProvider(ABC):
    """
    Abstract base class for germline database providers.

    All providers must implement these methods to provide
    germline sequence data from different sources.

    Attributes
    ----------
    name : str
        Provider name (e.g., "imgt", "ogrdb", "custom")
    data_dir : Path
        Base directory for this provider's data
    """

    def __init__(self, data_dir: Optional[Path] = None):
        """
        Initialize provider.

        Parameters
        ----------
        data_dir : Path, optional
            Base directory for provider data
        """
        if data_dir is None:
            # Default to sources/<provider_name>
            data_dir = Path(__file__).parent.parent / "sources"

        self.data_dir = data_dir
        self.name = self._get_provider_name()

    def _get_provider_name(self) -> str:
        """
        Get provider name from class name.

        Returns
        -------
        str
            Provider name in lowercase
        """
        class_name = self.__class__.__name__
        return class_name.replace("Provider", "").lower()

    @abstractmethod
    def fetch_genes(
        self,
        species: str,
        segment: str,
        chain: str
    ) -> List[GermlineGene]:
        """
        Fetch genes for given species/segment/chain.

        Parameters
        ----------
        species : str
            Species name (e.g., "human", "mouse")
        segment : str
            Segment type: "V", "D", or "J"
        chain : str
            Chain type: "H", "K", or "L"

        Returns
        -------
        List[GermlineGene]
            List of germline genes
        """
        pass

    @abstractmethod
    def fetch_gene_by_name(
        self,
        name: str,
        species: str
    ) -> Optional[GermlineGene]:
        """
        Fetch specific gene by name.

        Parameters
        ----------
        name : str
            Gene name (e.g., "IGHV1-69*01")
        species : str
            Species name

        Returns
        -------
        GermlineGene or None
            Gene if found, None otherwise
        """
        pass

    @abstractmethod
    def get_metadata(self) -> ProviderMetadata:
        """
        Get provider metadata (version, species, etc.).

        Returns
        -------
        ProviderMetadata
            Provider metadata
        """
        pass

    @abstractmethod
    def is_available(self, species: str) -> bool:
        """
        Check if data is available for species.

        Parameters
        ----------
        species : str
            Species name

        Returns
        -------
        bool
            True if data available
        """
        pass

    @abstractmethod
    def download(self, species: List[str]) -> None:
        """
        Download data for specified species.

        Parameters
        ----------
        species : List[str]
            Species to download

        Raises
        ------
        NotImplementedError
            If provider doesn't support downloads
        """
        pass

    def get_fasta_path(
        self,
        species: str,
        segment: str,
        chain: str
    ) -> Path:
        """
        Get path to FASTA file for segment.

        This is a helper method that provides a standard path
        structure. Providers can override if needed.

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
            Path to FASTA file
        """
        filename = f"IG{chain}{segment}.fasta"
        return self.data_dir / species / filename
