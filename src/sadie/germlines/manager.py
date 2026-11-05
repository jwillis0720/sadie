"""
Germline Manager - Priority-Based Database Access
==================================================

Manages multiple germline databases with priority-based lookup.
Priority works like Python dict: first provider wins on conflicts.

Design Principles (Zen of Python):
- Explicit is better than implicit (clear priority order)
- Simple is better than complex (dict-like priority)
- Flat is better than nested (no deep hierarchies)
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

from .models import GermlineGene
from .providers.base import GermlineProvider

logger = logging.getLogger(__name__)


class GermlineManager:
    """
    Manages multiple germline databases with priority-based lookup.

    Default priority: vdjbase > ogrdb > imgt > custom
    - VDJbase: Best for human/macaque (curated, validated alleles from population studies)
    - OGRDB: Good for mouse (community-curated novel alleles)
    - IMGT: Species diversity (comprehensive reference database)
    - Custom: Fill gaps (internal lab sequences for edge cases)

    Deduplication Logic:
    1. If gene names match exactly → use first provider's version
    2. If sequences match exactly → use first provider's version
    3. Otherwise keep both

    Examples
    --------
    >>> # Default priority: vdjbase, ogrdb, imgt, custom
    >>> manager = GermlineManager()
    >>> genes = manager.get_genes("human", "V", "H")
    >>>
    >>> # Custom priority: IMGT only
    >>> manager = GermlineManager(providers=["imgt"])
    >>> genes = manager.get_genes("human", "V", "H")
    """

    # Default provider priority order (highest to lowest):
    # 1. VDJbase: Best for human/macaque - curated, validated alleles from population studies
    # 2. OGRDB: Good for mouse - community-curated novel alleles
    # 3. IMGT: Species diversity - comprehensive reference database
    # 4. Custom: Fill gaps - internal lab sequences for edge cases
    DEFAULT_PROVIDERS = ["vdjbase", "ogrdb", "imgt", "custom"]

    def __init__(self, providers: Optional[List[str]] = None, data_dir: Optional[Path] = None):
        """
        Initialize manager with ordered list of providers.

        Parameters
        ----------
        providers : List[str], optional
            Ordered list of provider names.
            Default: ["vdjbase", "ogrdb", "imgt", "custom"]
            First provider has highest priority.
        data_dir : Path, optional
            Base directory for germline data.
            Default: module's sources/ directory
        """
        if providers is None:
            providers = self.DEFAULT_PROVIDERS

        if data_dir is None:
            data_dir = Path(__file__).parent / "sources"

        self.provider_names = providers
        self.data_dir = data_dir
        self.providers = self._initialize_providers(providers, data_dir)

    def _initialize_providers(self, provider_names: List[str], data_dir: Path) -> List[GermlineProvider]:
        """
        Initialize provider instances.

        Parameters
        ----------
        provider_names : List[str]
            Names of providers to initialize
        data_dir : Path
            Base data directory

        Returns
        -------
        List[GermlineProvider]
            Initialized provider instances
        """
        providers = []

        for name in provider_names:
            provider = self._create_provider(name, data_dir)
            if provider:
                providers.append(provider)
            else:
                logger.warning(f"Unknown provider: {name}")

        return providers

    def _create_provider(self, name: str, data_dir: Path) -> Optional[GermlineProvider]:
        """
        Create provider instance by name.

        Parameters
        ----------
        name : str
            Provider name
        data_dir : Path
            Base data directory

        Returns
        -------
        GermlineProvider or None
            Provider instance if known, None otherwise
        """
        if name == "custom":
            from .providers.custom import CustomProvider

            return CustomProvider(data_dir / "custom")

        if name == "imgt":
            from .providers.imgt import IMGTProvider

            return IMGTProvider(data_dir / "imgt")

        if name == "ogrdb":
            from .providers.ogrdb import OGRDBProvider

            return OGRDBProvider(data_dir / "ogrdb")

        if name == "vdjbase":
            from .providers.vdjbase import VDJbaseProvider

            return VDJbaseProvider(data_dir / "vdjbase")

        return None

    def get_genes(
        self, species: str, segment: str, chain: str, functional_only: bool = True, strict: bool = False
    ) -> List[GermlineGene]:
        """
        Get genes from all providers with priority-based deduplication.

        Single provider configuration applies to all V/D/J segments (FR-014).
        Per-segment provider mixing is not supported.

        Deduplication rules:
        1. Same gene name → first provider wins
        2. Same exact sequence → first provider wins
        3. Novel gene → include from any provider

        Parameters
        ----------
        species : str
            Species name (e.g., "human", "mouse")
        segment : str
            Segment type: "V", "D", or "J"
        chain : str
            Chain type: "H", "K", or "L"
        functional_only : bool
            Only return functional genes (default: True)
        strict : bool
            If True, raise ValueError when no genes found (default: False)

        Returns
        -------
        List[GermlineGene]
            Deduplicated genes in priority order

        Raises
        ------
        ValueError
            If strict=True and no providers have data for the species/segment/chain

        Examples
        --------
        >>> manager = GermlineManager()
        >>> genes = manager.get_genes("human", "V", "H")
        >>> print(f"Found {len(genes)} genes")
        """
        all_genes: Dict[str, GermlineGene] = {}
        seq_to_gene: Dict[str, str] = {}

        # Iterate providers in priority order
        for provider in self.providers:
            genes = self._fetch_from_provider(provider, species, segment, chain, functional_only)

            for gene in genes:
                if self._should_include_gene(gene, all_genes, seq_to_gene):
                    all_genes[gene.name] = gene
                    seq_to_gene[gene.sequence] = gene.name

        result = list(all_genes.values())

        if strict and len(result) == 0:
            available_species = self.get_available_species()
            raise ValueError(
                f"No germline data found for species='{species}', "
                f"segment='{segment}', chain='{chain}'. "
                f"Available species: {available_species}. "
                f"Run update_databases('{species}') to populate data."
            )

        return result

    def _fetch_from_provider(
        self, provider: GermlineProvider, species: str, segment: str, chain: str, functional_only: bool
    ) -> List[GermlineGene]:
        """
        Fetch genes from single provider with error handling.

        Parameters
        ----------
        provider : GermlineProvider
            Provider to fetch from
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type
        functional_only : bool
            Filter for functional genes

        Returns
        -------
        List[GermlineGene]
            Genes from provider (empty list on error)
        """
        try:
            genes = provider.fetch_genes(species, segment, chain)

            if functional_only:
                genes = [g for g in genes if g.is_functional]

            logger.debug(f"Provider {provider.name}: {len(genes)} genes for " f"{species} {chain}{segment}")

            return genes

        except Exception as e:
            logger.warning(f"Provider {provider.name} failed: {e}")
            return []

    def _should_include_gene(
        self, gene: GermlineGene, all_genes: Dict[str, GermlineGene], seq_to_gene: Dict[str, str]
    ) -> bool:
        """
        Determine if gene should be included based on deduplication rules.

        Parameters
        ----------
        gene : GermlineGene
            Gene to check
        all_genes : Dict[str, GermlineGene]
            Already included genes by name
        seq_to_gene : Dict[str, str]
            Sequence to gene name mapping

        Returns
        -------
        bool
            True if gene should be included
        """
        # Rule 1: Gene name conflict - first provider wins
        if gene.name in all_genes:
            logger.debug(f"Skipping {gene.name}: name conflict")
            return False

        # Rule 2: Sequence exact match - first provider wins
        if gene.sequence in seq_to_gene:
            existing_name = seq_to_gene[gene.sequence]
            logger.debug(f"Skipping {gene.name}: sequence matches {existing_name}")
            return False

        # Rule 3: Novel gene - include it
        return True

    def get_gene_by_name(self, name: str, species: str) -> Optional[GermlineGene]:
        """
        Get specific gene by name (first provider wins).

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
        for provider in self.providers:
            try:
                gene = provider.fetch_gene_by_name(name, species)
                if gene:
                    return gene
            except Exception as e:
                logger.debug(f"Provider {provider.name} lookup failed: {e}")
                continue

        return None

    def get_available_species(self) -> List[str]:
        """
        Get list of species with available data across all providers.

        Returns
        -------
        List[str]
            Unique species names
        """
        species_set = set()

        for provider in self.providers:
            try:
                metadata = provider.get_metadata()
                species_set.update(metadata.species_available)
            except Exception as e:
                logger.debug(f"Could not get species from {provider.name}: {e}")

        return sorted(species_set)

    def get_provider(self, name: str) -> GermlineProvider:
        """
        Get a configured provider by name.

        Parameters
        ----------
        name : str
            Provider name (e.g., "imgt", "ogrdb")

        Returns
        -------
        GermlineProvider
            Provider instance

        Raises
        ------
        ValueError
            If provider is not configured
        """
        for provider in self.providers:
            if provider.name == name:
                return provider

        raise ValueError(f"Provider '{name}' is not configured. Available providers: {self.provider_names}")

    def validate_species(self, name: str, species: str) -> None:
        """
        Validate provider availability for a species.

        Parameters
        ----------
        name : str
            Provider name (e.g., "imgt", "vdjbase")
        species : str
            Species name

        Raises
        ------
        ValueError
            If provider has no data for the species
        """
        provider = self.get_provider(name)

        if provider.is_available(species):
            return

        try:
            available_species = sorted(provider.get_metadata().species_available)
        except Exception as e:
            logger.debug(f"Could not get metadata for {provider.name}: {e}")
            available_species = []

        raise ValueError(
            f"Provider '{provider.name}' has no data for species '{species}'. "
            f"Available species: {available_species}."
        )
