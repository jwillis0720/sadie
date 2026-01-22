"""
Sadie Germlines Module
======================

Self-contained germline database management for immunoglobulin repertoires.

This module can be extracted as a standalone package for use in other projects.

Basic Usage:
    >>> from sadie.germlines import get_germline_genes, GermlineManager
    >>>
    >>> # Simple API - uses default priority (custom > imgt > ogrdb > vdjbase)
    >>> genes = get_germline_genes("human", "V", "H")
    >>>
    >>> # Advanced API with custom priority
    >>> manager = GermlineManager(providers=["custom", "imgt"])
    >>> genes = manager.get_genes("human", "V", "H")

Priority Logic:
    - Multiple databases are used by default (custom, IMGT, OGRDB, VDJbase)
    - First database in list has priority for conflicts
    - Conflicts resolved by: (1) gene name, (2) exact sequence match
    - Novel genes from any source are included

Architecture:
    - sources/ : Raw data from IMGT, OGRDB, custom user files
    - normalized/ : Merged and processed sequences (gapped + ungapped)
    - igblast/ : IgBLAST-ready databases and auxiliary files

Pipeline:
    sources/ → normalize → normalized/ → build → igblast/
"""

import logging
from typing import List, Optional

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _log_event(level: int, event: str, **kwargs) -> None:
    """
    Log event with structured key-value fields.

    Parameters
    ----------
    level : int
        Logging level
    event : str
        Event name
    **kwargs
        Key-value pairs to log
    """
    fields = " ".join(f"{k}={v}" for k, v in kwargs.items())
    logger.log(level, f"event={event} {fields}")

from .models import GermlineGene, ProviderMetadata
from .manager import GermlineManager
from .pipeline import GermlinePipeline
from .providers.vdjbase import VDJbaseProvider

# Global instances
_default_manager: Optional[GermlineManager] = None
_default_pipeline: Optional[GermlinePipeline] = None


def get_manager() -> GermlineManager:
    """Get or create default GermlineManager with default priority."""
    global _default_manager
    if _default_manager is None:
        _default_manager = GermlineManager()
    return _default_manager


def get_pipeline() -> GermlinePipeline:
    """Get or create default GermlinePipeline."""
    global _default_pipeline
    if _default_pipeline is None:
        from pathlib import Path
        base_dir = Path(__file__).parent
        _default_pipeline = GermlinePipeline(base_dir)
    return _default_pipeline


def get_germline_genes(
    species: str,
    segment: str,
    chain: str,
    providers: Optional[List[str]] = None,
    functional_only: bool = True
) -> List[GermlineGene]:
    """
    Get germline genes from all available databases.

    Parameters
    ----------
    species : str
        Species name (e.g., "human", "mouse")
    segment : str
        Segment type: "V", "D", or "J"
    chain : str
        Chain type: "H", "K", or "L"
    providers : List[str], optional
        Custom provider priority order. Default: ["custom", "imgt", "ogrdb", "vdjbase"]
    functional_only : bool
        Only return functional genes (default: True)

    Returns
    -------
    List[GermlineGene]
        Deduplicated genes (first provider wins on conflicts)

    Examples
    --------
    >>> # Get human heavy chain V genes
    >>> genes = get_germline_genes("human", "V", "H")
    >>>
    >>> # Custom priority: IMGT only
    >>> genes = get_germline_genes("human", "V", "H", providers=["imgt"])
    """
    if providers:
        manager = GermlineManager(providers=providers)
    else:
        manager = get_manager()

    return manager.get_genes(species, segment, chain, functional_only)


def get_gene_by_name(
    name: str,
    species: str,
    providers: Optional[List[str]] = None
) -> Optional[GermlineGene]:
    """
    Get specific gene by name (first provider that has it wins).

    Parameters
    ----------
    name : str
        Gene name (e.g., "IGHV1-69*01")
    species : str
        Species name
    providers : List[str], optional
        Custom provider priority order

    Returns
    -------
    GermlineGene or None
        Gene if found, None otherwise
    """
    if providers:
        manager = GermlineManager(providers=providers)
    else:
        manager = get_manager()

    return manager.get_gene_by_name(name, species)


def get_germlines_base_dir() -> "Path":
    """Get the base directory for the germlines module."""
    from pathlib import Path
    return Path(__file__).parent


def update_databases(species: str = "human", force: bool = False) -> None:
    """
    Update germline databases for species.

    Automatically detects changes and rebuilds as needed:
    - If sources/ changed → rebuild normalized/
    - If normalized/ changed → rebuild igblast/

    Parameters
    ----------
    species : str
        Species to update (default: "human")
    force : bool
        Force rebuild even if no changes detected
    """
    pipeline = get_pipeline()

    if force:
        pipeline.force_rebuild(species)
    else:
        pipeline.update(species)


__all__ = [
    # Core classes
    "GermlineManager",
    "GermlinePipeline",
    "GermlineGene",
    "ProviderMetadata",
    "VDJbaseProvider",

    # Public API functions
    "get_germline_genes",
    "get_gene_by_name",
    "get_manager",
    "get_pipeline",
    "get_germlines_base_dir",
    "update_databases",
]

__version__ = "1.0.0"
