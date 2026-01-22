"""
Tests for GermlineManager
==========================

Tests priority-based germline database management.
"""

import pytest
from pathlib import Path

# TODO: Implement tests after data is populated


def test_manager_initialization():
    """Test manager can be initialized."""
    from sadie.germlines.manager import GermlineManager

    manager = GermlineManager(providers=["imgt"])
    assert manager is not None
    assert len(manager.providers) == 1


def test_manager_default_providers():
    """Test default provider order."""
    from sadie.germlines.manager import GermlineManager

    manager = GermlineManager()
    assert manager.provider_names == ["custom", "imgt", "ogrdb", "vdjbase"]


def test_manager_custom_providers():
    """Test custom provider order."""
    from sadie.germlines.manager import GermlineManager

    manager = GermlineManager(providers=["ogrdb", "imgt"])
    assert manager.provider_names == ["ogrdb", "imgt"]


# TODO: Add more tests
# - test_fetch_genes()
# - test_priority_deduplication()
# - test_novel_genes_included()
# - test_get_gene_by_name()
