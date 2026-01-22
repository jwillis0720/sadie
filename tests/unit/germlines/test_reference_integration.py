"""Integration tests for Reference system with germlines module.

Tests verify that Reference system works with germlines backend and produces
G3 API-compatible output format.
"""

import pytest

from sadie.reference import Reference
from sadie.reference.models import GeneEntry


class TestReferenceIntegration:
    """Test Reference system with germlines module backend."""

    def test_reference_with_germlines_backend(self, monkeypatch):
        """T027: Test Reference system using germlines backend.

        Verifies that:
        - Reference initializes with use_germlines=True
        - Gene queries succeed using germlines module
        - Output format matches G3 API structure
        """
        # Enable germlines backend for Reference
        ref = Reference(use_germlines=True)

        # Verify germlines components initialized
        assert hasattr(ref, "use_germlines"), "Should have use_germlines attribute"
        assert ref.use_germlines is True, "Should be using germlines backend"
        assert hasattr(ref, "germline_manager"), "Should have germline manager"
        assert hasattr(ref, "g3_adapter"), "Should have G3 adapter"

        # Add a gene using germlines backend
        gene_dict = {
            "species": "human",
            "gene": "IGHV1-69*01",
            "source": "imgt"
        }
        ref.add_gene(gene_dict)

        # Verify gene was added
        assert len(ref.data) == 1, "Should have one gene in data"

        # Get dataframe and verify structure
        df = ref.get_dataframe()
        assert not df.empty, "DataFrame should not be empty"
        assert len(df) == 1, "Should have one row"

        # Verify G3 API-compatible fields exist
        expected_fields = ["source", "common", "gene", "sequence"]
        for field in expected_fields:
            assert field in df.columns, f"Should have G3 field: {field}"

        # Verify gene name matches
        assert df["gene"].iloc[0] == "IGHV1-69*01", "Gene name should match"

    def test_reference_with_g3_backend(self, monkeypatch):
        """T028: Test Reference system with G3 backend (default).

        Verifies that:
        - Reference works with use_germlines=False (default)
        - Backwards compatibility maintained
        - G3 API endpoint is used
        """
        # Use default G3 backend
        ref = Reference(use_germlines=False)

        # Verify G3 backend is used
        assert hasattr(ref, "use_germlines"), "Should have use_germlines attribute"
        assert ref.use_germlines is False, "Should be using G3 backend"
        assert hasattr(ref, "endpoint"), "Should have G3 endpoint"
        assert not hasattr(ref, "germline_manager"), "Should not have germline manager"

        # This test validates that the default behavior is preserved
        # Note: Actual G3 API calls would require network access

    def test_output_format_consistency(self, monkeypatch):
        """T029: Verify output format consistency between backends.

        Verifies that:
        - Germlines backend output matches G3 API structure
        - All required G3 fields are present
        - Data types are consistent
        """
        # Test with germlines backend
        ref_germlines = Reference(use_germlines=True)

        gene_dict = {
            "species": "human",
            "gene": "IGHV1-69*01",
            "source": "imgt"
        }
        ref_germlines.add_gene(gene_dict)

        df_germlines = ref_germlines.get_dataframe()

        # Verify G3-compatible structure
        required_fields = [
            "source",
            "common",
            "gene",
            "label",
            "gene_segment",
            "receptor",
            "sequence",
        ]

        for field in required_fields:
            assert field in df_germlines.columns, (
                f"Germlines output missing required G3 field: {field}"
            )

        # Verify nested IMGT structure exists
        imgt_fields = [
            "imgt.sequence_gapped",
            "imgt.sequence_gapped_aa",
            "imgt.imgt_functional",
        ]

        for field in imgt_fields:
            assert field in df_germlines.columns, (
                f"Germlines output missing IMGT field: {field}"
            )

        # Verify data types
        assert isinstance(df_germlines["gene"].iloc[0], str), "Gene should be string"
        assert isinstance(df_germlines["sequence"].iloc[0], str), "Sequence should be string"

    def test_add_genes_batch(self, monkeypatch):
        """Test adding multiple genes with germlines backend."""
        ref = Reference(use_germlines=True)

        # Add multiple genes
        genes = ["IGHV1-69*01", "IGHV1-2*02", "IGHJ4*01"]
        ref.add_genes(species="human", source="imgt", genes=genes)

        # Verify genes were added
        df = ref.get_dataframe()
        assert len(df) >= 1, "Should have at least one gene (some may not be found)"

        # Verify genes are in dataframe
        added_genes = df["gene"].tolist()
        for gene in added_genes:
            assert gene in genes, f"Unexpected gene: {gene}"

    def test_gene_not_found_error(self, monkeypatch):
        """Test error handling when gene not found."""
        ref = Reference(use_germlines=True)

        # Try to add non-existent gene
        from sadie.reference.reference import G3Error

        with pytest.raises(G3Error):
            ref.add_gene({
                "species": "human",
                "gene": "IGHV999-999*99",  # Non-existent gene
                "source": "imgt"
            })


class TestGermlineToG3Adapter:
    """Test G3 adapter functionality."""

    def test_adapter_initialization(self, monkeypatch):
        """Test G3 adapter can be initialized."""
        from sadie.germlines.g3_adapter import GermlineToG3Adapter

        adapter = GermlineToG3Adapter()
        assert adapter is not None, "Adapter should initialize"

    def test_adapter_transform(self, monkeypatch):
        """Test adapter transforms GermlineGene to G3 format."""
        from sadie.germlines import get_gene_by_name
        from sadie.germlines.g3_adapter import GermlineToG3Adapter

        # Get a gene from germlines
        gene = get_gene_by_name("IGHV1-69*01", "human")
        assert gene is not None, "Should find gene in germlines"

        # Transform to G3 format
        adapter = GermlineToG3Adapter()
        g3_dict = adapter.to_g3_format(gene)

        # Verify G3 structure
        assert "gene" in g3_dict, "Should have gene field"
        assert "source" in g3_dict, "Should have source field"
        assert "sequence" in g3_dict, "Should have sequence field"
        assert "imgt" in g3_dict, "Should have imgt nested structure"

        # Verify IMGT fields
        assert "sequence_gapped" in g3_dict["imgt"], "Should have gapped sequence"
        assert "imgt_functional" in g3_dict["imgt"], "Should have functionality"

    def test_adapter_batch_transform(self, monkeypatch):
        """Test adapter batch transformation."""
        from sadie.germlines import get_germline_genes
        from sadie.germlines.g3_adapter import GermlineToG3Adapter

        # Get multiple genes
        genes = get_germline_genes("human", "V", "H")
        assert len(genes) > 0, "Should find genes"

        # Transform batch
        adapter = GermlineToG3Adapter()
        g3_list = adapter.to_g3_format_batch(genes[:5])  # First 5 genes

        assert len(g3_list) == 5, "Should transform all genes"
        for g3_dict in g3_list:
            assert "gene" in g3_dict, "Each should have gene field"
            assert "imgt" in g3_dict, "Each should have IMGT structure"
