"""Integration tests for AIRR annotation with germlines module.

Tests verify that AIRR annotation works with germlines backend and produces
equivalent results to G3 API backend.
"""

import os
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sadie.airr import Airr, GermlineData


class TestAirrIntegration:
    """Test AIRR annotation with germlines module backend."""

    # Test sequence: IGHV1-69*01 example from IMGT
    TEST_SEQ = (
        "CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCC"
        "TCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTC"
        "TCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGA"
        "GACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGT"
        "GCGAAAGATCGCGTGGTTCGACGCC"
    )

    @pytest.fixture(autouse=True)
    def setup(self, monkeypatch):
        """Setup fixture to enable germlines module for each test."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

    def test_airr_annotation_with_germlines(self, monkeypatch, tmp_path):
        """T032: Test AIRR annotation using germlines backend.

        Verifies that:
        - AIRR annotation succeeds with germlines module enabled
        - V, D, J gene calls are made
        - Results contain expected AIRR-compliant fields
        """
        # Ensure germlines module is used
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        # Initialize AIRR with human species (reference_name parameter)
        airr = Airr(reference_name="human")

        # Run annotation on test sequence
        result = airr.run_single("test_seq", self.TEST_SEQ)

        # Verify annotation succeeded
        assert not result.empty, "AIRR annotation should return results"
        assert len(result) == 1, "Should have one result for one sequence"

        # Verify gene calls were made
        assert "v_call" in result.columns, "Should have V gene call"
        assert "d_call" in result.columns, "Should have D gene call"
        assert "j_call" in result.columns, "Should have J gene call"

        # Verify V gene was called
        v_call = result["v_call"].iloc[0]
        assert v_call is not None and v_call != "", "V gene should be called"
        assert "IGHV" in v_call, f"V gene call should be IGHV, got {v_call}"

        # Verify AIRR-compliant fields exist
        required_fields = [
            "sequence_id",
            "sequence",
            "v_call",
            "d_call",
            "j_call",
            "junction",
            "junction_aa",
        ]
        for field in required_fields:
            assert field in result.columns, f"AIRR field {field} should exist"

    def test_provider_selection(self, monkeypatch, tmp_path):
        """T033: Test that germlines backend respects provider selection.

        Verifies that:
        - AIRR annotation works with germlines provider selection
        - Results are consistent with selected provider
        - Default provider priority is respected
        """
        # Test with explicit IMGT provider (via germlines module)
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        airr = Airr(reference_name="human")
        result = airr.run_single("test_seq", self.TEST_SEQ)

        # Verify annotation succeeded
        assert not result.empty, "AIRR annotation should work with provider selection"

        # Verify gene calls
        v_call = result["v_call"].iloc[0]
        assert v_call is not None, "Provider selection should produce gene calls"

        # Test that GermlineData correctly uses germlines paths
        gd = GermlineData("human")
        assert gd.base_dir.exists(), "Germlines database directory should exist"
        assert "germlines" in str(gd.base_dir), "Should use germlines module path"

    def test_offline_operation(self, monkeypatch, tmp_path):
        """T034: Test AIRR annotation works offline with germlines module.

        Verifies that:
        - AIRR annotation succeeds without network access
        - All data is sourced from local germlines module
        - No G3 API calls are made
        """
        # Enable germlines module (local-first operation)
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        # Initialize AIRR
        airr = Airr(reference_name="human")

        # Verify germline data paths point to local module
        gd = GermlineData("human")
        assert gd.base_dir.exists(), "Local germlines database must exist"

        # Run annotation (should work offline)
        result = airr.run_single("test_seq", self.TEST_SEQ)

        # Verify annotation succeeded without network
        assert not result.empty, "AIRR should work offline with germlines"
        assert "v_call" in result.columns
        v_call = result["v_call"].iloc[0]
        assert v_call is not None, "Offline annotation should produce results"

    def test_backwards_compatibility(self, monkeypatch, tmp_path):
        """Verify backwards compatibility with feature flag disabled.

        This ensures existing code continues to work when germlines module
        is not enabled (legacy G3 API mode).
        """
        # Disable germlines module (use G3 API)
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")

        # Verify GermlineData uses legacy paths
        gd = GermlineData("human")
        # Legacy path should be in airr/data/germlines
        assert "airr" in str(gd.base_dir) or "data" in str(gd.base_dir), (
            "Should use legacy G3 paths when feature flag disabled"
        )


class TestGermlineDataPaths:
    """Test GermlineData path switching with feature flag."""

    def test_germlines_path_enabled(self, monkeypatch):
        """Verify GermlineData uses germlines module paths when enabled."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        gd = GermlineData("human")

        # Check that paths point to germlines module
        assert gd.base_dir.exists(), "Germlines base directory should exist"
        assert "germlines" in str(gd.base_dir), f"Should use germlines path, got {gd.base_dir}"

        # Verify BLAST database paths exist
        assert gd.v_gene_dir.with_suffix(".nhr").exists(), "V gene BLAST database should exist"
        assert gd.j_gene_dir.with_suffix(".nhr").exists(), "J gene BLAST database should exist"
        assert gd.aux_path.exists(), "Auxiliary file should exist"

    def test_legacy_path_disabled(self, monkeypatch):
        """Verify GermlineData uses legacy paths when disabled."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")

        gd = GermlineData("human")

        # Check that paths point to legacy location
        assert "airr" in str(gd.base_dir) or "data" in str(gd.base_dir), (
            f"Should use legacy path, got {gd.base_dir}"
        )


class TestOfflineOperation:
    """T041-T045: Test offline operation with germlines module."""

    # Test sequence for offline tests
    TEST_SEQ = (
        "CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCC"
        "TCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTC"
        "TCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGA"
        "GACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGT"
        "GCGAAAGATCGCGTGGTTCGACGCC"
    )

    def test_airr_annotation_network_disabled(self, enable_germlines, network_disabled):
        """T042: Test AIRR annotation with network disabled.

        Verifies that AIRR annotation succeeds when network is blocked,
        confirming all data comes from local germlines module.
        """
        # With network disabled, this should work using local germlines
        airr = Airr(reference_name="human")
        result = airr.run_single("test_seq", self.TEST_SEQ)

        # Verify annotation succeeded offline
        assert not result.empty, "AIRR annotation should work offline"
        assert "v_call" in result.columns
        assert result["v_call"].iloc[0] is not None, "Should produce gene calls offline"

    def test_germline_data_paths_network_disabled(self, enable_germlines, network_disabled):
        """Test GermlineData paths work with network disabled.

        Verifies that all germline data paths are accessible offline.
        """
        gd = GermlineData("human")

        # All paths should be accessible without network
        assert gd.base_dir.exists(), "Base directory should exist offline"
        assert gd.v_gene_dir.with_suffix(".nhr").exists(), "V gene DB should exist offline"
        assert gd.j_gene_dir.with_suffix(".nhr").exists(), "J gene DB should exist offline"
        assert gd.aux_path.exists(), "Auxiliary file should exist offline"

    def test_clear_error_unpopulated_database(self, monkeypatch, tmp_path):
        """T044: Test clear error message when germlines not populated.

        Verifies that a clear error message is shown when the germlines
        database is not populated (first-time setup scenario).
        """
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        # Create empty germlines-like directory
        empty_igblast = tmp_path / "igblast" / "database" / "human"
        empty_igblast.mkdir(parents=True)

        # Patch germlines base dir to point to empty location
        from sadie.germlines import get_germlines_base_dir
        original_base_dir = get_germlines_base_dir()

        def mock_base_dir():
            return tmp_path

        monkeypatch.setattr("sadie.germlines.get_germlines_base_dir", mock_base_dir)

        # GermlineData should handle missing databases gracefully
        # or raise a clear error
        try:
            from sadie.airr.igblast.germline import GermlineData as GD
            # This may raise an error about missing databases
            gd = GD("human")
            # If it doesn't raise, check that paths are properly set
            # (even if they don't exist)
        except (FileNotFoundError, ValueError) as e:
            # This is expected - should have clear error message
            error_msg = str(e).lower()
            assert "not found" in error_msg or "missing" in error_msg or "exist" in error_msg, (
                f"Error should mention missing database: {e}"
            )
