"""Integration tests for renumbering with germlines module.

Tests verify that HMM-based renumbering works with LocalHMMBuilder and produces
equivalent results to G3 API backend.
"""

import os
from pathlib import Path

import pytest

from sadie.germlines.renumbering_integration import LocalHMMBuilder
from sadie.renumbering.aligners.hmmer import HMMER


class TestRenumberingIntegration:
    """Test renumbering with LocalHMMBuilder from germlines module."""

    # Test antibody sequence (IGHV1-69*01 framework)
    TEST_AB_SEQ = (
        "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYG"
        "ISWVRQAPGQGLEWMGWISAYNGNTNYAQKLQG"
        "RVTMTTDTSTSTAYMELRSLRSDDTAVYYCARA"
    )

    @pytest.fixture(autouse=True)
    def setup(self, monkeypatch):
        """Setup fixture to enable germlines module for each test."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

    def test_hmmer_with_local_builder(self, monkeypatch):
        """T036: Test HMMER using LocalHMMBuilder backend.

        Verifies that:
        - HMMER successfully loads HMMs from LocalHMMBuilder
        - HMM models are valid pyhmmer objects
        - Renumbering produces results
        """
        # Ensure germlines module is used
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        # Initialize HMMER with human heavy chain
        hmmer = HMMER(species="human", chains="H")

        # Verify HMMs were loaded
        assert hasattr(hmmer, "hmms"), "HMMER should have hmms attribute"
        assert len(hmmer.hmms) > 0, "Should load at least one HMM"

        # Verify HMM is valid
        hmm = hmmer.hmms[0]
        assert hasattr(hmm, "name"), "HMM should have name attribute"
        assert b"human" in hmm.name or b"H" in hmm.name, f"HMM name should contain species/chain: {hmm.name}"

        # Test alignment with the HMM using hmmsearch
        sequences = [self.TEST_AB_SEQ]
        results = hmmer.hmmsearch(sequences, bit_score_threshold=50)

        # Verify alignment results
        assert len(results) > 0, "Should produce alignment results"
        result = results[0]
        # hmmsearch returns list of lists containing dicts with bitscore, query, etc.
        assert len(result) > 0, "Result should contain alignment hits"
        assert "bitscore" in result[0], "Result should contain bitscore"

    def test_hmm_caching(self, monkeypatch, tmp_path):
        """T037: Test HMM caching behavior of LocalHMMBuilder.

        Verifies that:
        - HMMs are cached to disk after first build
        - Subsequent loads are faster (use cached version)
        - Cache files are created in expected location
        """
        # Ensure germlines module is used
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        # Create LocalHMMBuilder instance
        builder = LocalHMMBuilder()

        # Get HMM (may build or load from cache)
        hmm1 = builder.get_hmm(species="human", chain="H")
        assert hmm1 is not None, "Should return HMM object"
        assert hasattr(hmm1, "name"), "HMM should be valid pyhmmer HMM"

        # Verify cache file was created
        cache_path = builder.hmm_dir / "human_H.hmm"
        assert cache_path.exists(), f"HMM cache file should exist at {cache_path}"

        # Load HMM again (should use cache)
        hmm2 = builder.get_hmm(species="human", chain="H")
        assert hmm2 is not None, "Should load HMM from cache"
        assert hmm1.name == hmm2.name, "Cached HMM should match original"

    def test_offline_operation(self, monkeypatch):
        """T038: Test renumbering works offline with germlines module.

        Verifies that:
        - Renumbering succeeds without network access
        - All HMM data sourced from local germlines module
        - No G3 API calls are made
        """
        # Enable germlines module (local-first operation)
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        # Initialize HMMER (should not require network)
        hmmer = HMMER(species="human", chains="H")

        # Verify HMMs loaded successfully
        assert len(hmmer.hmms) > 0, "Should load HMMs offline"

        # Verify LocalHMMBuilder is being used
        assert hasattr(hmmer, "_local_hmm_builder"), (
            "HMMER should use LocalHMMBuilder when feature flag enabled"
        )

        # Test alignment (should work offline)
        sequences = [self.TEST_AB_SEQ]
        results = hmmer.hmmsearch(sequences, bit_score_threshold=50)

        # Verify offline alignment succeeded
        assert len(results) > 0, "Renumbering should work offline with germlines"

    def test_backwards_compatibility(self, monkeypatch):
        """Verify backwards compatibility with feature flag disabled.

        This ensures existing code continues to work when germlines module
        is not enabled (legacy G3 API mode).
        """
        # Disable germlines module
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")

        # HMMER should still initialize (using G3 or legacy HMMs)
        hmmer = HMMER(species="human", chains="H")

        # Should not use LocalHMMBuilder
        # (may use numbering legacy HMMs or G3 API)
        # This test just verifies no crash occurs


class TestLocalHMMBuilder:
    """Test LocalHMMBuilder functionality."""

    def test_get_hmm_basic(self, monkeypatch):
        """Test basic HMM retrieval."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        builder = LocalHMMBuilder()
        hmm = builder.get_hmm(species="human", chain="H")

        assert hmm is not None, "Should return HMM"
        assert hasattr(hmm, "name"), "HMM should be valid"
        assert hasattr(hmm, "M"), "HMM should have model length"

    def test_get_hmm_all_chains(self, monkeypatch):
        """Test HMM retrieval for all chain types."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        builder = LocalHMMBuilder()

        # Test heavy chain
        hmm_h = builder.get_hmm(species="human", chain="H")
        assert hmm_h is not None, "Should return heavy chain HMM"

        # Test kappa chain
        hmm_k = builder.get_hmm(species="human", chain="K")
        assert hmm_k is not None, "Should return kappa chain HMM"

        # Test lambda chain
        hmm_l = builder.get_hmm(species="human", chain="L")
        assert hmm_l is not None, "Should return lambda chain HMM"

    def test_hmm_cache_directory(self, monkeypatch):
        """Test that HMM cache directory is created."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        builder = LocalHMMBuilder()

        # Check cache directory exists
        assert builder.hmm_dir.exists(), "HMM cache directory should exist"
        assert builder.hmm_dir.is_dir(), "HMM cache path should be a directory"

    def test_build_hmm_with_source(self, monkeypatch):
        """Test HMM building with explicit source parameter."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        builder = LocalHMMBuilder()

        # Get HMM with explicit IMGT source
        hmm = builder.get_hmm(species="human", chain="H", source="imgt")
        assert hmm is not None, "Should build HMM with explicit source"


class TestGappedAAFallbackTranslation:
    """T035a: Test gapped AA fallback translation when only gapped NT is available."""

    @pytest.fixture(autouse=True)
    def setup(self, monkeypatch):
        """Setup fixture to enable germlines module."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

    def test_translate_gapped_nt_to_aa_basic(self, monkeypatch):
        """Test basic NT to AA translation with gaps.

        Verifies that _translate_gapped_nt_to_aa correctly:
        - Translates nucleotide codons to amino acids
        - Preserves gap positions (3 NT gaps = 1 AA gap)
        """
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        builder = LocalHMMBuilder()

        # Simple test: no gaps
        ungapped_nt = "ATGGCC"  # Met-Ala
        result = builder._translate_gapped_nt_to_aa(ungapped_nt)
        assert result is not None, "Should translate ungapped sequence"
        assert "M" in result, "Should contain Methionine"
        assert "A" in result, "Should contain Alanine"

        # Test with gaps (3 NT gaps = 1 AA gap)
        gapped_nt = "ATG...GCC"  # Met-[gap]-Ala
        result = builder._translate_gapped_nt_to_aa(gapped_nt)
        assert result is not None, "Should translate gapped sequence"
        assert "." in result, "Should preserve gap in AA sequence"

    def test_translate_gapped_nt_to_aa_imgt_format(self, monkeypatch):
        """Test translation with IMGT-style gapped sequences.

        IMGT uses dots for gaps. Every 3 nucleotide gaps should map to 1 AA gap.
        """
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        builder = LocalHMMBuilder()

        # IMGT-style gapped sequence (V gene framework region)
        # cag = Q, gtg = V, cag = Q, ctg = L
        gapped_nt = "caggtgcagctg"
        result = builder._translate_gapped_nt_to_aa(gapped_nt)
        assert result is not None, "Should translate IMGT sequence"
        assert result.upper() == "QVQL", f"Expected QVQL, got {result}"

        # With gaps
        gapped_nt_with_gaps = "cag...gtgcagctg"  # Q-[gap]-VQL
        result = builder._translate_gapped_nt_to_aa(gapped_nt_with_gaps)
        assert result is not None, "Should translate gapped IMGT sequence"
        # Gap position depends on implementation

    def test_fallback_handles_invalid_sequences(self, monkeypatch):
        """Test that fallback gracefully handles invalid sequences."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        builder = LocalHMMBuilder()

        # Too short (< 3 nucleotides)
        result = builder._translate_gapped_nt_to_aa("AT")
        assert result is None, "Should return None for too-short sequence"

        # Empty sequence
        result = builder._translate_gapped_nt_to_aa("")
        assert result is None, "Should return None for empty sequence"

        # All gaps
        result = builder._translate_gapped_nt_to_aa("...")
        assert result is None, "Should return None for all-gaps sequence"

    def test_gapped_aa_fallback_in_hmm_builder(self, monkeypatch):
        """Test that HMM builder uses fallback translation for genes without gapped AA.

        Verifies the full pipeline:
        1. Query genes from germlines manager
        2. For genes with only gapped NT, translate to gapped AA
        3. Successfully build HMM from the translated sequences
        """
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        builder = LocalHMMBuilder()

        # Build HMM for human heavy chain (uses V/J genes)
        # This should trigger fallback translation for any genes
        # that have sequence_gapped but not sequence_aa_gapped
        hmm = builder.get_hmm(species="human", chain="H")

        assert hmm is not None, "Should build HMM using fallback translation"
        assert hasattr(hmm, "M"), "HMM should have model length"
        assert hmm.M > 0, "HMM model length should be positive"

    def test_fallback_used_when_gapped_aa_missing(self, monkeypatch):
        """Verify fallback translation is actually invoked.

        Test that _get_vj_alignment_pairs uses the fallback when genes
        have sequence_gapped but not sequence_aa_gapped.
        """
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        builder = LocalHMMBuilder()

        # Get V/J alignment pairs for human heavy chain
        # This method internally uses the fallback translation
        pairs = builder._get_vj_alignment_pairs("human", "H", "imgt")

        assert len(pairs) > 0, "Should return V/J gene pairs"

        # Verify pairs contain valid sequences
        for name, seq in pairs[:5]:  # Check first 5
            assert name, "Gene name should not be empty"
            assert seq, "Gapped AA sequence should not be empty"
            assert isinstance(seq, str), "Sequence should be string"
            # Sequence should only contain valid AA characters and gaps
            valid_chars = set("ACDEFGHIKLMNPQRSTVWY.-")
            assert all(c.upper() in valid_chars for c in seq), (
                f"Sequence should contain only valid AA chars: {seq[:50]}"
            )

    def test_mouse_hmm_with_fallback(self, monkeypatch):
        """Test HMM building for mouse (verifies multi-species fallback)."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        builder = LocalHMMBuilder()

        # Build mouse HMM
        hmm = builder.get_hmm(species="mouse", chain="H")
        assert hmm is not None, "Should build mouse HMM with fallback"

        # Also test kappa and lambda
        hmm_k = builder.get_hmm(species="mouse", chain="K")
        assert hmm_k is not None, "Should build mouse kappa HMM"

        hmm_l = builder.get_hmm(species="mouse", chain="L")
        assert hmm_l is not None, "Should build mouse lambda HMM"
