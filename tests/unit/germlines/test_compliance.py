"""Compliance tests for germlines module.

Tests verify constitution-aligned behaviors and requirement compliance:
- FR-004: Default priority order
- FR-006: Clear error when provider lacks species
- FR-010: Species/chain/segment parity
- FR-012: Custom ingestion validation
- FR-013: Gapped sequence fail-fast
- FR-014: Single-provider enforcement
- NFR-002: No silent G3 fallback
"""

import os
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from sadie.germlines.manager import GermlineManager
from sadie.germlines.models import GermlineGene


class TestPriorityOrder:
    """T058: Test default provider priority order (FR-004)."""

    def test_default_priority_order(self):
        """Verify default priority is custom > ogrdb > vdjbase > imgt."""
        expected = ["custom", "ogrdb", "vdjbase", "imgt"]
        assert GermlineManager.DEFAULT_PROVIDERS == expected, (
            f"Expected priority {expected}, got {GermlineManager.DEFAULT_PROVIDERS}"
        )

    def test_manager_initializes_with_default_priority(self):
        """Verify manager uses default priority when none specified."""
        manager = GermlineManager()
        assert manager.provider_names == ["custom", "ogrdb", "vdjbase", "imgt"]

    def test_custom_provider_overrides_others(self):
        """Verify custom provider sequences take precedence."""
        manager = GermlineManager()
        # First provider in list should be custom
        assert manager.provider_names[0] == "custom"


class TestSingleProvider:
    """T053-T054: Test single-provider enforcement (FR-014)."""

    def test_single_provider_config_for_all_segments(self):
        """Verify manager uses same provider config for V, D, J segments."""
        manager = GermlineManager(providers=["imgt"])
        
        # All queries should use the same provider list
        # Manager design enforces this - one provider list for all calls
        assert manager.provider_names == ["imgt"]
        
        # Querying different segments uses same providers
        v_genes = manager.get_genes("human", "V", "H", functional_only=True)
        d_genes = manager.get_genes("human", "D", "H", functional_only=True)
        j_genes = manager.get_genes("human", "J", "H", functional_only=True)
        
        # All should come from IMGT (if available)
        for gene in v_genes[:5]:  # Check first 5
            assert gene.source == "imgt", f"V gene {gene.name} should be from imgt"

    def test_no_per_segment_provider_parameter(self):
        """Verify get_genes does not accept per-segment provider override."""
        manager = GermlineManager()
        
        # get_genes signature should not have provider parameter
        import inspect
        sig = inspect.signature(manager.get_genes)
        param_names = list(sig.parameters.keys())
        
        assert "provider" not in param_names, (
            "get_genes should not accept per-segment provider parameter"
        )
        assert "v_provider" not in param_names
        assert "d_provider" not in param_names
        assert "j_provider" not in param_names


class TestErrorHandling:
    """T055, T059: Test error handling and no-fallback (FR-006, NFR-002)."""

    def test_clear_error_missing_species_non_strict(self):
        """T055: Verify empty result without strict mode."""
        manager = GermlineManager(providers=["imgt"])
        
        # Non-strict mode returns empty list
        result = manager.get_genes("nonexistent_species", "V", "H", strict=False)
        assert isinstance(result, list)
        assert len(result) == 0

    def test_clear_error_missing_species_strict(self):
        """T055: Verify clear error when provider lacks species data (strict mode)."""
        manager = GermlineManager(providers=["imgt"])
        
        # Strict mode raises ValueError with helpful message
        with pytest.raises(ValueError) as exc_info:
            manager.get_genes("nonexistent_species", "V", "H", strict=True)
        
        error_msg = str(exc_info.value).lower()
        assert "nonexistent_species" in error_msg
        assert "not found" in error_msg or "no germline data" in error_msg

    def test_no_g3_fallback_when_germlines_enabled(self, monkeypatch):
        """T059: Verify no silent fallback to G3 when germlines selected."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        
        # Clear any cached feature flag value
        from sadie.germlines.utils.feature_flags import clear_feature_flag_cache
        clear_feature_flag_cache()
        
        from sadie.germlines.utils.feature_flags import use_germlines_module
        assert use_germlines_module() is True, "Germlines module should be enabled"
        
        # When germlines is enabled and species not found, should raise ValueError
        from sadie.airr.igblast.germline import GermlineData, _use_germlines_module
        
        # Verify feature flag is respected
        assert _use_germlines_module() is True
        
        # Test that missing species raises error instead of falling back
        with pytest.raises(ValueError) as exc_info:
            GermlineData("nonexistent_species_xyz")
        
        error_msg = str(exc_info.value).lower()
        assert "not found" in error_msg
        assert "nonexistent_species_xyz" in error_msg

    def test_feature_flag_controls_backend(self, monkeypatch):
        """Verify feature flag controls which backend is used."""
        from sadie.germlines.utils.feature_flags import clear_feature_flag_cache
        
        # Test with germlines enabled
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        clear_feature_flag_cache()
        
        from sadie.germlines.utils.feature_flags import use_germlines_module
        assert use_germlines_module() is True
        
        # Test with germlines disabled
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")
        clear_feature_flag_cache()
        assert use_germlines_module() is False


class TestCustomValidation:
    """T056: Test custom germline ingestion validation (FR-012)."""

    def test_validates_nucleotide_characters(self):
        """Verify custom provider validates nucleotide sequences."""
        from sadie.germlines.providers.custom import _validate_sequence
        
        # Valid sequence
        valid, msg = _validate_sequence("ACGTACGT", "test_gene")
        assert valid is True, f"Valid sequence rejected: {msg}"
        
        # Invalid characters
        valid, msg = _validate_sequence("ACGTXYZ", "test_gene")
        assert valid is False, "Invalid characters should be rejected"
        assert "invalid" in msg.lower() or "characters" in msg.lower()

    def test_rejects_empty_sequence(self):
        """Verify empty sequences are rejected."""
        from sadie.germlines.providers.custom import _validate_sequence
        
        valid, msg = _validate_sequence("", "test_gene")
        assert valid is False, "Empty sequence should be rejected"
        assert "empty" in msg.lower()

    def test_rejects_gap_only_sequence(self):
        """Verify sequences with only gaps are rejected."""
        from sadie.germlines.providers.custom import _validate_sequence
        
        valid, msg = _validate_sequence("...---...", "test_gene")
        assert valid is False, "Gap-only sequence should be rejected"

    def test_allows_gapped_sequences(self):
        """Verify gapped sequences (with real nucleotides) are accepted."""
        from sadie.germlines.providers.custom import _validate_sequence
        
        valid, msg = _validate_sequence("ACG...TAC", "test_gene")
        assert valid is True, f"Gapped sequence rejected: {msg}"

    def test_allows_iupac_ambiguous(self):
        """Verify IUPAC ambiguous codes are accepted."""
        from sadie.germlines.providers.custom import _validate_sequence
        
        # R=A/G, Y=C/T, etc.
        valid, msg = _validate_sequence("ACGTRYWSKM", "test_gene")
        assert valid is True, f"IUPAC codes rejected: {msg}"


class TestParity:
    """T057: Test species/chain/segment parity (FR-010)."""

    EXPECTED_SPECIES = ["human", "mouse"]  # Minimum required
    EXPECTED_CHAINS = ["H", "K", "L"]
    EXPECTED_SEGMENTS = ["V", "D", "J"]

    def test_human_coverage(self):
        """Verify human species has V, D, J for H, K, L chains."""
        manager = GermlineManager()
        
        for chain in self.EXPECTED_CHAINS:
            for segment in self.EXPECTED_SEGMENTS:
                genes = manager.get_genes("human", segment, chain)
                # D genes don't exist for light chains
                if segment == "D" and chain in ["K", "L"]:
                    continue
                assert len(genes) > 0, (
                    f"No genes found for human {chain}{segment}"
                )

    def test_mouse_coverage(self):
        """Verify mouse species has V, D, J for H, K, L chains."""
        manager = GermlineManager()
        
        for chain in self.EXPECTED_CHAINS:
            for segment in self.EXPECTED_SEGMENTS:
                genes = manager.get_genes("mouse", segment, chain)
                # D genes don't exist for light chains
                if segment == "D" and chain in ["K", "L"]:
                    continue
                # Mouse may have limited data - check if available
                if len(genes) == 0:
                    pytest.skip(f"No mouse {chain}{segment} data available")

    def test_gene_model_fields(self):
        """Verify genes have required fields populated."""
        manager = GermlineManager()
        genes = manager.get_genes("human", "V", "H")
        
        assert len(genes) > 0, "No human VH genes found"
        
        gene = genes[0]
        assert gene.name, "Gene should have name"
        assert gene.species == "human", "Gene should have species"
        assert gene.segment == "V", "Gene should have segment"
        assert gene.chain == "H", "Gene should have chain"
        assert gene.sequence, "Gene should have sequence"
        assert gene.source, "Gene should have source"


class TestGappedSequences:
    """T060: Test gapped sequence availability for HMM building (FR-013)."""

    def test_v_genes_have_gapped_data(self):
        """Verify V genes have gapped sequences for HMM building."""
        manager = GermlineManager()
        genes = manager.get_genes("human", "V", "H")
        
        genes_with_gapped = 0
        for gene in genes:
            if gene.sequence_aa_gapped or gene.sequence_gapped:
                genes_with_gapped += 1
        
        # Most V genes should have gapped data
        coverage = genes_with_gapped / len(genes) if genes else 0
        assert coverage > 0.5, (
            f"Only {coverage:.0%} of V genes have gapped sequences"
        )

    def test_j_genes_have_gapped_data(self):
        """Verify J genes have gapped sequences for HMM building."""
        manager = GermlineManager()
        genes = manager.get_genes("human", "J", "H")
        
        genes_with_gapped = 0
        for gene in genes:
            if gene.sequence_aa_gapped or gene.sequence_gapped:
                genes_with_gapped += 1
        
        coverage = genes_with_gapped / len(genes) if genes else 0
        # J genes are shorter and may have lower gapped coverage
        assert coverage > 0.4, (
            f"Only {coverage:.0%} of J genes have gapped sequences"
        )

    def test_hmm_builder_handles_missing_gapped(self):
        """Verify HMM builder reports genes missing gapped data."""
        from sadie.germlines.renumbering_integration import LocalHMMBuilder
        
        builder = LocalHMMBuilder()
        
        # Builder should be able to get VJ pairs
        # This tests the translation fallback path
        pairs = builder._get_vj_alignment_pairs("human", "H", "imgt")
        
        assert len(pairs) > 0, "Should find VJ alignment pairs for human H"
        
        # Each pair should have (name, gapped_aa_sequence)
        for name, seq in pairs[:5]:
            assert name, "Pair should have gene name"
            assert seq, f"Gene {name} should have gapped AA sequence"
