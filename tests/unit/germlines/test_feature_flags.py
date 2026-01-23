"""Tests for feature flag utilities."""

import os
import pytest
from sadie.germlines.utils.feature_flags import use_germlines_module, clear_feature_flag_cache


class TestUseGermlinesModule:
    @pytest.fixture(autouse=True)
    def clear_cache(self):
        """Clear the feature flag cache before each test."""
        clear_feature_flag_cache()
        yield
        clear_feature_flag_cache()
    """Test suite for use_germlines_module() feature flag function."""

    def test_default_behavior_returns_true(self, monkeypatch):
        """Test that default behavior uses germlines module (per FR-016)."""
        # Remove env var to test default
        monkeypatch.delenv("SADIE_USE_GERMLINES_MODULE", raising=False)

        assert use_germlines_module() is True

    def test_explicit_true_returns_true(self, monkeypatch):
        """Test that SADIE_USE_GERMLINES_MODULE=true uses germlines module (per FR-016a)."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        assert use_germlines_module() is True

    def test_explicit_false_returns_false(self, monkeypatch):
        """Test that SADIE_USE_GERMLINES_MODULE=false uses G3 API (per FR-016b)."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")

        assert use_germlines_module() is False

    def test_case_insensitive_true_variants(self, monkeypatch):
        """Test that various true values are recognized (case-insensitive)."""
        true_variants = ["true", "True", "TRUE", "1", "yes", "Yes", "YES", "on", "On", "ON"]

        for variant in true_variants:
            monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", variant)
            assert use_germlines_module() is True, f"Failed for variant: {variant}"

    def test_case_insensitive_false_variants(self, monkeypatch):
        """Test that false and unrecognized values return False."""
        false_variants = ["false", "False", "FALSE", "0", "no", "off", "invalid", ""]

        for variant in false_variants:
            monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", variant)
            assert use_germlines_module() is False, f"Failed for variant: {variant}"

    def test_deprecation_warning_when_false(self, monkeypatch, caplog):
        """Test that G3 mode logs deprecation warning (per FR-019a)."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")

        use_germlines_module()

        # Check that warning was logged
        assert any("G3 API is deprecated" in record.message for record in caplog.records)
        assert any("deprecated" in record.message.lower() for record in caplog.records)

    def test_no_warning_when_true(self, monkeypatch, caplog):
        """Test that germlines mode does not log deprecation warning."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        use_germlines_module()

        # Check that no deprecation warning was logged
        assert not any("deprecated" in record.message.lower() for record in caplog.records)

    def test_idempotent_calls(self, monkeypatch):
        """Test that multiple calls with same env var return consistent results."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

        result1 = use_germlines_module()
        result2 = use_germlines_module()
        result3 = use_germlines_module()

        assert result1 == result2 == result3 == True

    def test_env_var_changes_between_calls(self, monkeypatch):
        """Test that changing env var between calls updates behavior after cache clear."""
        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        assert use_germlines_module() is True

        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")
        clear_feature_flag_cache()
        assert use_germlines_module() is False

        monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
        clear_feature_flag_cache()
        assert use_germlines_module() is True
