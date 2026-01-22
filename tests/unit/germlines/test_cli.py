"""
Tests for germlines CLI commands.

Tests the `sadie germlines populate` and `sadie germlines status` CLI commands.
"""

import json
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
from click.testing import CliRunner

from sadie.app import sadie


@pytest.fixture
def cli_runner():
    """Create CLI runner for testing."""
    return CliRunner()


class TestGermlinesPopulateHelp:
    """Tests for CLI help and basic functionality."""

    def test_germlines_help(self, cli_runner):
        """Test germlines group help."""
        result = cli_runner.invoke(sadie, ["germlines", "--help"])
        assert result.exit_code == 0
        assert "Germline database management commands" in result.output

    def test_populate_help(self, cli_runner):
        """Test populate command help."""
        result = cli_runner.invoke(sadie, ["germlines", "populate", "--help"])
        assert result.exit_code == 0
        assert "--provider" in result.output
        assert "--species" in result.output
        assert "--force" in result.output
        assert "--dry-run" in result.output

    def test_status_help(self, cli_runner):
        """Test status command help."""
        result = cli_runner.invoke(sadie, ["germlines", "status", "--help"])
        assert result.exit_code == 0
        assert "status" in result.output.lower()


class TestGermlinesPopulateDryRun:
    """Tests for dry-run functionality."""

    def test_populate_dry_run_all(self, cli_runner):
        """Test dry-run mode shows what would be downloaded."""
        result = cli_runner.invoke(sadie, ["germlines", "populate", "--dry-run"])
        assert result.exit_code == 0
        assert "DRY RUN" in result.output
        assert "Would download" in result.output

    def test_populate_dry_run_imgt(self, cli_runner):
        """Test dry-run mode for IMGT provider."""
        result = cli_runner.invoke(
            sadie, ["germlines", "populate", "--provider", "imgt", "--dry-run"]
        )
        assert result.exit_code == 0
        assert "imgt" in result.output.lower()
        assert "Would download" in result.output

    def test_populate_dry_run_single_species(self, cli_runner):
        """Test dry-run mode for single species."""
        result = cli_runner.invoke(
            sadie,
            ["germlines", "populate", "--provider", "imgt", "--species", "human", "--dry-run"],
        )
        assert result.exit_code == 0
        assert "human" in result.output


class TestGermlinesStatus:
    """Tests for status command."""

    def test_status_command(self, cli_runner):
        """Test status command shows table."""
        result = cli_runner.invoke(sadie, ["germlines", "status"])
        assert result.exit_code == 0
        assert "imgt" in result.output.lower()
        assert "ogrdb" in result.output.lower()
        assert "vdjbase" in result.output.lower()


class TestVersionTracking:
    """Tests for version tracking functionality."""

    def test_get_local_version_missing(self, tmp_path):
        """Test get_local_version returns None when no version file."""
        from sadie.germlines.cli import get_local_version
        
        with patch("sadie.germlines.cli.get_provider_data_dir", return_value=tmp_path):
            result = get_local_version("imgt")
            assert result is None

    def test_get_local_version_exists(self, tmp_path):
        """Test get_local_version reads version file."""
        from sadie.germlines.cli import get_local_version, save_version
        
        with patch("sadie.germlines.cli.get_provider_data_dir", return_value=tmp_path):
            save_version("imgt", "release-202601", 29)
            result = get_local_version("imgt")
            
            assert result is not None
            assert result["version"] == "release-202601"
            assert result["species_count"] == 29

    def test_is_up_to_date(self, tmp_path):
        """Test is_up_to_date checks version."""
        from sadie.germlines.cli import is_up_to_date, save_version, get_current_version_string
        
        with patch("sadie.germlines.cli.get_provider_data_dir", return_value=tmp_path):
            assert not is_up_to_date("imgt")
            
            current = get_current_version_string()
            save_version("imgt", current, 29)
            
            assert is_up_to_date("imgt")


class TestCheckpointHandling:
    """Tests for checkpoint/resume functionality."""

    def test_checkpoint_save_load(self, tmp_path):
        """Test checkpoint save and load."""
        from sadie.germlines.cli import save_checkpoint, load_checkpoint, clear_checkpoint
        
        with patch("sadie.germlines.cli.get_provider_data_dir", return_value=tmp_path):
            completed = {"human", "mouse", "rat"}
            save_checkpoint("imgt", completed, 10)
            
            loaded = load_checkpoint("imgt")
            assert loaded == completed
            
            clear_checkpoint("imgt")
            assert load_checkpoint("imgt") == set()

    def test_checkpoint_empty_on_missing(self, tmp_path):
        """Test checkpoint returns empty set when no checkpoint exists."""
        from sadie.germlines.cli import load_checkpoint
        
        with patch("sadie.germlines.cli.get_provider_data_dir", return_value=tmp_path):
            result = load_checkpoint("imgt")
            assert result == set()


class TestProviderHelpers:
    """Tests for provider helper functions."""

    def test_get_provider_imgt(self):
        """Test get_provider returns IMGT provider."""
        from sadie.germlines.cli import get_provider
        
        provider = get_provider("imgt")
        assert provider.name == "imgt"

    def test_get_provider_ogrdb(self):
        """Test get_provider returns OGRDB provider."""
        from sadie.germlines.cli import get_provider
        
        provider = get_provider("ogrdb")
        assert provider.name == "ogrdb"

    def test_get_provider_vdjbase(self):
        """Test get_provider returns VDJbase provider."""
        from sadie.germlines.cli import get_provider
        
        provider = get_provider("vdjbase")
        assert provider.name == "vdjbase"

    def test_get_provider_invalid(self):
        """Test get_provider raises for invalid provider."""
        from sadie.germlines.cli import get_provider
        
        with pytest.raises(ValueError, match="Unknown provider"):
            get_provider("invalid")

    def test_get_all_provider_species_imgt(self):
        """Test get_all_provider_species returns IMGT species."""
        from sadie.germlines.cli import get_all_provider_species
        
        species = get_all_provider_species("imgt")
        assert "human" in species
        assert "mouse" in species
        assert len(species) > 10

    def test_get_all_provider_species_ogrdb(self):
        """Test get_all_provider_species returns OGRDB species."""
        from sadie.germlines.cli import get_all_provider_species
        
        species = get_all_provider_species("ogrdb")
        assert "human" in species
        assert "mouse" in species


class TestDataValidation:
    """Tests for data validation."""

    def test_validate_provider_data_empty(self, tmp_path):
        """Test validation fails for empty provider."""
        from sadie.germlines.cli import validate_provider_data
        
        with patch("sadie.germlines.cli.get_provider") as mock_get:
            mock_provider = MagicMock()
            mock_provider.get_metadata.return_value = MagicMock(species_available=[])
            mock_get.return_value = mock_provider
            
            result = validate_provider_data("imgt")
            assert result is False
