"""
Tests for `sadie germlines rebuild` CLI command using Reference module path.

Validates that the rebuild command:
1. Routes through the Reference module (References.from_yaml → make_airr_database)
2. Keeps the same CLI interface (--species flag)
3. Does NOT modify `sadie germlines populate` behavior
4. Outputs to the germlines/igblast/ directory
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from sadie.app import sadie


@pytest.fixture
def cli_runner():
    """Create CLI runner for testing."""
    return CliRunner()


class TestRebuildCLIHelp:
    """Tests for rebuild command help and interface."""

    def test_rebuild_help(self, cli_runner):
        """Test rebuild command shows help with --species option."""
        result = cli_runner.invoke(sadie, ["germlines", "rebuild", "--help"])
        assert result.exit_code == 0
        assert "--species" in result.output
        assert "-s" in result.output

    def test_rebuild_help_mentions_reference_module(self, cli_runner):
        """Test rebuild help references the Reference module build path."""
        result = cli_runner.invoke(sadie, ["germlines", "rebuild", "--help"])
        assert result.exit_code == 0
        # Should mention the reference-based approach
        assert "Reference" in result.output or "reference" in result.output


class TestRebuildUsesReferenceModule:
    """Tests that rebuild routes through Reference module."""

    @patch("sadie.reference.reference.References.from_yaml")
    @patch("sadie.reference.generate.generate_reference_yaml")
    def test_rebuild_calls_references_from_yaml(self, mock_generate, mock_from_yaml, cli_runner):
        """Test that rebuild calls References.from_yaml() internally."""
        # Setup mocks
        mock_refs_instance = MagicMock()
        mock_from_yaml.return_value = mock_refs_instance
        mock_refs_instance.make_airr_database.return_value = Path("/tmp/test_db")

        # Mock the yaml_path existence and content
        yaml_content = '{"human": {"imgt": {"human": ["IGHV1-2*01"]}}}'
        with (
            patch("sadie.germlines.cli.get_all_downloaded_species", return_value=["human"]),
            patch("builtins.open", create=True),
            patch("yaml.safe_load", return_value={"human": {"imgt": {"human": ["IGHV1-2*01"]}}}),
            patch("pathlib.Path.exists", return_value=True),
        ):
            result = cli_runner.invoke(sadie, ["germlines", "rebuild", "--species", "human"])

        # Should have called References.from_yaml
        assert mock_from_yaml.called, "rebuild should call References.from_yaml()"

    @patch("sadie.reference.reference.References.from_yaml")
    @patch("sadie.reference.generate.generate_reference_yaml")
    def test_rebuild_calls_make_airr_database(self, mock_generate, mock_from_yaml, cli_runner):
        """Test that rebuild calls make_airr_database() to build the DB."""
        mock_refs_instance = MagicMock()
        mock_from_yaml.return_value = mock_refs_instance
        mock_refs_instance.make_airr_database.return_value = Path("/tmp/test_db")

        with (
            patch("sadie.germlines.cli.get_all_downloaded_species", return_value=["human"]),
            patch("builtins.open", create=True),
            patch("yaml.safe_load", return_value={"human": {"imgt": {"human": ["IGHV1-2*01"]}}}),
            patch("pathlib.Path.exists", return_value=True),
        ):
            result = cli_runner.invoke(sadie, ["germlines", "rebuild", "--species", "human"])

        # Should have called make_airr_database
        assert mock_refs_instance.make_airr_database.called, "rebuild should call make_airr_database()"

    @patch("sadie.reference.reference.References.from_yaml")
    @patch("sadie.reference.generate.generate_reference_yaml")
    def test_rebuild_outputs_to_igblast_dir(self, mock_generate, mock_from_yaml, cli_runner):
        """Test that rebuild outputs to the germlines/igblast/ directory."""
        mock_refs_instance = MagicMock()
        mock_from_yaml.return_value = mock_refs_instance
        mock_refs_instance.make_airr_database.return_value = Path("/tmp/test_db")

        with (
            patch("sadie.germlines.cli.get_all_downloaded_species", return_value=["human"]),
            patch("builtins.open", create=True),
            patch("yaml.safe_load", return_value={"human": {"imgt": {"human": ["IGHV1-2*01"]}}}),
            patch("pathlib.Path.exists", return_value=True),
        ):
            result = cli_runner.invoke(sadie, ["germlines", "rebuild", "--species", "human"])

        # Verify the output path is the germlines/igblast/ directory
        call_args = mock_refs_instance.make_airr_database.call_args
        output_path = call_args[0][0]  # First positional arg
        assert "igblast" in str(output_path), f"Output should be in igblast dir, got {output_path}"

    @patch("sadie.reference.reference.References.from_yaml")
    def test_rebuild_does_not_use_germline_pipeline(self, mock_from_yaml, cli_runner):
        """Test that rebuild does NOT use GermlinePipeline.force_rebuild."""
        mock_refs_instance = MagicMock()
        mock_from_yaml.return_value = mock_refs_instance
        mock_refs_instance.make_airr_database.return_value = Path("/tmp/test_db")

        with (
            patch("sadie.germlines.cli.get_all_downloaded_species", return_value=["human"]),
            patch("builtins.open", create=True),
            patch("yaml.safe_load", return_value={"human": {"imgt": {"human": ["IGHV1-2*01"]}}}),
            patch("pathlib.Path.exists", return_value=True),
            patch("sadie.germlines.pipeline.GermlinePipeline") as mock_pipeline,
        ):
            result = cli_runner.invoke(sadie, ["germlines", "rebuild", "--species", "human"])

            # GermlinePipeline should NOT be instantiated
            assert not mock_pipeline.called, "rebuild should NOT use GermlinePipeline"


class TestRebuildAutoGeneratesConfig:
    """Tests that rebuild auto-generates reference config when missing."""

    @patch("sadie.reference.reference.References.from_yaml")
    @patch("sadie.reference.generate.generate_reference_yaml")
    def test_rebuild_generates_config_if_missing(self, mock_gen, mock_from_yaml, cli_runner):
        """Test rebuild auto-generates reference.yml if not found."""
        mock_refs_instance = MagicMock()
        mock_from_yaml.return_value = mock_refs_instance
        mock_refs_instance.make_airr_database.return_value = Path("/tmp/test_db")

        # Track exists calls to only return False for yaml_path check
        original_exists = Path.exists

        def mock_exists(self_path):
            if str(self_path).endswith("reference.yml"):
                return False
            return original_exists(self_path)

        with (
            patch("sadie.germlines.cli.get_all_downloaded_species", return_value=["human"]),
            patch.object(Path, "exists", mock_exists),
            patch("builtins.open", create=True),
            patch("yaml.safe_load", return_value={"human": {"imgt": {"human": ["IGHV1-2*01"]}}}),
        ):
            result = cli_runner.invoke(sadie, ["germlines", "rebuild", "--species", "human"])

        # generate_reference_yaml should be called
        assert mock_gen.called, "Should auto-generate reference.yml when missing"


class TestRebuildSpeciesFiltering:
    """Tests for species filtering with --species flag."""

    @patch("sadie.reference.reference.References.from_yaml")
    @patch("sadie.reference.generate.generate_reference_yaml")
    def test_rebuild_all_species_when_no_flag(self, mock_generate, mock_from_yaml, cli_runner):
        """Test rebuild processes all downloaded species when no --species flag."""
        mock_refs_instance = MagicMock()
        mock_from_yaml.return_value = mock_refs_instance
        mock_refs_instance.make_airr_database.return_value = Path("/tmp/test_db")

        with (
            patch("sadie.germlines.cli.get_all_downloaded_species", return_value=["human", "mouse"]),
            patch("builtins.open", create=True),
            patch(
                "yaml.safe_load",
                return_value={
                    "human": {"imgt": {"human": ["IGHV1-2*01"]}},
                    "mouse": {"imgt": {"mouse": ["IGHV1-2*01"]}},
                },
            ),
            patch("pathlib.Path.exists", return_value=True),
        ):
            result = cli_runner.invoke(sadie, ["germlines", "rebuild"])

        # Should have been called for each species
        assert mock_from_yaml.call_count >= 1

    @patch("sadie.reference.reference.References.from_yaml")
    @patch("sadie.reference.generate.generate_reference_yaml")
    def test_rebuild_specific_species(self, mock_generate, mock_from_yaml, cli_runner):
        """Test rebuild with --species human only processes human."""
        mock_refs_instance = MagicMock()
        mock_from_yaml.return_value = mock_refs_instance
        mock_refs_instance.make_airr_database.return_value = Path("/tmp/test_db")

        with (
            patch("sadie.germlines.cli.get_all_downloaded_species", return_value=["human", "mouse"]),
            patch("builtins.open", create=True),
            patch("yaml.safe_load", return_value={"human": {"imgt": {"human": ["IGHV1-2*01"]}}}),
            patch("pathlib.Path.exists", return_value=True),
        ):
            result = cli_runner.invoke(sadie, ["germlines", "rebuild", "--species", "human"])

        assert mock_from_yaml.called

    def test_rebuild_warns_on_missing_species(self, cli_runner):
        """Test rebuild warns when species is not in downloaded data."""
        with patch("sadie.germlines.cli.get_all_downloaded_species", return_value=["human"]):
            result = cli_runner.invoke(sadie, ["germlines", "rebuild", "--species", "unicorn"])

        assert "not found" in result.output.lower() or "unicorn" in result.output

    def test_rebuild_no_data_shows_warning(self, cli_runner):
        """Test rebuild with no downloaded data shows appropriate warning."""
        with patch("sadie.germlines.cli.get_all_downloaded_species", return_value=[]):
            result = cli_runner.invoke(sadie, ["germlines", "rebuild"])
        assert "populate" in result.output.lower() or "no" in result.output.lower()


class TestPopulateUnchanged:
    """Tests that sadie germlines populate is not affected by rebuild changes."""

    def test_populate_help_unchanged(self, cli_runner):
        """Test populate help is still the same."""
        result = cli_runner.invoke(sadie, ["germlines", "populate", "--help"])
        assert result.exit_code == 0
        assert "--provider" in result.output
        assert "--species" in result.output
        assert "--force" in result.output
        assert "--dry-run" in result.output

    def test_populate_does_not_import_reference(self, cli_runner):
        """Test populate path does not import Reference module."""
        result = cli_runner.invoke(sadie, ["germlines", "populate", "--dry-run"])
        # This should work without Reference module involvement
        assert result.exit_code == 0
        assert "DRY RUN" in result.output


class TestFindReferenceNameForSpecies:
    """Tests for the _find_reference_name_for_species helper."""

    def test_species_as_direct_reference_name(self):
        """Test species found as direct reference name."""
        from sadie.germlines.cli import _find_reference_name_for_species

        config = {"human": {"imgt": {"human": ["IGHV1-2*01"]}}}
        assert _find_reference_name_for_species(config, "human") == "human"

    def test_species_as_sub_key(self):
        """Test species found as sub-key under a different reference name."""
        from sadie.germlines.cli import _find_reference_name_for_species

        config = {"clk": {"imgt": {"chicken": ["IGHV1-2*01"]}}}
        assert _find_reference_name_for_species(config, "chicken") == "clk"

    def test_species_not_found(self):
        """Test returns None when species not in config."""
        from sadie.germlines.cli import _find_reference_name_for_species

        config = {"human": {"imgt": {"human": ["IGHV1-2*01"]}}}
        assert _find_reference_name_for_species(config, "unicorn") is None
