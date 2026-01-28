"""
 Working directly with reference functions to create custom or trimmed databases. Also tests G3 intereaction
"""
import glob
import io
import os
import sys
from typing import List

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from click.testing import CliRunner
from pydantic import ValidationError

from sadie import app
from sadie.reference.models import GeneEntries, GeneEntry
from sadie.reference.reference import G3Error, Reference, References
from sadie.reference.util import write_out_fasta
from sadie.reference.yaml import YamlRef
from tests.conftest import SadieFixture


def test_yaml(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    # load the default yaml file
    yaml_object = YamlRef()
    assert yaml_object.get_names() == {"clk", "dog", "human", "mouse", "rabbit", "se09", "rat", "macaque"}
    assert len(yaml_object.get_genes("human", "imgt", "human")) == 479
    v_genes: List[str] = yaml_object.get_gene_segment("human", "imgt", "human", "V")
    assert all([x[3] == "V" for x in v_genes])
    assert isinstance(yaml_object.get_yaml_as_dataframe(), pd.DataFrame)
    assert yaml_object.__repr__()
    assert set([i for i in yaml_object]) == {"clk", "dog", "human", "mouse", "rabbit", "se09", "rat", "macaque"}
    assert yaml_object["human"]
    assert len(yaml_object) == 4977

    with pytest.raises(ValueError):
        YamlRef(fixture_setup.get_duplicated_yaml())
    with pytest.raises(ValueError):
        YamlRef(fixture_setup.get_duplicated_diff_source_yaml())


def test_check_models(fixture_setup: SadieFixture) -> None:
    """Coverage for all models including validations exceptions for bonehead entries"""
    entry = {"species": "human", "gene": "IGHV3-10*01", "source": "custom"}
    GeneEntry(**entry)
    GeneEntries(species="human", genes=["IGHV1-68*01", "IGHV1-69*01"], source="custom")
    with pytest.raises(ValidationError):
        # bad third postion
        entry = {"species": "human", "gene": "IGHZ3-10*01", "source": "custom"}
        GeneEntry(**entry)
    with pytest.raises(ValidationError):
        # don't call it a database
        entry = {"species": "human", "gene": "IGHV1-69*01", "database": "custom"}
        GeneEntry(**entry)
    with pytest.raises(ValidationError):
        # Bad database
        entry = {"species": "human", "gene": "IGHZ3-10*01", "database": "jordan_personal_stash"}
        GeneEntry(**entry)


def test_source_validation_all_providers() -> None:
    """Test that all four germline sources are valid (v1.2 expansion)"""
    # All four sources should be valid
    for source in ["imgt", "ogrdb", "vdjbase", "custom"]:
        entry = GeneEntry(species="human", gene="IGHV1-69*01", source=source)
        assert entry.source == source

        entries = GeneEntries(species="human", genes=["IGHV1-69*01"], source=source)
        assert entries.source == source

    # Invalid source should raise ValidationError
    with pytest.raises(ValidationError) as exc_info:
        GeneEntry(species="human", gene="IGHV1-69*01", source="invalid_source")
    assert "not a valid source" in str(exc_info.value)
    assert "choices" in str(exc_info.value)  # verify typo fix

    with pytest.raises(ValidationError) as exc_info:
        GeneEntries(species="human", genes=["IGHV1-69*01"], source="g3")
    assert "not a valid source" in str(exc_info.value)


def test_reference_use_germlines() -> None:
    """Test Reference class with use_germlines=True (v1.2 integration)"""
    # Create reference with germlines enabled
    ref = Reference(use_germlines=True)
    assert ref.use_germlines is True

    # Add a gene using germlines module
    ref.add_gene({"species": "human", "gene": "IGHV1-69*01", "source": "imgt"})
    assert len(ref.data) == 1

    # Verify _id field is present (INT-03)
    df = ref.get_dataframe()
    assert "_id" in df.columns
    assert len(df["_id"].iloc[0]) == 24  # SHA-256 hash truncated to 24 chars

    # Verify same gene produces same _id (deterministic)
    ref2 = Reference(use_germlines=True)
    ref2.add_gene({"species": "human", "gene": "IGHV1-69*01", "source": "imgt"})
    df2 = ref2.get_dataframe()
    assert df["_id"].iloc[0] == df2["_id"].iloc[0]


def test_reference_use_germlines_missing_species() -> None:
    """Test germlines source/species validation for missing provider data."""
    ref = Reference(use_germlines=True)

    with pytest.raises(ValueError) as exc_info:
        ref.add_gene({"species": "mouse", "gene": "IGHV1-69*01", "source": "vdjbase"})

    message = str(exc_info.value)
    assert "vdjbase" in message
    assert "mouse" in message
    assert "Available species" in message


def test_references_from_yaml_use_germlines(fixture_setup: "SadieFixture") -> None:
    """Test References.from_yaml with use_germlines parameter (v1.2 INT-01)"""
    shortened_yaml = fixture_setup.get_shortened_yaml()

    # Load with use_germlines=True
    refs = References.from_yaml(shortened_yaml, use_germlines=True)

    # Verify references were created
    assert len(refs.references) > 0

    # Verify _id fields present in dataframe
    df = refs.get_dataframe()
    assert "_id" in df.columns
    assert df["_id"].notna().all()  # All rows have _id


def test_missing_imgt_positions_fail_fast(tmp_path_factory: pytest.TempPathFactory) -> None:
    from sadie.germlines import get_gene_by_name
    from sadie.germlines.g3_adapter import GermlineToG3Adapter

    gene = get_gene_by_name("IGHV1-69*01", "human")
    assert gene is not None

    adapter = GermlineToG3Adapter()
    g3_dict = adapter.to_g3_format(gene)

    for key in [
        "fwr1_start",
        "fwr1_end",
        "cdr1_start",
        "cdr1_end",
        "fwr2_start",
        "fwr2_end",
        "cdr2_start",
        "cdr2_end",
        "fwr3_start",
        "fwr3_end",
    ]:
        g3_dict["imgt"][key] = None

    ref = Reference(use_germlines=True)
    ref.data = [g3_dict]

    refs = References()
    refs.add_reference("test", ref)

    outpath = tmp_path_factory.mktemp("missing_imgt_positions")
    with pytest.raises(ValueError) as exc_info:
        refs.make_airr_database(outpath)

    message = str(exc_info.value)
    assert "Missing IMGT V-region positions" in message
    assert "IGHV1-69*01" in message


def test_util_methods(tmp_path_factory: pytest.TempPathFactory) -> None:
    seq = "AAAAA"
    file = tmp_path_factory.mktemp("test_private_methods").joinpath("test.fasta")
    write_out_fasta([SeqRecord(Seq(seq), id="test", name="test")], file)


def test_check_default_reference_df(fixture_setup: SadieFixture) -> None:
    """Test if we have default reference"""
    ref_api = References()
    df = ref_api.get_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4977


def test_load_reference_from_yml(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    shortened_yaml = fixture_setup.get_shortened_yaml()
    references: References = References().from_yaml(shortened_yaml)
    outpath = tmp_path_factory.mktemp("test_load_reference_from_yml")
    output_db = references.make_airr_database(outpath)
    # HMM generation adds hmms/ and stockholms/ directories
    assert sorted([i.name for i in output_db.glob("*")]) == sorted(
        [".references_dataframe.csv.gz", "aux_db", "Ig", "hmms", "stockholms"]
    )

    # test we can get a dataframe

    df = references.get_dataframe()
    assert isinstance(df, pd.DataFrame)

    # make a new reference from the dataframe
    new_ref = References().from_dataframe(df)
    pd._testing.assert_frame_equal(new_ref.get_dataframe(), df)


def test_reference_class(tmp_path_factory: pytest.TempPathFactory) -> None:
    """Test if we can JIT reference class."""
    ref_class = Reference()

    # Add one gene at a time
    ref_class.add_gene({"species": "human", "gene": "IGHV1-69*01", "source": "imgt"})
    ref_class.add_gene({"species": "human", "gene": "IGHD3-3*01", "source": "imgt"})
    ref_class.add_gene({"species": "human", "gene": "IGHJ6*01", "source": "imgt"})
    with pytest.raises(G3Error):
        # G3 does not have Gene
        ref_class.add_gene({"species": "human", "gene": "IGHV111-69*01", "source": "imgt"})
    with pytest.raises(ValidationError):
        # don't call it a database
        ref_class.add_gene({"species": "human", "gene": "IGHV111-69*01", "database": "imgt"})
    # Add many genes at a time
    genes: List[str] = []
    genes.append("IGHV1-69*01")
    genes.append("IGHD3-3*01")
    genes.append("IGHJ6*01")
    ref_class.add_genes("human", "imgt", genes)

    references = References()

    references.add_reference("human", ref_class)
    dataframe = references.get_dataframe()
    assert len(dataframe) == 3


def test_make_from_empty(tmp_path_factory: pytest.TempPathFactory) -> None:
    """Test when we create reference without pasing any data, it uses yaml"""
    tmpdir = tmp_path_factory.mktemp("test_creation_from_empty_reference")
    ref_class = References()
    output = ref_class.make_airr_database(tmpdir)
    assert sorted([i.name for i in output.glob("*")]) == sorted([".references_dataframe.csv.gz", "aux_db", "Ig"])


def test_G3_errors() -> None:
    """Test G3 errors"""
    with pytest.raises(G3Error):
        Reference(endpoint="https://mock.codes/202")


def test_missing_makeblast_df(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    """Test the makeblast_df function with a missing blastdb"""
    fasta = fixture_setup.get_catnap_light_nt()
    bogus_file = fixture_setup.get_card()
    from sadie.reference.util import write_blast_db

    tmpdir = tmp_path_factory.mktemp("test_missing_makeblast_df")
    with pytest.raises(ValueError):
        write_blast_db(fasta, tmpdir.joinpath("missing.fasta"), "some_bogus_makeblastdb")
    with pytest.raises(RuntimeError):
        write_blast_db(bogus_file, tmpdir.joinpath("missing.fasta"))


def test_make_hmm_files_creates_directories(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    """Test that _make_hmm_files creates stockholms/ and hmms/ directories."""
    shortened_yaml = fixture_setup.get_shortened_yaml()
    references: References = References().from_yaml(shortened_yaml, use_germlines=True)
    outpath = tmp_path_factory.mktemp("test_hmm_directories")

    # Build the database (which calls _make_hmm_files internally)
    output_db = references.make_airr_database(outpath)

    # Verify directories exist
    stockholms_dir = output_db / "stockholms"
    hmms_dir = output_db / "hmms"
    assert stockholms_dir.exists(), "stockholms directory should be created"
    assert hmms_dir.exists(), "hmms directory should be created"


def test_make_hmm_files_generates_hmm(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    """Test that _make_hmm_files generates .hmm files for available chains."""
    import pyhmmer

    shortened_yaml = fixture_setup.get_shortened_yaml()
    references: References = References().from_yaml(shortened_yaml, use_germlines=True)
    outpath = tmp_path_factory.mktemp("test_hmm_generation")

    # Build the database
    output_db = references.make_airr_database(outpath)

    hmms_dir = output_db / "hmms"
    stockholms_dir = output_db / "stockholms"

    # Check that at least one HMM was created
    hmm_files = list(hmms_dir.glob("*.hmm"))
    assert len(hmm_files) > 0, "At least one HMM file should be generated"

    # Check that corresponding Stockholm files were created
    sto_files = list(stockholms_dir.glob("*.sto"))
    assert len(sto_files) > 0, "At least one Stockholm file should be generated"

    # Verify HMM files are valid by loading them with pyhmmer
    for hmm_file in hmm_files:
        with pyhmmer.plan7.HMMFile(hmm_file) as hmm_reader:
            hmm = next(hmm_reader)
            assert hmm is not None, f"HMM file {hmm_file} should be loadable"
            assert hmm.M > 0, f"HMM should have positions (M > 0)"


def test_make_hmm_files_handles_missing_gapped_sequences(tmp_path_factory: pytest.TempPathFactory) -> None:
    """Test that _make_hmm_files handles missing gapped sequences gracefully."""
    from sadie.germlines import get_gene_by_name
    from sadie.germlines.g3_adapter import GermlineToG3Adapter

    # Create a reference with genes that have no gapped AA sequences
    gene = get_gene_by_name("IGHV1-69*01", "human")
    assert gene is not None

    adapter = GermlineToG3Adapter()
    g3_dict = adapter.to_g3_format(gene)

    # Remove gapped AA sequence to simulate missing data
    g3_dict["imgt"]["sequence_gapped_aa"] = None

    ref = Reference(use_germlines=True)
    ref.data = [g3_dict]

    refs = References()
    refs.add_reference("test_missing", ref)

    outpath = tmp_path_factory.mktemp("test_hmm_missing_gapped")

    # Build should complete without error, even with missing gapped sequences
    # (make_airr_database may fail due to missing D gene, so we call _make_hmm_files directly)
    refs._make_hmm_files(outpath)

    # Directories should still be created
    stockholms_dir = outpath / "stockholms"
    hmms_dir = outpath / "hmms"
    assert stockholms_dir.exists(), "stockholms directory should be created"
    assert hmms_dir.exists(), "hmms directory should be created"


def test_cli(tmp_path_factory: pytest.TempPathFactory):
    """Confirm the CLI works as expected This runs the entire generation pipeline that ships with SADIE and checks that the file structure is exactly the same"""
    # Create runner
    runner = CliRunner()

    # these are the expected file structures
    # make a hierarchy of directories
    tmpdir = tmp_path_factory.mktemp("igblast_dir")

    # run the entire pipeline via CLICK cli
    # Pass string path to avoid any Path object serialization issues
    # Capture output to avoid I/O closed file errors in CI
    stdout_capture = io.StringIO()
    stderr_capture = io.StringIO()

    old_stdout = sys.stdout
    old_stderr = sys.stderr

    try:
        # Only redirect if we're in CI environment
        if os.environ.get("CI") == "true" or os.environ.get("GITHUB_ACTIONS") == "true":
            sys.stdout = stdout_capture
            sys.stderr = stderr_capture

        result = runner.invoke(app.make_igblast_reference, ["--outpath", str(tmpdir)], catch_exceptions=True)
    finally:
        # Restore original stdout/stderr
        sys.stdout = old_stdout
        sys.stderr = old_stderr

    if result.exit_code != 0:
        # Check if there's an exception before accessing output
        if result.exception:
            import traceback

            tb = "".join(
                traceback.format_exception(type(result.exception), result.exception, result.exception.__traceback__)
            )
            assert False, f"Command failed with exception:\n{tb}"
        else:
            # Only try to access output if no exception
            # Safely handle potential I/O errors
            try:
                output = result.output if hasattr(result, "output") else "No output available"
            except (ValueError, IOError):
                output = "Output unavailable due to I/O error"
            assert result.exit_code == 0, f"Command failed with output:\n{output}"

    # was the file actually output?
    assert os.path.exists(tmpdir)

    # assert we made an imgt and custom directory, but still don't know if anything is in it
    directories_created = glob.glob(str(tmpdir) + "/*")
    assert sorted(directories_created) == sorted([f"{tmpdir}/aux_db", f"{tmpdir}/Ig"])
