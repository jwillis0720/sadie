"""
 Working directly with reference functions to create custom or trimmed databases. Also tests G3 intereaction
"""
from typing import List

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydantic import ValidationError

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
    assert len(yaml_object) == 4976

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


def test_util_methods(tmp_path_factory: pytest.TempPathFactory) -> None:
    seq = "AAAAA"
    file = tmp_path_factory.mktemp("test_private_methods").joinpath("test.fasta")
    write_out_fasta([SeqRecord(Seq(seq), id="test", name="test")], file)


def test_check_default_reference_df(fixture_setup: SadieFixture) -> None:
    """Test if we have default reference"""
    ref_api = References()
    df = ref_api.get_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4976


def test_load_reference_from_yml(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    shortened_yaml = fixture_setup.get_shortened_yaml()
    references: References = References().from_yaml(shortened_yaml)
    outpath = tmp_path_factory.mktemp("test_load_reference_from_yml")
    output_db = references.make_airr_database(outpath)
    assert sorted([i.name for i in output_db.glob("*")]) == sorted([".references_dataframe.csv.gz", "aux_db", "Ig"])

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
