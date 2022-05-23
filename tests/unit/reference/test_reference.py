"""
 Working directly with reference functions to create custom or trimmed databases. Also tests G3 intereaction
"""

# from pathlib import Path, PosixPath
from typing import List

import pandas as pd

# import yaml
from sadie.reference.reference import References

# from sadie.reference.util import write_out_fasta
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
import pytest
from sadie.reference.yaml import YamlRef
from tests.conftest import SadieFixture


def test_yaml(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    # load the default yaml file
    yaml_object = YamlRef()
    assert yaml_object.get_names() == {"clk", "dog", "human", "mouse", "rabbit", "se09", "rat", "macaque"}
    assert len(yaml_object.get_genes("human", "imgt", "human")) == 465
    v_genes: List[str] = yaml_object.get_gene_segment("human", "imgt", "human", "V")
    assert all([x[3] == "V" for x in v_genes])
    assert isinstance(yaml_object.get_yaml_as_dataframe(), pd.DataFrame)
    assert yaml_object.__repr__()
    assert set([i for i in yaml_object]) == {"clk", "dog", "human", "mouse", "rabbit", "se09", "rat", "macaque"}
    assert yaml_object["human"]
    assert len(yaml_object) == 4689

    with pytest.raises(ValueError):
        YamlRef(fixture_setup.get_duplicated_yaml())
    with pytest.raises(ValueError):
        YamlRef(fixture_setup.get_duplicated_diff_source_yaml())


def test_load_reference_from_yml(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    shortened_yaml = fixture_setup.get_shortened_yaml()
    references: References = References().from_yaml(shortened_yaml)
    outpath = tmp_path_factory.mktemp("test_load_reference_from_yml")
    references._make_igblast_ref_database(outpath)
    references._make_internal_annotaion_file(outpath)
    references._make_auxillary_file(outpath)


# def test_reference_class(tmp_path_factory: pytest.TempPathFactory) -> None:
#     """Test if we can JIT reference class."""
#     ref_class = Reference()

#     # Add one gene at a time
#     ref_class.add_gene({"species": "human", "gene": "IGHV1-69*01", "database": "imgt"})
#     ref_class.add_gene({"species": "human", "gene": "IGHD3-3*01", "database": "imgt"})
#     ref_class.add_gene({"species": "human", "gene": "IGHJ6*01", "database": "imgt"})
#     with pytest.raises(G3Error):
#         ref_class.add_gene({"species": "human", "gene": "IGHV111-69*01", "database": "imgt"})
#     # Add many genes at a time
#     genes: List[str] = []
#     genes.append("IGHV1-69*01")
#     genes.append("IGHD3-3*01")
#     genes.append("IGHJ6*01")
#     ref_class.add_genes("human", "imgt", genes)

#     references = References()
#     references.add_reference("human", ref_class)
#     dataframe = references.get_dataframe()
#     assert len(dataframe) == 3


# def test_references_class_from_yaml(tmp_path_factory: pytest.TempPathFactory) -> None:
#     references = References.from_yaml()
#     df: pd.DataFrame = references.get_dataframe()
#     print(df)


# def test_util_methods(tmpdir_factory):
#     # query = "https://g3.jordanrwillis.com/api/v1/genes?source=imgt&common=human&gene=IGHV1-69%2A01"
#     seq = "AAAAA"
#     file = Path(tmpdir_factory.mktemp("test_private_methods") + "/test.fasta")
#     write_out_fasta([SeqRecord(Seq(seq), id="test", name="test")], file)


# def test_creation_from_empty_reference(tmpdir_factory):
#     """Test when we create reference without pasing any data, it uses yaml"""
#     tmpdir = tmpdir_factory.mktemp("test_creation_from_empty_reference")
#     ref_class = Reference()
#     ref_class.make_airr_database(Path(tmpdir))


# def test_load_ref_from_df(fixture_setup, tmpdir_factory):
#     """Test if we can statically load a reference csv"""
#     ref_class = Reference.read_file(fixture_setup.get_reference_dataset_csv())
#     assert ref_class.data
#     outpath = tmpdir_factory.mktemp("test_load_ref_from_df")
#     outfile = pd.read_csv(fixture_setup.get_reference_dataset_csv(), index_col=0)
#     outfile.to_csv(outpath + "/test.csv")
#     outfile.to_feather(outpath + "/test.feather")  # pylint: disable=maybe-no-member
#     outfile.to_json(outpath + "/test.json", orient="records")  # pylint: disable=maybe-no-member
#     ref_class = Reference.read_file(outpath + "/test.json", type="json")
#     ref_class = Reference.read_file(outpath + "/test.csv")
#     ref_class = Reference.read_file(outpath + "/test.feather", type="feather")
#     with pytest.raises(ValueError):
#         ref_class = Reference.read_file(outpath + "/test.feather", type="oinga")


# def test_make_reference_class_from_yaml():
#     """Test reference class."""
#     ref_class = Reference.parse_yaml()
#     assert isinstance(ref_class, Reference)
#     ref_class_data = ref_class.get_dataframe()
#     assert isinstance(ref_class_data, pd.DataFrame)


# def test_check_models(tmpdir, fixture_setup):
#     """Coverage for all models including validations exceptions for bonehead entries"""
#     ref = Reference()
#     entry = {"species": "human", "sub_species": "human", "gene": "IGHV3-10*01", "database": "custom"}
#     GeneEntry(**entry)
#     entry = {"species": "human", "gene": "IGHV3-10*01", "database": "custom"}
#     GeneEntry(**entry)
#     with pytest.raises(ValidationError):
#         # Bad third postion
#         entry = {"species": "human", "gene": "IGHZ3-10*01", "database": "custom"}
#         GeneEntry(**entry)
#     with pytest.raises(ValidationError):
#         # Bad database
#         entry = {"species": "human", "gene": "IGHZ3-10*01", "database": "jordan_personal_stash"}
#         GeneEntry(**entry)

#     with pytest.raises(ValidationError):
#         # Bad third postion
#         entry = {"species": "human", "gene": ["IGHZ3-10*01", "IGHZ3-20*02"], "database": "custom"}
#         GeneEntries(**entry)
#     with pytest.raises(ValidationError):
#         # Bad database
#         entry = {"species": "human", "gene": ["IGHV3-10*01", "IGHV3-20*02"], "database": "jordan_personal_stash"}
#         GeneEntries(**entry)

#     with pytest.raises(ValueError):
#         ref._get_gene(GeneEntries)
#     with pytest.raises(ValueError):
#         ref._get_genes(GeneEntry)
#     entry = {"species": "mouse", "sub_species": "mouse", "gene": "IGHV2-6-3*01", "database": "imgt"}
#     ref._get_gene(GeneEntry(**entry))


# def test_yaml(tmpdir, fixture_setup):
#     """Test yaml module"""
#     yaml = YamlRef()
#     assert yaml.get_database_types() == {"imgt", "custom"}
#     imgt_keys = ["clk", "dog", "human", "macaque", "mouse", "rabbit", "rat", "se09"]

#     custom_keys = ["macaque"]
#     assert yaml.get_species_keys("imgt") == imgt_keys
#     assert yaml.get_species_keys("custom") == custom_keys

#     for key in imgt_keys:
#         sub_key = yaml.get_sub_species("imgt", key)
#         assert sub_key
#         for sub in sub_key:
#             genes = yaml.get_genes("imgt", key, sub)
#             yaml.get_gene_segment("imgt", key, sub, "V")
#             yaml.get_gene_segment("imgt", key, sub, "D")
#             yaml.get_gene_segment("imgt", key, sub, "J")
#             if not genes:
#                 raise ValueError(f"{key} {sub_key} {genes}")
#     for key in custom_keys:
#         sub_key = yaml.get_sub_species("custom", key)
#         assert sub_key
#         for sub in sub_key:
#             genes = yaml.get_genes("custom", key, sub)
#             yaml.get_gene_segment("custom", key, sub, "V")
#             yaml.get_gene_segment("custom", key, sub, "D")
#             yaml.get_gene_segment("custom", key, sub, "J")
#             if not genes:
#                 raise ValueError(f"{key} {sub_key} {genes}")
#     assert yaml.__repr__()


# def test_G3_errors():
#     """Test G3 errors"""
#     with pytest.raises(G3Error):
#         Reference(endpoint="https://mock.codes/202")


# def test_missing_makeblast_df(tmpdir, fixture_setup):
#     """Test the makeblast_df function with a missing blastdb"""
#     fasta = fixture_setup.get_catnap_light_nt()
#     bogus_file = fixture_setup.get_card()
#     from sadie.reference.blast import write_blast_db

#     with pytest.raises(ValueError):
#         write_blast_db(fasta, os.path.join(tmpdir, "missing.fasta"), "some_bogus_makeblastdb")
#     with pytest.raises(RuntimeError):
#         write_blast_db(bogus_file, os.path.join(tmpdir, "missing.fasta"))
