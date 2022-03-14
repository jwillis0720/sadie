import glob
import os
from pydantic import ValidationError
import pandas as pd
import pytest
from click.testing import CliRunner
from sadie import app
from sadie.reference import Reference
from sadie.reference.reference import G3Error, get_database, get_loaded_database, get_species_from_database
from sadie.reference.models import GeneEntry, GeneEntries
from sadie.reference.yaml import YamlRef

known_aux_exceptions = {
    ("mouse", "imgt", "IGLJ4*01"): "igblast has wrong number of c-term remaining",
    ("rabbit", "imgt", "IGHJ3*02"): "igblast has wrong reading frame",
}


def test_reference_class():
    """Test if we can JIT reference class."""
    ref_class = Reference()
    ref_class.add_gene({"species": "human", "gene": "IGHV1-69*01", "database": "imgt"})
    ref_class.add_gene({"species": "human", "gene": "IGHD3-3*01", "database": "imgt"})
    with pytest.raises(G3Error):
        ref_class.add_gene({"species": "human", "gene": "IGHV111-69*01", "database": "imgt"})
    assert len(ref_class.get_dataframe()) == 2


def test_load_ref_from_df(fixture_setup):
    """Test if we can statically load a reference csv"""
    ref_class = Reference.read_file(fixture_setup.get_reference_dataset_csv())
    assert ref_class.data


def test_make_reference_class_from_yaml():
    """Test reference class."""
    ref_class = Reference.parse_yaml()
    assert isinstance(ref_class, Reference)
    ref_class_data = ref_class.get_dataframe()
    assert isinstance(ref_class_data, pd.DataFrame)


def _test_auxilary_file_structure(tmpdir, fixture_setup):
    # Send aux and internal to compare against IGBLAST internal data and aux data
    my_aux = []
    generated_aux_path_files = glob.glob(f"{tmpdir}/*/**/*.aux", recursive=True)
    for file in generated_aux_path_files:
        df = pd.read_csv(
            file,
            skip_blank_lines=True,
            delimiter="\t",
            header=None,
            names=["gene", "reading_frame", "segment", "cdr3_end", "left_over"],
        )
        df.insert(0, "common", os.path.basename(file).split("_")[0])
        df.insert(1, "db_type", file.split("/")[-3])
        my_aux.append(df)
    my_aux = pd.concat(my_aux).reset_index(drop=True).set_index(["common", "db_type", "gene"])

    igblast_aux = []
    igblast_aux_files = fixture_setup.get_aux_files()
    for file in igblast_aux_files:
        df = pd.read_csv(
            file,
            skip_blank_lines=True,
            delim_whitespace=True,
            skiprows=2,
            header=None,
            names=["gene", "reading_frame", "segment", "cdr3_end", "left_over"],
        )
        df.insert(0, "common", os.path.basename(file).split("_")[0])
        df.insert(1, "db_type", "imgt")
        igblast_aux.append(df)
    igblast_aux = pd.concat(igblast_aux).reset_index(drop=True).set_index(["common", "db_type", "gene"])
    common_index = my_aux.index.intersection(igblast_aux.index)
    # sadie_missing_aux = igblast_aux.index.difference(my_aux.index)
    # igblast_missing_sadie = my_aux.index.difference(igblast_aux.index)

    my_aux_common_index = my_aux.loc[common_index]
    igblast_common_index = igblast_aux.loc[common_index]
    for index in my_aux_common_index.index:
        try:
            pd._testing.assert_series_equal(igblast_common_index.loc[index], my_aux_common_index.loc[index])
        except AssertionError:
            if index in known_aux_exceptions.keys():
                print(
                    index,
                    "is known exception exception",
                    known_aux_exceptions[index],
                    "skipping",
                )
                continue
            else:
                # raise again since pandas gives way better info
                pd._testing.assert_series_equal(
                    igblast_common_index.loc[index],
                    my_aux_common_index.loc[index],
                    obj=index,
                )
    return True


def _test_internal_data_file_structure(tmpdir, fixture_setup):
    # what we have made
    internal_path = glob.glob(f"{tmpdir}/imgt/**/*.imgt", recursive=True)
    reference_internal_path = fixture_setup.get_internal_files()
    my_internal_path_df = []
    for file in internal_path:
        df = pd.read_csv(
            file,
            skip_blank_lines=True,
            delimiter="\t",
            header=None,
            names=[
                "gene",
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
                "segment",
                "weird",
            ],
        )
        df.insert(0, "common", os.path.basename(file).split(".ndm")[0])
        df.insert(1, "db_type", "imgt")
        my_internal_path_df.append(df)

    my_internal_path_df = (
        pd.concat(my_internal_path_df).reset_index(drop=True).groupby(["common", "db_type", "gene"]).head(1)
    )
    my_internal_path_df = my_internal_path_df.set_index(["common", "db_type", "gene"])
    ref_internal_path_df = []
    for file in reference_internal_path:
        df = pd.read_csv(
            file,
            skip_blank_lines=True,
            delimiter="\t",
            skiprows=2,
            header=None,
            names=[
                "gene",
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
                "segment",
                "weird",
            ],
        )
        common_name = os.path.basename(file).split(".ndm")[0]
        if common_name == "rhesus_monkey":
            print("here")
            common_name = "macaque"

        df.insert(0, "common", common_name)
        df.insert(1, "db_type", "imgt")
        ref_internal_path_df.append(df)

    ref_internal_path_df = pd.concat(ref_internal_path_df).reset_index(drop=True)
    ref_internal_path_df = ref_internal_path_df.groupby(["common", "db_type", "gene"]).head(1)
    ref_internal_path_df = ref_internal_path_df.set_index(["common", "db_type", "gene"])
    common_index = my_internal_path_df.index.intersection(ref_internal_path_df.index)
    assert not common_index.empty
    my_internal_path_df_common = my_internal_path_df.loc[common_index]
    ref_internal_path_df_common = ref_internal_path_df.loc[common_index]

    known_internal_db_exceptions = [
        ["rat", "imgt", "IGHV1S62*01"],
        ["rat", "imgt", "IGHV2S1*01"],
        ["rat", "imgt", "IGHV2S12*01"],
        ["rat", "imgt", "IGHV2S13*01"],
        ["rat", "imgt", "IGHV2S18*01"],
        ["rat", "imgt", "IGHV2S30*01"],
        ["rat", "imgt", "IGHV2S35*01"],
        ["rat", "imgt", "IGHV2S54*01"],
        ["rat", "imgt", "IGHV2S61*01"],
        ["rat", "imgt", "IGHV2S63*01"],
        ["rat", "imgt", "IGHV2S75*01"],
        ["rat", "imgt", "IGHV2S78*01"],
        ["rat", "imgt", "IGHV2S8*01"],
        ["rat", "imgt", "IGHV5S10*01"],
        ["rat", "imgt", "IGHV5S11*01"],
        ["rat", "imgt", "IGHV5S13*01"],
        ["rat", "imgt", "IGHV5S14*01"],
        ["rat", "imgt", "IGHV5S23*01"],
        ["rat", "imgt", "IGHV5S47*01"],
        ["rat", "imgt", "IGHV5S54*01"],
        ["rat", "imgt", "IGHV5S8*01"],
        ["rat", "imgt", "IGHV8S18*01"],
        ["rat", "imgt", "IGHV9S3*01"],
        ["human", "imgt", "IGHV2-70*02"],
        ["human", "imgt", "IGHV2-70*03"],
        ["human", "imgt", "IGHV2-70*06"],
        ["human", "imgt", "IGHV2-70*07"],
        ["human", "imgt", "IGHV2-70*08"],
    ]
    known_internal_db_exceptions = set(map(lambda x: tuple(x), known_internal_db_exceptions))
    for index in my_internal_path_df_common.index:
        try:
            pd._testing.assert_series_equal(
                my_internal_path_df_common.loc[index],
                ref_internal_path_df_common.loc[index],
                obj=index,
            )
        except AssertionError:

            if index in known_internal_db_exceptions:
                print(index, "is known exception")
                continue
            else:
                pd._testing.assert_series_equal(
                    my_internal_path_df_common.loc[index],
                    ref_internal_path_df_common.loc[index],
                    obj=index,
                )
    return True


def test_make_igblast_reference(fixture_setup, tmpdir_factory):
    """Confirm the CLI works as expecte"""
    runner = CliRunner(echo_stdin=True)
    expected_blast_dir = fixture_setup.get_known_blast_dir_structure()
    expected_aux = fixture_setup.get_known_aux_dir_structure()
    expected_internal = fixture_setup.get_known_internal_dir_structure()
    expected_nhd = fixture_setup.get_known_nhd_dir_structure()
    tmpdir = tmpdir_factory.mktemp("igblast_dir")
    result = runner.invoke(app.make_igblast_reference, ["--outpath", tmpdir], catch_exceptions=True)
    if result.exit_code != 0:
        print(result)
        assert result.exit_code == 0
    assert os.path.exists(tmpdir)

    directories_created = glob.glob(str(tmpdir) + "/*")
    assert sorted(directories_created) == sorted([f"{tmpdir}/imgt", f"{tmpdir}/custom"])
    imgt_blast_dir = [
        i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/blastdb/*.fasta", recursive=True)
    ]
    made_diff = set(imgt_blast_dir).difference(set(expected_blast_dir))
    expected_diff = set(expected_blast_dir).difference(set(imgt_blast_dir))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a blast dbs {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected a blast db entri {sorted(expected_diff)} that was not made")
    internal = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/*.imgt", recursive=True)]
    made_diff = set(internal).difference(set(expected_internal))
    expected_diff = set(expected_internal).difference(set(internal))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a internal dbs {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected a internal db entri {sorted(expected_diff)} that was not made")
    aux = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/*.aux", recursive=True)]
    made_diff = set(aux).difference(set(expected_aux))
    expected_diff = set(expected_aux).difference(set(aux))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a aux structure {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected aux_structure entri {sorted(expected_diff)} that was not made")

    nhd = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/*.nhd", recursive=True)]
    made_diff = set(nhd).difference(set(expected_nhd))
    expected_diff = set(expected_nhd).difference(set(nhd))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a nhd structure {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected nhd entri {sorted(expected_diff)} that was not made")

    # test auxillary file building
    assert _test_auxilary_file_structure(tmpdir, fixture_setup)

    # test internal dat file
    assert _test_internal_data_file_structure(tmpdir, fixture_setup)


def test_reference_functions(tmpdir, fixture_setup):
    """Test functions that interact directly G3"""
    with pytest.raises(G3Error):
        try:
            get_database("blah")
        except G3Error as e:
            # for coverage
            e.__str__()
            raise e
    both_databases = get_loaded_database()
    imgt_database_json_species = get_species_from_database(both_databases["imgt"])
    custom_database_json_species = get_species_from_database(both_databases["custom"])
    assert imgt_database_json_species
    assert custom_database_json_species


def test_missing_makeblast_df(tmpdir, fixture_setup):
    """Test the makeblast_df function with a missing blastdb"""
    fasta = fixture_setup.get_catnap_light_nt()
    bogus_file = fixture_setup.get_card()
    from sadie.reference.blast import write_blast_db

    with pytest.raises(ValueError):
        write_blast_db(fasta, os.path.join(tmpdir, "missing.fasta"), "some_bogus_makeblastdb")
    with pytest.raises(RuntimeError):
        write_blast_db(bogus_file, os.path.join(tmpdir, "missing.fasta"))


def test_check_models(tmpdir, fixture_setup):
    """Coverage for all models including validations exceptions for bonehead entries"""
    entry = {"species": "human", "sub_species": "human", "gene": "IGHV3-10*01", "database": "custom"}
    GeneEntry(**entry)
    entry = {"species": "human", "gene": "IGHV3-10*01", "database": "custom"}
    GeneEntry(**entry)
    with pytest.raises(ValidationError):
        # Bad third postion
        entry = {"species": "human", "gene": "IGHZ3-10*01", "database": "custom"}
        GeneEntry(**entry)
    with pytest.raises(ValidationError):
        # Bad database
        entry = {"species": "human", "gene": "IGHZ3-10*01", "database": "jordan_personal_stash"}
        GeneEntry(**entry)

    with pytest.raises(ValidationError):
        # Bad third postion
        entry = {"species": "human", "gene": ["IGHZ3-10*01", "IGHZ3-20*02"], "database": "custom"}
        GeneEntries(**entry)
    with pytest.raises(ValidationError):
        # Bad database
        entry = {"species": "human", "gene": ["IGHV3-10*01", "IGHV3-20*02"], "database": "jordan_personal_stash"}
        GeneEntries(**entry)


def test_yaml(tmpdir, fixture_setup):
    """Test yaml module"""
    yaml = YamlRef()
    assert yaml.get_database_types() == {"imgt", "custom"}
    imgt_keys = [
        "se6156",
        "sa684",
        "bat64",
        "clk",
        "dog",
        "hugl18",
        "human",
        "macaque",
        "mouse",
        "rabbit",
        "rat",
        "se09",
        "se0916",
        "se16",
    ]
    custom_keys = ["cat", "dog", "macaque"]
    assert yaml.get_species_keys("imgt") == imgt_keys
    assert yaml.get_species_keys("custom") == custom_keys

    for key in imgt_keys:
        sub_key = yaml.get_sub_species("imgt", key)
        assert sub_key
        for sub in sub_key:
            genes = yaml.get_genes("imgt", key, sub)
            yaml.get_gene_segment("imgt", key, sub, "V")
            yaml.get_gene_segment("imgt", key, sub, "D")
            yaml.get_gene_segment("imgt", key, sub, "J")
            if not genes:
                raise ValueError(f"{key} {sub_key} {genes}")
    for key in custom_keys:
        sub_key = yaml.get_sub_species("custom", key)
        assert sub_key
        for sub in sub_key:
            genes = yaml.get_genes("custom", key, sub)
            yaml.get_gene_segment("custom", key, sub, "V")
            yaml.get_gene_segment("custom", key, sub, "D")
            yaml.get_gene_segment("custom", key, sub, "J")
            if not genes:
                raise ValueError(f"{key} {sub_key} {genes}")
    assert yaml.__repr__()
