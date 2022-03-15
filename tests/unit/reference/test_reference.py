import glob
import os
from pathlib import Path
from pydantic import ValidationError
import pandas as pd
import pytest
from click.testing import CliRunner
from sadie import app
from sadie.reference import Reference
from sadie.reference.reference import G3Error, get_database, get_loaded_database, get_species_from_database
from sadie.reference.models import GeneEntry, GeneEntries
from sadie.reference.yaml import YamlRef
from sadie.reference.reference import _write_out_fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _test_internal_data_file_structure(tmpdir, fixture_setup):
    """test pipeline and the expected and internal data. Will skip expected missing"""
    internal_path = glob.glob(f"{tmpdir}/imgt/**/*.imgt", recursive=True)
    reference_internal_path = fixture_setup.get_internal_files()
    my_internal_path_df = []

    # read each .imgt file and turn it into a dataframe
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

    # read in expected internal db files
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

        # change rhesus to macaque
        if common_name == "rhesus_monkey":
            common_name = "macaque"
        df.insert(0, "common", common_name)
        df.insert(1, "db_type", "imgt")
        ref_internal_path_df.append(df)

    ref_internal_path_df = pd.concat(ref_internal_path_df).reset_index(drop=True)
    # only get first occurance
    ref_internal_path_df = ref_internal_path_df.groupby(["common", "db_type", "gene"]).head(1)
    ref_internal_path_df = ref_internal_path_df.set_index(["common", "db_type", "gene"])

    # get common index
    common_index = my_internal_path_df.index.intersection(ref_internal_path_df.index)
    assert not common_index.empty
    my_internal_path_df_common = my_internal_path_df.loc[common_index]
    ref_internal_path_df_common = ref_internal_path_df.loc[common_index]
    known_internal_db_exceptions = fixture_setup.get_internal_db_excetions()
    known_internal_db_exceptions = set(map(lambda x: tuple(x), known_internal_db_exceptions))

    # go through each file as a series and compare. This gives more verbose output
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


def _test_auxilary_file_structure(tmpdir, fixture_setup):
    """test pipeline and the expected and aux data. Will skip expected missing"""
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

    # make the aux file structure into a dataframe
    my_aux = pd.concat(my_aux).reset_index(drop=True).set_index(["common", "db_type", "gene"])

    # get what it should look like and it's conent
    igblast_aux = []
    igblast_aux_files = fixture_setup.get_aux_files()

    # go through every file and add it to the dataframe
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

    # here is what we should get
    igblast_aux = pd.concat(igblast_aux).reset_index(drop=True).set_index(["common", "db_type", "gene"])

    # make sure there is proper intersection
    common_index = my_aux.index.intersection(igblast_aux.index)
    my_aux_common_index = my_aux.loc[common_index]
    igblast_common_index = igblast_aux.loc[common_index]

    # here are known exceptions
    known_aux_exceptions = fixture_setup.get_aux_exceptions()

    # go through index one by one and see that the series are equal. This is the best to see specifics of each series
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


# begin tests
def test_reference_class(tmpdir_factory):
    """Test if we can JIT reference class."""
    ref_class = Reference()
    ref_class.add_gene({"species": "human", "gene": "IGHV1-69*01", "database": "imgt"})
    ref_class.add_gene({"species": "human", "gene": "IGHD3-3*01", "database": "imgt"})
    ref_class.add_gene({"species": "human", "gene": "IGHJ6*01", "database": "imgt"})
    with pytest.raises(G3Error):
        ref_class.add_gene({"species": "human", "gene": "IGHV111-69*01", "database": "imgt"})
    assert len(ref_class.get_dataframe()) == 3
    output = tmpdir_factory.mktemp("test_reference_class")
    ref_class.make_airr_database(output)


def test_private_methods(tmpdir_factory):
    # query = "https://g3.jordanrwillis.com/api/v1/genes?source=imgt&common=human&gene=IGHV1-69%2A01"
    seq = "AAAAA"
    file = Path(tmpdir_factory.mktemp("test_private_methods") + "/test.fasta")
    _write_out_fasta([SeqRecord(Seq(seq), id="test", name="test")], file)


def test_creation_from_empty_reference(tmpdir_factory):
    """Test when we create reference without pasing any data, it uses yaml"""
    tmpdir = tmpdir_factory.mktemp("test_creation_from_empty_reference")
    ref_class = Reference()
    ref_class.make_airr_database(Path(tmpdir))


def test_load_ref_from_df(fixture_setup, tmpdir_factory):
    """Test if we can statically load a reference csv"""
    ref_class = Reference.read_file(fixture_setup.get_reference_dataset_csv())
    assert ref_class.data
    outpath = tmpdir_factory.mktemp("test_load_ref_from_df")
    outfile = pd.read_csv(fixture_setup.get_reference_dataset_csv(), index_col=0)
    outfile.to_csv(outpath + "/test.csv")
    outfile.to_feather(outpath + "/test.feather")
    outfile.to_json(outpath + "/test.json", orient="records")
    ref_class = Reference.read_file(outpath + "/test.json", type="json")
    ref_class = Reference.read_file(outpath + "/test.csv")
    ref_class = Reference.read_file(outpath + "/test.feather", type="feather")
    with pytest.raises(ValueError):
        ref_class = Reference.read_file(outpath + "/test.feather", type="oinga")


def test_make_reference_class_from_yaml():
    """Test reference class."""
    ref_class = Reference.parse_yaml()
    assert isinstance(ref_class, Reference)
    ref_class_data = ref_class.get_dataframe()
    assert isinstance(ref_class_data, pd.DataFrame)


def test_make_igblast_reference(fixture_setup, tmpdir_factory):
    """Confirm the CLI works as expected This runs the entire generation pipeline that ships with SADIE and checks that the file structure is exactly the same"""
    runner = CliRunner(echo_stdin=True)

    # these are the expected file structures
    expected_blast_dir = fixture_setup.get_known_blast_dir_structure()
    expected_aux = fixture_setup.get_known_aux_dir_structure()
    expected_internal = fixture_setup.get_known_internal_dir_structure()
    expected_nhd = fixture_setup.get_known_nhd_dir_structure()

    # make a hierarchy of directories
    tmpdir = tmpdir_factory.mktemp("igblast_dir")

    # run the entire pipeline via CLICK cli
    result = runner.invoke(app.make_igblast_reference, ["--outpath", tmpdir], catch_exceptions=True)
    if result.exit_code != 0:
        print(result)
        assert result.exit_code == 0

    # was the file actually output?
    assert os.path.exists(tmpdir)

    # assert we made an imgt and custom directory, but still don't know if anything is in it
    directories_created = glob.glob(str(tmpdir) + "/*")
    assert sorted(directories_created) == sorted([f"{tmpdir}/imgt", f"{tmpdir}/custom"])

    # for the blast directory, let's check if all the fasta files are there
    imgt_blast_dir = [
        i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/blastdb/*.fasta", recursive=True)
    ]

    # even though this could be done with a symmetric diff, using diff tells us which on is missing, expected or made
    made_diff = set(imgt_blast_dir).difference(expected_blast_dir)
    expected_diff = set(expected_blast_dir).difference(set(imgt_blast_dir))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a blast dbs {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected a blast db entri {sorted(expected_diff)} that was not made")

    # do the same with the .imgt files in the internal directory. This doesn't check content, just that the files are there
    internal = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/*.imgt", recursive=True)]
    made_diff = set(internal).difference(set(expected_internal))
    expected_diff = set(expected_internal).difference(set(internal))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a internal dbs {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected a internal db entri {sorted(expected_diff)} that was not made")

    # do the same with .aux files in the aux directory. This doesn't check content, just that the files are there
    aux = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/*.aux", recursive=True)]
    made_diff = set(aux).difference(set(expected_aux))
    expected_diff = set(expected_aux).difference(set(aux))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a aux structure {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected aux_structure entri {sorted(expected_diff)} that was not made")

    # finally do the same with nhd files, the nhd files are made by makeblastdb, and are binary so we still aren't checking the conent
    nhd = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/*.nhd", recursive=True)]
    made_diff = set(nhd).difference(set(expected_nhd))
    expected_diff = set(expected_nhd).difference(set(nhd))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a nhd structure {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected nhd entri {sorted(expected_diff)} that was not made")

    # these next two functions actually check content of internal and aux
    # test auxillary file content
    assert _test_auxilary_file_structure(tmpdir, fixture_setup)

    # test internal dat file content
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
    ref = Reference()
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

    with pytest.raises(ValueError):
        ref._get_gene(GeneEntries)
    with pytest.raises(ValueError):
        ref._get_genes(GeneEntry)
    entry = {"species": "mouse", "sub_species": "mouse", "gene": "IGHV2-6-3*01", "database": "imgt"}
    ref._get_gene(GeneEntry(**entry))


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


def test_G3_errors():
    """Test G3 errors"""
    with pytest.raises(G3Error):
        Reference(endpoint="https://mock.codes/202")
