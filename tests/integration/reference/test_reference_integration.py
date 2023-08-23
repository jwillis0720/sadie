import glob
import os
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner

from sadie import app
from tests.conftest import SadieFixture


def _test_internal_data_file_structure(tmpdir: Path, fixture_setup: SadieFixture):
    """test pipeline and the expected and internal data. Will skip expected missing"""
    internal_path = glob.glob(f"{tmpdir}/Ig/internal_data/**/*.imgt", recursive=True)
    reference_internal_path = fixture_setup.get_internal_files()
    my_internal_path_df = []

    # read each .imgt file and turn it into a dataframe
    _names: Any = [
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
    ]
    for file in internal_path:
        df = pd.read_csv(
            file,
            skip_blank_lines=True,
            delimiter="\t",
            names=_names,
            header=None,
        )
        df.insert(0, "name", os.path.basename(file).split(".ndm")[0])
        # df.insert(1, "db_type", "imgt")  # pylint: disable=maybe-no-member
        my_internal_path_df.append(df)

    my_internal_path_df = pd.concat(my_internal_path_df).reset_index(drop=True).groupby(["name", "gene"]).head(1)
    my_internal_path_df = my_internal_path_df.set_index(["name", "gene"])

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
        df.insert(0, "name", common_name)
        ref_internal_path_df.append(df)

    ref_internal_path_df = pd.concat(ref_internal_path_df).reset_index(drop=True)
    # only get first occurance
    ref_internal_path_df = ref_internal_path_df.groupby(["name", "gene"]).head(1)
    ref_internal_path_df = ref_internal_path_df.set_index(["name", "gene"])

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


def _test_auxilary_file_structure(tmpdir: Path, fixture_setup: SadieFixture) -> bool:
    """test pipeline and the expected and aux data. Will skip expected missing"""
    my_aux = []
    generated_aux_path_files = glob.glob(f"{tmpdir}/aux_db/**/*.aux", recursive=True)
    for file in generated_aux_path_files:
        df = pd.read_csv(
            file,
            skip_blank_lines=True,
            delimiter="\t",
            header=None,
            names=["gene", "reading_frame", "segment", "cdr3_end", "left_over"],
        )
        df.insert(0, "name", os.path.basename(file).split("_")[0])
        my_aux.append(df)

    # make the aux file structure into a dataframe
    my_aux = pd.concat(my_aux).reset_index(drop=True).set_index(["name", "gene"])

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
        df.insert(0, "name", os.path.basename(file).split("_")[0])
        # df.insert(1, "db_type", "imgt")  # pylint: disable=maybe-no-member
        igblast_aux.append(df)

    # here is what we should get
    igblast_aux = pd.concat(igblast_aux).reset_index(drop=True).set_index(["name", "gene"])

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


def test_make_igblast_reference(fixture_setup: SadieFixture, tmp_path_factory: pytest.TempPathFactory):
    """Confirm the CLI works as expected This runs the entire generation pipeline that ships with SADIE and checks that the file structure is exactly the same"""
    runner = CliRunner(echo_stdin=True)

    # these are the expected file structures
    expected_blast_dir = fixture_setup.get_known_blast_dir_structure()
    expected_aux = fixture_setup.get_known_aux_dir_structure()
    expected_internal = fixture_setup.get_known_internal_dir_structure()
    # make a hierarchy of directories
    tmpdir = tmp_path_factory.mktemp("igblast_dir")

    # run the entire pipeline via CLICK cli
    result = runner.invoke(app.make_igblast_reference, ["--outpath", tmpdir], catch_exceptions=True)
    if result.exit_code != 0:
        print(result)
        assert result.exit_code == 0

    # was the file actually output?
    assert os.path.exists(tmpdir)

    # assert we made an imgt and custom directory, but still don't know if anything is in it
    directories_created = glob.glob(str(tmpdir) + "/*")
    assert sorted(directories_created) == sorted([f"{tmpdir}/aux_db", f"{tmpdir}/Ig"])

    # for the blast directory, let's check if all the fasta files are there
    blast_files = glob.glob(f"{tmpdir}/Ig/blastdb/**/*.*", recursive=True)
    blast_files = [i.split(os.path.basename(tmpdir) + "/")[-1] for i in blast_files]

    # even though this could be done with a symmetric diff, using diff tells us which on is missing, expected or made
    made_diff = set(blast_files).difference(set(expected_blast_dir))
    expected_diff = set(expected_blast_dir).difference(set(blast_files))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a blast dbs {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected a blast db entri {sorted(expected_diff)} that was not made")

    # do the same with .aux files in the aux directory. This doesn't check content, just that the files are there
    aux = [i.split(os.path.basename(tmpdir) + "/")[-1] for i in glob.glob(f"{tmpdir}/**/*.aux", recursive=True)]
    made_diff = set(aux).difference(set(expected_aux))
    expected_diff = set(expected_aux).difference(set(aux))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a aux structure {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected aux_structure entri {sorted(expected_diff)} that was not made")
    # do the same with the .imgt files in the internal directory. This doesn't check content, just that the files are there
    # ** is the species
    internal = [
        i.split(os.path.basename(tmpdir) + "/")[-1]
        for i in glob.glob(f"{tmpdir}/Ig/internal_data/**/*", recursive=True)
        if not Path(i).is_dir()
    ]
    made_diff = set(internal).difference(set(expected_internal))
    expected_diff = set(expected_internal).difference(set(internal))
    if made_diff or expected_diff:
        if made_diff:
            raise AssertionError(f"We made a internal dbs {sorted(made_diff)} that was not expected")
        if expected_diff:
            raise AssertionError(f"We expected a internal db entri {sorted(expected_diff)} that was not made")

    # these next two functions actually check content of internal and aux that come shipped with IgBlast. This ensures we are using the same file
    # test auxillary file content
    assert _test_auxilary_file_structure(tmpdir, fixture_setup)

    # test internal dat file content
    assert _test_internal_data_file_structure(tmpdir, fixture_setup)
