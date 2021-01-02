import glob
import os
import shutil
import tempfile
import pandas as pd
import json
from click.testing import CliRunner
from pkg_resources import resource_filename
from sadie.reference import app


split_aux = [
    "/imgt/aux_db/mouse_gl.aux",
    "/imgt/aux_db/cat_gl.aux",
    "/imgt/aux_db/rabbit_gl.aux",
    "/imgt/aux_db/alpaca_gl.aux",
    "/imgt/aux_db/dog_gl.aux",
    "/imgt/aux_db/human_gl.aux",
    "/imgt/aux_db/macaque_gl.aux",
    "/imgt/aux_db/rat_gl.aux",
    "/custom/aux_db/cat_gl.aux",
    "/custom/aux_db/macaque_gl.aux",
]

split_internal = [
    "/imgt/all/Ig/internal_data/cat/cat.ndm.imgt",
    "/imgt/all/Ig/internal_data/macaque/macaque.ndm.imgt",
    "/imgt/all/Ig/internal_data/dog/dog.ndm.imgt",
    "/imgt/all/Ig/internal_data/alpaca/alpaca.ndm.imgt",
    "/imgt/all/Ig/internal_data/rat/rat.ndm.imgt",
    "/imgt/all/Ig/internal_data/mouse/mouse.ndm.imgt",
    "/imgt/all/Ig/internal_data/rabbit/rabbit.ndm.imgt",
    "/imgt/all/Ig/internal_data/human/human.ndm.imgt",
    "/custom/all/Ig/internal_data/cat/cat.ndm.imgt",
    "/custom/all/Ig/internal_data/macaque/macaque.ndm.imgt",
]
split_nhd = [
    "/imgt/all/Ig/blastdb/cat_V.nhd",
    "/imgt/all/Ig/blastdb/human_D.nhd",
    "/imgt/all/Ig/blastdb/mouse_D.nhd",
    "/imgt/all/Ig/blastdb/rat_J.nhd",
    "/imgt/all/Ig/blastdb/human_V.nhd",
    "/imgt/all/Ig/blastdb/rabbit_J.nhd",
    "/imgt/all/Ig/blastdb/mouse_V.nhd",
    "/imgt/all/Ig/blastdb/alpaca_V.nhd",
    "/imgt/all/Ig/blastdb/dog_D.nhd",
    "/imgt/all/Ig/blastdb/macaque_J.nhd",
    "/imgt/all/Ig/blastdb/dog_V.nhd",
    "/imgt/all/Ig/blastdb/alpaca_D.nhd",
    "/imgt/all/Ig/blastdb/macaque_D.nhd",
    "/imgt/all/Ig/blastdb/alpaca_J.nhd",
    "/imgt/all/Ig/blastdb/dog_J.nhd",
    "/imgt/all/Ig/blastdb/macaque_V.nhd",
    "/imgt/all/Ig/blastdb/rat_V.nhd",
    "/imgt/all/Ig/blastdb/rabbit_D.nhd",
    "/imgt/all/Ig/blastdb/cat_J.nhd",
    "/imgt/all/Ig/blastdb/mouse_J.nhd",
    "/imgt/all/Ig/blastdb/human_J.nhd",
    "/imgt/all/Ig/blastdb/rabbit_V.nhd",
    "/imgt/all/Ig/blastdb/rat_D.nhd",
    "/imgt/all/Ig/internal_data/cat/cat_V.nhd",
    "/imgt/all/Ig/internal_data/macaque/macaque_V.nhd",
    "/imgt/all/Ig/internal_data/dog/dog_V.nhd",
    "/imgt/all/Ig/internal_data/alpaca/alpaca_V.nhd",
    "/imgt/all/Ig/internal_data/rat/rat_V.nhd",
    "/imgt/all/Ig/internal_data/mouse/mouse_V.nhd",
    "/imgt/all/Ig/internal_data/rabbit/rabbit_V.nhd",
    "/imgt/all/Ig/internal_data/human/human_V.nhd",
    "/custom/all/Ig/blastdb/cat_V.nhd",
    "/custom/all/Ig/blastdb/macaque_J.nhd",
    "/custom/all/Ig/blastdb/macaque_D.nhd",
    "/custom/all/Ig/blastdb/macaque_V.nhd",
    "/custom/all/Ig/blastdb/cat_J.nhd",
    "/custom/all/Ig/internal_data/cat/cat_V.nhd",
    "/custom/all/Ig/internal_data/macaque/macaque_V.nhd",
]

known_aux_exceptions = {
    ("mouse", "imgt", "IGLJ4*01"): "igblast has wrong number of c-term remaining",
    ("rabbit", "imgt", "IGHJ3*02"): "igblast has wrong reading frame",
}


def fixture_dir(dir):
    """Helper method for test execution."""
    return resource_filename(__name__, "fixtures/{}".format(dir))


def fixture_path(dir):
    """Helper method for test execution."""
    return resource_filename(__name__, "fixtures/{}".format(dir))


def _test_auxilary_file_structure(tmpdir):
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
    igblast_aux_files = glob.glob(fixture_dir("igblast_aux") + "/**.aux")
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


def _test_internal_data_file_structure(tmpdir):
    # what we have made
    internal_path = glob.glob(f"{tmpdir}/**/*.imgt", recursive=True)
    reference_internal_path = glob.glob(fixture_dir("igblast_internal") + "/**.imgt")
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
        df.insert(1, "db_type", file.split("/")[-6])
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


def test_make_igblast_reference():
    """Confirm the CLI works as expecte"""
    runner = CliRunner(echo_stdin=True)
    with tempfile.TemporaryDirectory(suffix="igblast_dir") as tmpdir:
        result = runner.invoke(app.make_igblast_reference, ["--outpath", tmpdir], catch_exceptions=True)
        if result.exit_code != 0:
            print(result)
        assert result.exit_code == 0
        assert os.path.exists(tmpdir)

        directories_created = glob.glob(tmpdir + "/*")
        assert sorted(directories_created) == sorted([f"{tmpdir}/imgt", f"{tmpdir}/custom"])
        imgt_blast_dir = [
            i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/blastdb/*.fasta", recursive=True)
        ]
        assert sorted(imgt_blast_dir) == sorted(json.load(open(fixture_path("blast_dir.json"))))
        internal = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/*.imgt", recursive=True)]
        assert sorted(internal) == sorted(json.load(open(fixture_path("internal.json"))))
        aux = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/*.aux", recursive=True)]
        assert sorted(aux) == sorted(json.load(open(fixture_path("aux.json"))))
        nhd = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/**/*.nhd", recursive=True)]
        assert sorted(nhd) == sorted(json.load(open(fixture_path("nhd.json"))))

        # test auxillary file building
        assert _test_auxilary_file_structure(tmpdir)

        # test internal dat file
        assert _test_internal_data_file_structure(tmpdir)

    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    assert not os.path.exists(tmpdir)
