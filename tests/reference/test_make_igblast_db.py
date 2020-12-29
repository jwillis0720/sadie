import glob
import os
import shutil
import tempfile
import pandas as pd
from click.testing import CliRunner
from pkg_resources import resource_filename
from sadie.reference import app

split_fastas = [
    "/imgt/Ig/blastdb/dog_V.fasta",
    "/imgt/Ig/blastdb/cat_V.fasta",
    "/imgt/Ig/blastdb/rabbit_D.fasta",
    "/imgt/Ig/blastdb/alpaca_V.fasta",
    "/imgt/Ig/blastdb/mouse_D.fasta",
    "/imgt/Ig/blastdb/macaque_J.fasta",
    "/imgt/Ig/blastdb/rat_D.fasta",
    "/imgt/Ig/blastdb/human_V.fasta",
    "/imgt/Ig/blastdb/human_J.fasta",
    "/imgt/Ig/blastdb/alpaca_J.fasta",
    "/imgt/Ig/blastdb/macaque_V.fasta",
    "/imgt/Ig/blastdb/cat_J.fasta",
    "/imgt/Ig/blastdb/dog_J.fasta",
    "/imgt/Ig/blastdb/dog_D.fasta",
    "/imgt/Ig/blastdb/rabbit_V.fasta",
    "/imgt/Ig/blastdb/rat_V.fasta",
    "/imgt/Ig/blastdb/human_D.fasta",
    "/imgt/Ig/blastdb/alpaca_D.fasta",
    "/imgt/Ig/blastdb/mouse_V.fasta",
    "/imgt/Ig/blastdb/mouse_J.fasta",
    "/imgt/Ig/blastdb/macaque_D.fasta",
    "/imgt/Ig/blastdb/rat_J.fasta",
    "/imgt/Ig/blastdb/rabbit_J.fasta",
    "/imgt/Ig/internal_data/cat/cat_V.fasta",
    "/imgt/Ig/internal_data/macaque/macaque_V.fasta",
    "/imgt/Ig/internal_data/dog/dog_V.fasta",
    "/imgt/Ig/internal_data/alpaca/alpaca_V.fasta",
    "/imgt/Ig/internal_data/rat/rat_V.fasta",
    "/imgt/Ig/internal_data/mouse/mouse_V.fasta",
    "/imgt/Ig/internal_data/rabbit/rabbit_V.fasta",
    "/imgt/Ig/internal_data/human/human_V.fasta",
    "/custom/Ig/blastdb/cat_V.fasta",
    "/custom/Ig/blastdb/macaque_J.fasta",
    "/custom/Ig/blastdb/macaque_V.fasta",
    "/custom/Ig/blastdb/cat_J.fasta",
    "/custom/Ig/blastdb/macaque_D.fasta",
    "/custom/Ig/internal_data/cat/cat_V.fasta",
    "/custom/Ig/internal_data/macaque/macaque_V.fasta",
]

split_aux = [
    "/imgt/Ig/aux_db/mouse_gl.aux",
    "/imgt/Ig/aux_db/cat_gl.aux",
    "/imgt/Ig/aux_db/rabbit_gl.aux",
    "/imgt/Ig/aux_db/alpaca_gl.aux",
    "/imgt/Ig/aux_db/dog_gl.aux",
    "/imgt/Ig/aux_db/human_gl.aux",
    "/imgt/Ig/aux_db/macaque_gl.aux",
    "/imgt/Ig/aux_db/rat_gl.aux",
    "/custom/Ig/aux_db/cat_gl.aux",
    "/custom/Ig/aux_db/macaque_gl.aux",
]

split_internal = [
    "/imgt/Ig/internal_data/cat/cat.ndm.imgt",
    "/imgt/Ig/internal_data/macaque/macaque.ndm.imgt",
    "/imgt/Ig/internal_data/dog/dog.ndm.imgt",
    "/imgt/Ig/internal_data/alpaca/alpaca.ndm.imgt",
    "/imgt/Ig/internal_data/rat/rat.ndm.imgt",
    "/imgt/Ig/internal_data/mouse/mouse.ndm.imgt",
    "/imgt/Ig/internal_data/rabbit/rabbit.ndm.imgt",
    "/imgt/Ig/internal_data/human/human.ndm.imgt",
    "/custom/Ig/internal_data/cat/cat.ndm.imgt",
    "/custom/Ig/internal_data/macaque/macaque.ndm.imgt",
]
split_nhd = [
    "/imgt/Ig/blastdb/cat_V.nhd",
    "/imgt/Ig/blastdb/human_D.nhd",
    "/imgt/Ig/blastdb/mouse_D.nhd",
    "/imgt/Ig/blastdb/rat_J.nhd",
    "/imgt/Ig/blastdb/human_V.nhd",
    "/imgt/Ig/blastdb/rabbit_J.nhd",
    "/imgt/Ig/blastdb/mouse_V.nhd",
    "/imgt/Ig/blastdb/alpaca_V.nhd",
    "/imgt/Ig/blastdb/dog_D.nhd",
    "/imgt/Ig/blastdb/macaque_J.nhd",
    "/imgt/Ig/blastdb/dog_V.nhd",
    "/imgt/Ig/blastdb/alpaca_D.nhd",
    "/imgt/Ig/blastdb/macaque_D.nhd",
    "/imgt/Ig/blastdb/alpaca_J.nhd",
    "/imgt/Ig/blastdb/dog_J.nhd",
    "/imgt/Ig/blastdb/macaque_V.nhd",
    "/imgt/Ig/blastdb/rat_V.nhd",
    "/imgt/Ig/blastdb/rabbit_D.nhd",
    "/imgt/Ig/blastdb/cat_J.nhd",
    "/imgt/Ig/blastdb/mouse_J.nhd",
    "/imgt/Ig/blastdb/human_J.nhd",
    "/imgt/Ig/blastdb/rabbit_V.nhd",
    "/imgt/Ig/blastdb/rat_D.nhd",
    "/imgt/Ig/internal_data/cat/cat_V.nhd",
    "/imgt/Ig/internal_data/macaque/macaque_V.nhd",
    "/imgt/Ig/internal_data/dog/dog_V.nhd",
    "/imgt/Ig/internal_data/alpaca/alpaca_V.nhd",
    "/imgt/Ig/internal_data/rat/rat_V.nhd",
    "/imgt/Ig/internal_data/mouse/mouse_V.nhd",
    "/imgt/Ig/internal_data/rabbit/rabbit_V.nhd",
    "/imgt/Ig/internal_data/human/human_V.nhd",
    "/custom/Ig/blastdb/cat_V.nhd",
    "/custom/Ig/blastdb/macaque_J.nhd",
    "/custom/Ig/blastdb/macaque_D.nhd",
    "/custom/Ig/blastdb/macaque_V.nhd",
    "/custom/Ig/blastdb/cat_J.nhd",
    "/custom/Ig/internal_data/cat/cat_V.nhd",
    "/custom/Ig/internal_data/macaque/macaque_V.nhd",
]

known_aux_exceptions = {
    ("mouse", "imgt", "IGLJ4*01"): "igblast has wrong number of c-term remaining",
    ("rabbit", "imgt", "IGHJ3*02"): "igblast has wrong reading frame",
}


def fixture_dir(dir):
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
        df.insert(1, "db_type", file.split("/")[-4])
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
            i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/*/**/*.fasta", recursive=True)
        ]
        assert sorted(imgt_blast_dir) == sorted(split_fastas)
        internal = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/*/**/*.imgt", recursive=True)]
        assert sorted(internal) == sorted(split_internal)
        aux = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/*/**/*.aux", recursive=True)]
        assert sorted(aux) == sorted(split_aux)
        nhd = [i.split(os.path.basename(tmpdir))[-1] for i in glob.glob(f"{tmpdir}/*/**/*.nhd", recursive=True)]
        assert sorted(nhd) == sorted(split_nhd)

        assert _test_auxilary_file_structure(tmpdir)

    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    assert not os.path.exists(tmpdir)
