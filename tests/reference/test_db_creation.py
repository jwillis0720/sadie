"""Unit tests for database creation interface."""
import glob
import logging
import os
import shutil
import tempfile

import pandas as pd
import pytest
from pybody.reference import app
from Bio import SeqIO
from click.testing import CliRunner
from pandas.testing import assert_frame_equal
from pkg_resources import resource_filename
from sqlalchemy import create_engine

logger = logging.getLogger()


def fixture_file(file):
    """Helper method for test execution."""
    return resource_filename(__name__, "fixtures/{}".format(file))


def _parse_file_structure(path):
    file_structure = glob.glob(f"{path}/*/*/*.fasta")

    _df = []
    for file in file_structure:
        file_name = os.path.basename(file)
        receptor = file.split("/")[-2]
        species = file.split("/")[-3]
        for index in SeqIO.parse(file, "fasta"):
            _df.append(
                {
                    "file_name": file_name,
                    "receptor": receptor,
                    "species": species,
                    "name": str(index.id),
                    "seq": str(index.seq),
                }
            )
    return pd.DataFrame(_df)


@pytest.mark.skip(reason="they will ban our ip")
def test_make_fasta_files():
    """
    Confirm the CLI can pull IMGT fasta files
    """
    runner = CliRunner()
    file_structure_df = _parse_file_structure(fixture_file("test_make_fasta_files"))
    with tempfile.TemporaryDirectory() as tmpdir:
        result = runner.invoke(app.make_imgt_fasta_files, ["--outdir", tmpdir])
        assert result.exit_code == 0
        assert os.path.exists(tmpdir)
        new_file_structure = _parse_file_structure(tmpdir)

        file_structure_df = file_structure_df.sort_values(["species", "receptor", "file_name", "name"]).reset_index(
            drop=True
        )
        # file_structure_df.to_csv("old_file_structure.csv")
        new_file_structure = new_file_structure.sort_values(["species", "receptor", "file_name", "name"]).reset_index(
            drop=True
        )
        # new_file_structure.to_csv("new_file_structure.csv")
        assert_frame_equal(file_structure_df, new_file_structure)


def test_local_imgt_database():
    """test that fixture and local copy are the same database

    this can be done by not passing an output database and skipping build

    """

    runner = CliRunner()

    fixture_engine = create_engine(
        f"sqlite:////{os.path.abspath(fixture_file('test_make_imgt_database/imgt_database_fixture.db'))}"
    )
    table_names = {
        "imgt_all": ["imgt_designation", "gene", "latin", "common", "functional"],
        "imgt_unique": ["imgt_designation", "gene", "latin", "common", "functional"],
        "j_segment_imgt": ["common", "latin", "gene"],
        "v_segment_abm": ["common", "latin", "gene"],
        "v_segment_chothia": ["common", "latin", "gene"],
        "v_segment_contact": ["common", "latin", "gene"],
        "v_segment_kabat": ["common", "latin", "gene"],
        "v_segment_imgt": ["common", "latin", "gene"],
        "v_segment_scdr": ["common", "latin", "gene"],
        "v_segment_imgt_dropped": ["common", "latin", "gene"],
    }
    assert sorted(set(table_names.keys())) == sorted(set(fixture_engine.table_names()))
    result = runner.invoke(app.make_imgt_database, ["--skip-imgt", "--skip-v", "--skip-j"])
    if result.exit_code != 0:
        print(result)
    target_engine = create_engine(
        f"sqlite:////{os.path.abspath(os.path.join(os.path.dirname(app.__file__),'data/IMGT_REFERENCE.db'))}"
    )
    assert result.exit_code == 0
    # we made a new database in tmp file
    assert set(target_engine.table_names()) == set(fixture_engine.table_names())
    # fixture and target tables should be exactly the same
    for table in table_names:
        target_df = pd.read_sql(table, con=target_engine, index_col="index")
        fixture_df = pd.read_sql(table, con=fixture_engine, index_col="index")
        # explicit assert
        assert_frame_equal(target_df, fixture_df)

        # non index matching
        assert_frame_equal(
            target_df.sort_values(table_names[table]).reset_index(drop=True),
            fixture_df.sort_values(table_names[table]).reset_index(drop=True),
        )


def test_make_imgt_database():
    """Confirm the CLI to confirm we can create the imgt database"""
    runner = CliRunner()

    fixture_engine = create_engine(
        f"sqlite:////{os.path.abspath(fixture_file('test_make_imgt_database/imgt_database_fixture.db'))}"
    )
    table_names = {
        "imgt_all": ["imgt_designation", "gene", "latin", "common", "functional"],
        "imgt_unique": ["imgt_designation", "gene", "latin", "common", "functional"],
        "j_segment_imgt": ["common", "latin", "gene"],
        "v_segment_abm": ["common", "latin", "gene"],
        "v_segment_chothia": ["common", "latin", "gene"],
        "v_segment_contact": ["common", "latin", "gene"],
        "v_segment_kabat": ["common", "latin", "gene"],
        "v_segment_imgt": ["common", "latin", "gene"],
        "v_segment_scdr": ["common", "latin", "gene"],
        "v_segment_imgt_dropped": ["common", "latin", "gene"],
    }
    assert set(table_names.keys()) == set(fixture_engine.table_names())
    with tempfile.NamedTemporaryFile(dir=".", suffix=".db", delete=True) as tmpout:
        result = runner.invoke(app.make_imgt_database, ["--outpath", tmpout.name])

        # create a brand new database
        target_engine = create_engine(f"sqlite:////{os.path.abspath(tmpout.name)}")
        if result.exit_code != 0:
            # show me the exception
            print(result)
            raise result.exception
        assert result.exit_code == 0
        # we made a new database in tmp file
        assert set(target_engine.table_names()) == set(fixture_engine.table_names())
        # fixture and target tables  may not be exactly the same depending on OS. I have no clue why
        for table in table_names:
            target_df = pd.read_sql(table, con=target_engine, index_col="index")
            fixture_df = pd.read_sql(table, con=fixture_engine, index_col="index")
            if "file" in target_df.columns:
                target_df = target_df.drop("file", axis=1)
                fixture_df = fixture_df.drop("file", axis=1)
            # # explicit assert
            # assert_frame_equal(target_df, fixture_df)
            # non index matching
            assert_frame_equal(
                target_df.sort_values(table_names[table]).reset_index(drop=True),
                fixture_df.sort_values(table_names[table]).reset_index(drop=True),
            )
