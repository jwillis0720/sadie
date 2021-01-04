"""Unit tests for analysis interface."""
import os
from click.testing import CliRunner
import logging
import glob
import pytest
from itertools import product
import tempfile
import pandas as pd
from sadie.airr import Airr, BadDataSet, AirrTable
from sadie.airr import app

logger = logging.getLogger()


def get_file(file):
    """Helper method for test execution."""
    _file = os.path.join(os.path.abspath(os.path.dirname(__file__)), f"fixtures/{file}")
    if not os.path.exists(_file):
        raise FileNotFoundError(_file)
    return _file


def test_airr_init():
    # show each of these can run
    for species in ["human", "mouse", "rat", "dog"]:
        air_api = Airr(species)
        air_api.get_available_datasets()
        assert isinstance(air_api, Airr)
    # show we can catch bad species inputs
    with pytest.raises(BadDataSet):
        for species in ["robot", "scarecrow"]:
            air_api = Airr(species)


def test_airr_single_sequence():
    pg9_seq = """
        CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGT
        CCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGT
        GGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCC
        AGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATT
        TTTGTGTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTA
        TTATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC""".replace(
        "\n", ""
    )
    air_api = Airr("human")
    airr_table = air_api.run_single("PG9", pg9_seq)
    airr_entry = airr_table.iloc[0]
    cdr3_ = airr_entry["cdr3_aa"]
    cdr2_ = airr_entry["cdr2_aa"]
    cdr1_ = airr_entry["cdr1_aa"]
    fw1_ = airr_entry["fwr1_aa"]
    fw2_ = airr_entry["fwr2_aa"]
    fw3_ = airr_entry["fwr3_aa"]
    fw4_ = airr_entry["fwr4_aa"]
    assert fw1_ == "QRLVESGGGVVQPGSSLRLSCAAS"
    assert fw2_ == "MHWVRQAPGQGLEWVAF"
    assert fw3_ == "YHADSVWGRLSISRDNSKDTLYLQMNSLRVEDTATYFC"
    assert fw4_ == "WGKGTTVTVSS"
    assert cdr1_ == "GFDFSRQG"
    assert cdr2_ == "IKYDGSEK"
    assert cdr3_ == "VREAGGPDYRNGYNYYDFYDGYYNYHYMDV"


def test_airr_from_dataframe():
    """Test we can pass a dataframe to runtime"""
    dog_df = pd.read_csv(get_file("airr_tables/dog_igh.csv.gz"), index_col=0)
    airr_api = Airr("dog")
    unjoined_df = airr_api.run_dataframe(dog_df, "sequence_id", "sequence")
    assert isinstance(unjoined_df, AirrTable)
    joined_df = airr_api.run_dataframe(dog_df, "sequence_id", "sequence", return_join=True)
    assert isinstance(joined_df, pd.DataFrame)


def _run_cli(args, tmpfile):
    runner = CliRunner()
    result = runner.invoke(app.run_airr, args)
    assert result.exit_code == 0
    assert os.path.exists(tmpfile.name)
    return True


def test_cli():
    """Confirm the CLI works as expecte"""
    quereies = glob.glob(get_file("fasta_inputs/") + "*")
    species = ["dog", "rat", "human", "mouse", "macaque"]
    ft = ["csv", "json"]
    functions = ["--all", "--functional"]
    products = product(quereies, species, ["imgt"], ft, functions)

    with tempfile.NamedTemporaryFile() as tmpfile:
        for p_tuple in products:
            cli_input = [
                "--query",
                p_tuple[0],
                "-o",
                tmpfile.name,
                "-s",
                p_tuple[1],
                "--database",
                p_tuple[2],
                "-f",
                p_tuple[3],
                p_tuple[4],
            ]
            logger.debug(f"CLI input {' '.join(cli_input)}")
            assert _run_cli(cli_input, tmpfile)
    quereies = glob.glob(get_file("fasta_inputs/") + "*")
    species = ["cat", "macaque"]
    ft = ["csv", "json"]
    functions = ["--all", "--functional"]
    products = product(quereies, species, ["custom"], ft, functions)
    with tempfile.NamedTemporaryFile() as tmpfile:
        for p_tuple in products:
            cli_input = [
                "--query",
                p_tuple[0],
                "-o",
                tmpfile.name,
                "-s",
                p_tuple[1],
                "--database",
                p_tuple[2],
                "-f",
                p_tuple[3],
                p_tuple[4],
            ]
            logger.debug(f"CLI input {' '.join(cli_input)}")
            assert _run_cli(cli_input, tmpfile)
    quereies = [get_file("fasta_inputs/scfv.fasta")]
    species = ["human"]
    ft = ["csv", "json"]
    functions = ["--all", "--functional"]
    products = product(quereies, species, ["imgt"], ft, functions)
    with tempfile.NamedTemporaryFile() as tmpfile:
        for p_tuple in products:
            cli_input = [
                "--query",
                p_tuple[0],
                "-o",
                tmpfile.name,
                "-s",
                p_tuple[1],
                "--database",
                p_tuple[2],
                "-f",
                p_tuple[3],
                p_tuple[4],
                "--scfv",
            ]
            logger.debug(f"CLI input {' '.join(cli_input)}")
            assert _run_cli(cli_input, tmpfile)
    assert not os.path.exists(tmpfile.name)
