"""Unit tests for analysis interface."""
import os
from pkg_resources import resource_filename
from click.testing import CliRunner
import logging
import glob
import pytest
import tempfile
from sadie.airr import Airr, BadSpecies
from sadie.airr import app

logger = logging.getLogger()


def fixture_file(file):
    """Helper method for test execution."""
    return resource_filename(__name__, "fixtures/{}".format(file))


def get_fasta_inputs(wildcard):
    fasta_input_dir = fixture_file("fasta_inputs")
    return glob.glob(fasta_input_dir + "/*/{wildcard}".format(wildcard=wildcard))


def test_airr_init():
    # show each of these can run
    for species in ["human", "mouse", "rat", "cat", "dog"]:
        air_api = Airr(species)
        # logger.debug(air_api.__repr__())
        assert isinstance(air_api, Airr)
    # show we can catch bad species inputs
    with pytest.raises(BadSpecies):
        for species in ["robot", "scarecrow"]:
            air_api = Airr(species)
            # logger.debug(air_api.__repr__())


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


def test_cli():
    """Confirm the CLI works as expecte"""
    runner = CliRunner()
    fastas = get_fasta_inputs("*")
    for query_file in fastas:
        with tempfile.NamedTemporaryFile() as tmpfile:
            logger.debug(tmpfile.name)
            result = runner.invoke(app.run_airr, ["--query", query_file, "-o", tmpfile.name])
            assert result.exit_code == 0
            assert os.path.exists(tmpfile.name)
    assert not os.path.exists(tmpfile.name)
    for query_file in fastas:
        with tempfile.NamedTemporaryFile() as tmpfile:
            logger.debug(tmpfile.name)
            result = runner.invoke(app.run_airr, ["--query", query_file, "-s", "cat", "-o", tmpfile.name])
            assert result.exit_code == 0
            assert os.path.exists(tmpfile.name)
        assert not os.path.exists(tmpfile.name)
