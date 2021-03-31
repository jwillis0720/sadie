"""Unit tests for analysis interface."""
import glob
import logging
import os
import tempfile

import gzip
from itertools import product

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner
from numpy import isnan
from sadie.airr import Airr, AirrTable, BadDataSet, BadRequstedFileType, GermlineData, ScfvAirrTable, app

logger = logging.getLogger()


def get_file(file):
    """Helper method for test execution."""
    _file = os.path.join(os.path.abspath(os.path.dirname(__file__)), f"fixtures/{file}")
    if not os.path.exists(_file):
        raise FileNotFoundError(_file)
    return _file


def test_germline_init():
    """Germline Data in airr module"""
    gd = GermlineData("human")
    with pytest.raises(FileNotFoundError):
        gd.base_dir = "/non/existant/path"
    with pytest.raises(FileNotFoundError):
        gd.blast_dir = "/non/existant/path"
    with pytest.raises(FileNotFoundError):
        gd.v_gene_dir = "/non/existant/path"
    with pytest.raises(FileNotFoundError):
        gd.j_gene_dir = "/non/existant/path"
    with pytest.raises(FileNotFoundError):
        gd.igdata = "/non/existant/path"
    with pytest.raises(FileNotFoundError):
        gd.aux_path = "/non/existant/path"
    with pytest.warns(UserWarning):
        gd.d_gene_dir = "/non/existant/path"


def test_airr_init():
    # show each of these can run
    for species in ["human", "mouse", "rat", "dog"]:
        air_api = Airr(species)
        air_api.get_available_datasets()
        assert isinstance(air_api, Airr)
    # show we can catch bad species inputs
    with pytest.raises(BadDataSet) as execinfo:
        air_api = Airr("robot")
    assert execinfo.value.__str__()

    # have an already created directory
    air_api = Airr("human", temp_directory="a_new_tempfile")
    os.rmdir("a_new_tempfile")

    assert air_api.__repr__() == air_api.__str__()


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
    v_mutation = airr_entry["v_mutation"]
    d_mutation = airr_entry["d_mutation"]
    j_mutation = airr_entry["j_mutation"]

    assert fw1_ == "QRLVESGGGVVQPGSSLRLSCAAS"
    assert fw2_ == "MHWVRQAPGQGLEWVAF"
    assert fw3_ == "YHADSVWGRLSISRDNSKDTLYLQMNSLRVEDTATYFC"
    assert fw4_ == "WGKGTTVTVSS"
    assert cdr1_ == "GFDFSRQG"
    assert cdr2_ == "IKYDGSEK"
    assert cdr3_ == "VREAGGPDYRNGYNYYDFYDGYYNYHYMDV"
    assert v_mutation == 14.0
    assert d_mutation == 17.875
    assert j_mutation == 11.3125
    with pytest.raises(TypeError):
        # id must be str
        airr_table = air_api.run_single(9, pg9_seq)

    # super edge case, two heavy or two light scfv
    pg9_seq_heavy = """
    CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGT
    CCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGT
    GGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCC
    AGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATT
    TTTGTGTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTA
    TTATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC""".replace(
        "\n", ""
    )
    pg9_seq_light = """
    CAGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCAATGGAACCAGCAA
    TGATGTTGGTGGCTATGAATCTGTCTCCTGGTACCAACAACATCCCGGCAAAGCCCCCAAAGTCGTGATTTATGATGTCA
    GTAAACGGCCCTCAGGGGTTTCTAATCGCTTCTCTGGCTCCAAGTCCGGCAACACGGCCTCCCTGACCATCTCTGGGCTC
    CAGGCTGAGGACGAGGGTGACTATTACTGCAAGTCTCTGACAAGCACGAGACGTCGGGTTTTCGGCACTGGGACCAAGCT
    GACCGTTCTA""".replace(
        "\n", ""
    )
    air_api.run_single("two_heavy", pg9_seq_heavy + pg9_seq_heavy, scfv=True)
    air_api.run_single("two_light", pg9_seq_light + pg9_seq_light, scfv=True)


def test_run_multiple():
    airr = Airr("human")
    pg9_seq = """
    CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGT
    CCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGT
    GGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCC
    AGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATT
    TTTGTGTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTA
    TTATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC""".replace(
        "\n", ""
    )
    seq = Seq(pg9_seq)
    seq_record = SeqRecord(seq=seq)
    airr.run_multiple([seq_record, seq_record])
    with pytest.raises(TypeError):
        airr.run_multiple(seq_record)


def test_run_multiple_scfv():
    scfv = get_file("fastq_inputs/long_scfv.fastq.gz")
    airr = Airr("human")
    list_to_run = list(SeqIO.parse(gzip.open(scfv, "rt"), "fastq"))
    results = airr.run_multiple(list_to_run, scfv=True)
    assert isinstance(results, ScfvAirrTable)


def test_airr_from_dataframe():
    """Test we can pass a dataframe to runtime"""
    dog_df = pd.read_csv(get_file("airr_tables/dog_igh.csv.gz"), index_col=0)
    airr_api = Airr("dog")
    unjoined_df = airr_api.run_dataframe(dog_df, "sequence_id", "sequence")
    assert isinstance(unjoined_df, AirrTable)
    joined_df = airr_api.run_dataframe(dog_df, "sequence_id", "sequence", return_join=True)
    assert isinstance(joined_df, pd.DataFrame)


def test_airr_from_file():
    """Test we can pass a dataframe to runtime"""
    f = get_file("fasta_inputs/PG9_H_multiple.fasta.gz")
    airr_api = Airr("human")
    result = airr_api.run_file(f)
    assert isinstance(result, AirrTable)

    # see if we can take in paths
    result = airr_api.run_file(Path(f))
    assert isinstance(result, AirrTable)

    f = get_file("card.png")
    with pytest.raises(BadRequstedFileType) as execinfo:
        airr_api.run_file(f)
    assert execinfo.value.__str__()

    f = get_file("fasta_inputs/scfv.fasta.gz")
    airr_api = Airr("human", temp_directory="a_non_existant_directory")
    airr_api.run_file(f, scfv=True)


def test_adaptable_penalty():
    test_sequence = """
    GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACC
    ATCACTTGCCAGGCGAGTCAGGACATTAGCAACTATTTAAATTGGTATCAGCAGAAACCA
    GGGAAAGCCCCTAAGCTCCTGATCTACGATGCATCCAATTTGGAAACAGGGGTCCCATCA
    AGGTTCAGTGGAAGTGGATCTGGGACAGATTTTACTTTCACCATCAGCAGCCTGCAGCCT
    GAAGATATTGCAACATATTACTGTCAACAGTATGATAATTTCGGCGGAGGGACCAAGGTG
    GACATCAAAC""".replace(
        "\n", ""
    )

    # Check what happens
    air_api = Airr("human", adaptable=False)
    airr_table = air_api.run_single("adaptable_seq", test_sequence)
    airr_entry = airr_table.iloc[0]
    cdr3_ = airr_entry["cdr3_aa"]
    cdr2_ = airr_entry["cdr2_aa"]
    cdr1_ = airr_entry["cdr1_aa"]
    fw1_ = airr_entry["fwr1_aa"]
    fw2_ = airr_entry["fwr2_aa"]
    fw3_ = airr_entry["fwr3_aa"]
    fw4_ = airr_entry["fwr4_aa"]
    assert fw1_ == "DIQMTQSPSSLSASVGDRVTITCQAS"
    assert fw2_ == "LNWYQQKPGKAPKLLIY"
    assert fw3_ == "NLETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYC"
    assert isnan(fw4_)
    assert cdr1_ == "QDISNY"
    assert cdr2_ == "DAS"
    assert isnan(cdr3_)

    air_api = Airr("human", adaptable=True)
    airr_table = air_api.run_single("adaptable_seq", test_sequence)
    airr_entry = airr_table.iloc[0]
    cdr3_ = airr_entry["cdr3_aa"]
    cdr2_ = airr_entry["cdr2_aa"]
    cdr1_ = airr_entry["cdr1_aa"]
    fw1_ = airr_entry["fwr1_aa"]
    fw2_ = airr_entry["fwr2_aa"]
    fw3_ = airr_entry["fwr3_aa"]
    fw4_ = airr_entry["fwr4_aa"]
    assert fw1_ == "DIQMTQSPSSLSASVGDRVTITCQAS"
    assert fw2_ == "LNWYQQKPGKAPKLLIY"
    assert fw3_ == "NLETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYC"
    assert fw4_ == "FGGGTKVDIK"
    assert cdr1_ == "QDISNY"
    assert cdr2_ == "DAS"
    assert cdr3_ == "QQYDN"

    scfv = """GGCGGCCGAGCTCGACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACCGTGTCACCATCACTT
              GCCGGGGCAAGTCAAATCATTAGCAACTATCTAAACTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTA
              TGCTACATCCAATTTGCACAGTGGGGTCCCATCACGCTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCA
              GCAGTCTGCAACCTGAAGATTTTGCAACTTACTACTGTCAACAGAGTTACAGTTTCCCGCCCACTTTCGGCGGAGGTACC
              AAGGTGGAGATCAAACGTACGGTTGCCGCTCCTTCTGTATTCATATTTCCGCCCTCCGATGAACAGCTTAAATCGGGCAC
              TGCTTCGGTAGTCTGCCTTCTGAATAATTTCTATCCCCGCGAGGCCAAGGTGCAATGGAAAGTCGACAATGCACTGCAAA
              GTGGAAACTCGCAAGAAAGCGTCACCGAACAGGACAGTAAAGATTCCACCTATAGCCTGTCATCGACACTTACCCTGAGT
              AAGGCTGATTACGAAAAGCACAAGGTTTACGCTTGCGAAGTAACTCACCAGGGCCTCTCAAGCCCTGTTACAAAGTCATT
              TAACAGAGGGGAATGCTAATTCTAGATAATTATCAAGGAGACAGTCATAATGAAATACCTATTGCCTACGGCAGCCGCTG
              GATTGTTATTACTCGCTGCCCAACCAGCCATGGCCGAGGTGCAGCTGCTCGAGGTGCAGCTGTTGGAGTCTGGGGGAGGC
              TTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGACTGACGATCTCATCCTACGCAATGTCATGGGT
              CCGCCAGGCTCCAGGGAANGGGGCTGGAGT"""
    air_api = Airr("human", adaptable=True)
    airr_table = air_api.run_single("scfv", scfv.replace("\n", "").replace(" ", ""), scfv=True)
    assert airr_table.table.fwr1_aa_heavy.iloc[0] == "EVQLLESGGGLVQPGGSLRLSCAAS"
    assert airr_table.table.fwr2_aa_heavy.iloc[0] == "MSWVRQAPGXGAGV"
    assert isinstance(airr_table, ScfvAirrTable)


# def test_mutational_analysis():
#     integration_file = "tests/integration/airr/catnap_nt_heavy.fasta.gz"
#     airr_api = Airr("human")
#     airrtable = airr_api.run_file(integration_file)


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
