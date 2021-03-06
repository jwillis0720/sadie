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
from sadie.airr import Airr, AirrTable, BadDataSet, BadRequstedFileType, GermlineData, LinkedAirrTable
from sadie.airr import methods as airr_methods
from sadie.app import airr as sadie_airr

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


def test_custom_mice_init():
    Airr("hugl18")
    Airr("se09")
    Airr("se16")
    Airr("se684")


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
    air_api = Airr("human", adaptable=False)
    airr_table = air_api.run_single("PG9", pg9_seq)
    airr_entry = airr_table.iloc[0]
    seq_id = airr_entry["sequence_id"]
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
    assert seq_id == "PG9"

    # will def change based on penalties, so be careful
    assert round(v_mutation, 2) == 14.0
    assert round(d_mutation, 2) == 17.88
    assert round(j_mutation, 2) == 11.31
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
    airr = Airr("human", adaptable=False)
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
    airr.run_records([seq_record, seq_record])
    with pytest.raises(TypeError):
        airr.run_records(seq_record)


def test_run_multiple_scfv():
    scfv = get_file("fastq_inputs/long_scfv.fastq.gz")
    airr = Airr("human", adaptable=False)
    list_to_run = list(SeqIO.parse(gzip.open(scfv, "rt"), "fastq"))
    results = airr.run_records(list_to_run, scfv=True)
    assert isinstance(results, LinkedAirrTable)


def test_vdj_overlap():
    pg9_seq = """
        CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGT
        CCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGT
        GGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCC
        AGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATT
        TTTGTGTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTA
        TTATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC""".replace(
        "\n", ""
    )
    air_api = Airr("human", allow_vdj_overlap=True, adaptable=False)
    air_api.run_single("PG9", pg9_seq)
    assert air_api.igblast.allow_vdj_overlap.value is True


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
    f = get_file("fasta_inputs/PG9_H_multiple.fasta")
    airr_api = Airr("human")
    result = airr_api.run_fasta(f)
    assert isinstance(result, AirrTable)

    # see if we can take in paths
    result = airr_api.run_fasta(Path(f))
    assert isinstance(result, AirrTable)

    f = get_file("card.png")
    with pytest.raises(BadRequstedFileType) as execinfo:
        airr_api.run_fasta(f)
    assert execinfo.value.__str__()

    f = get_file("fasta_inputs/scfv.fasta")
    airr_api = Airr("human", temp_directory="a_non_existant_directory")
    airr_api.run_fasta(f, scfv=True)

    # cleanup just in case test fails
    if Path("a_non_existant_directory").exists():
        Path("a_non_existant_directory").rmdir()


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
    assert airr_table.fwr1_aa_heavy.iloc[0] == "EVQLLESGGGLVQPGGSLRLSCAAS"
    assert airr_table.fwr2_aa_heavy.iloc[0] == "MSWVRQAPGXGAGV"
    assert isinstance(airr_table, pd.DataFrame)
    assert isinstance(airr_table, AirrTable)
    assert isinstance(airr_table, LinkedAirrTable)

    integration_file = "tests/integration/airr/fixtures/catnap_nt_heavy.fasta"
    air_api = Airr("human", adaptable=True)
    liable = air_api.run_fasta(integration_file)
    assert not liable["liable"].any()
    integration_file = "tests/integration/airr/fixtures/catnap_nt_light.fasta"
    liable = air_api.run_fasta(integration_file)
    assert not liable["liable"].any()


def test_mutational_analysis():
    integration_file = "tests/integration/airr/fixtures/catnap_nt_heavy.fasta"
    airr_api = Airr("human")
    airrtable_heavy = airr_api.run_fasta(integration_file)
    airrtable_heavy["cellid"] = airrtable_heavy["sequence_id"].str.split("_").str.get(0)
    airrtable_with_analysis = airr_methods.run_mutational_analysis(airrtable_heavy, "kabat")
    assert "mutations" in airrtable_with_analysis.columns

    integration_file = "tests/integration/airr/fixtures/catnap_nt_light.fasta"
    airr_api = Airr("human")
    airrtable_light = airr_api.run_fasta(integration_file)
    airrtable_light["cellid"] = airrtable_light["sequence_id"].str.split("_").str.get(0)
    airrtable_with_analysis = airr_methods.run_mutational_analysis(airrtable_light, "kabat")
    assert "mutations" in airrtable_with_analysis.columns
    link_table = airrtable_heavy.merge(airrtable_light, on="cellid", suffixes=["_heavy", "_light"])

    joined_airr_table = LinkedAirrTable(link_table, key_column="cellid")
    joined_airr_table_with_analysis = airr_methods.run_mutational_analysis(joined_airr_table, "kabat")
    assert "mutations_heavy" in joined_airr_table_with_analysis.columns
    assert "mutations_light" in joined_airr_table_with_analysis.columns


def test_igl_assignment():
    integration_file = "tests/integration/airr/fixtures/catnap_nt_heavy.fasta"
    airr_api = Airr("human")
    airrtable_heavy = airr_api.run_fasta(integration_file)
    airrtable_heavy["cellid"] = airrtable_heavy["sequence_id"].str.split("_").str.get(0)
    airrtable_with_igl = airr_methods.run_igl_assignment(airrtable_heavy)
    assert "iGL" in airrtable_with_igl.columns

    integration_file = "tests/integration/airr/fixtures/catnap_nt_light.fasta"
    airr_api = Airr("human")
    airrtable_light = airr_api.run_fasta(integration_file)
    airrtable_light["cellid"] = airrtable_light["sequence_id"].str.split("_").str.get(0)
    airrtable_with_igl = airr_methods.run_igl_assignment(airrtable_light)
    assert "iGL" in airrtable_with_igl.columns

    link_table = airrtable_heavy.merge(airrtable_light, on="cellid", suffixes=["_heavy", "_light"])
    joined_airr_table = LinkedAirrTable(link_table, key_column="cellid")
    joined_airr_table_with_analysis = airr_methods.run_igl_assignment(joined_airr_table)
    assert "iGL_heavy" in joined_airr_table_with_analysis.columns
    assert "iGL_light" in joined_airr_table_with_analysis.columns


def _run_cli(args, tmpfile):
    runner = CliRunner(echo_stdin=True)
    result = runner.invoke(sadie_airr, args)
    if result.exit_code != 0:
        logger.error(f"Exception - {result.exception.__str__()}")
        logger.error(f"{args} produces error")
    assert result.exit_code == 0
    assert os.path.exists(tmpfile.name)
    return True


def test_cli(caplog):
    """Confirm the CLI works as expecte"""

    # we need this fixture because... well i don't know
    # caplog.set_level(20)
    # IMGT DB
    quereies = glob.glob(get_file("fasta_inputs/") + "*")
    species = ["dog", "rat", "human", "mouse", "macaque", "se09"]
    # ft = ["csv", "json"]
    products = product(species, ["imgt"], quereies)

    with tempfile.NamedTemporaryFile(suffix=".csv") as tmpfile:
        for p_tuple in products:
            cli_input = [
                "-v",
                "--species",
                p_tuple[0],
                "--db-type",
                p_tuple[1],
                p_tuple[2],
                tmpfile.name,
                "--skip-mutation",
            ]
            print(f"CLI input {' '.join(cli_input)}")
            test_success = _run_cli(cli_input, tmpfile)
            assert test_success


def test_cli_custom():
    quereies = glob.glob(get_file("fasta_inputs/") + "*")
    species = ["cat", "macaque", "dog"]
    products = product(species, ["custom"], quereies)

    with tempfile.NamedTemporaryFile(suffix=".csv") as tmpfile:
        for p_tuple in products:
            cli_input = [
                "-v",
                "--species",
                p_tuple[0],
                "--db-type",
                p_tuple[1],
                p_tuple[2],
                tmpfile.name,
                "--skip-mutation",
            ]
            print(f"CLI input {' '.join(cli_input)}")
            test_success = _run_cli(cli_input, tmpfile)
            assert test_success


def test_cli_scfv():
    quereies = [get_file("fasta_inputs/scfv.fasta.gz")]
    species = ["human"]
    products = product(species, ["imgt"], quereies)
    with tempfile.NamedTemporaryFile(suffix=".csv") as tmpfile:
        for p_tuple in products:
            cli_input = [
                "-v",
                "--species",
                p_tuple[0],
                "--db-type",
                p_tuple[1],
                p_tuple[2],
                tmpfile.name,
                "--skip-mutation",
            ]
            print(f"CLI input {' '.join(cli_input)}")
            test_success = _run_cli(cli_input, tmpfile)
            assert test_success


def test_edge_cases():
    airr_api = Airr("macaque", database="custom", adaptable=False)
    airr_api.run_single(
        "bad_seq",
        "TCCAGTCCCTGCAGGCCGGGAGGCAGGTGACCTCTGCCTCAGACCCCCACTCCAGACACCAGACAGAGGGGCAGGCCCCCCAGAACCAAAGTGGAGGGACGACCCGTCAAGGACAAACCAGACCAAGGGACACTGAGCCCAGCACGGGAAGGTCCCCAGATAGACCAGGAGGTTTCTGGAGGTGTCTGTGCCACAGTGGGGTATAGCAGCAGATCCGACTACGGTAGCAACTTTTGGGACTACTGGGGCCAGGGAGTCCTGGTCACCGTCTCCTCAGCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTCCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGC",
    )
