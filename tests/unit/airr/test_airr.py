"""Unit tests for analysis interface."""
import os
import tempfile

from itertools import product

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner
from numpy import isnan
from sadie.airr import Airr, AirrTable, BadDataSet, BadRequstedFileType, GermlineData, LinkedAirrTable
from sadie.airr import methods as airr_methods
from sadie.app import airr as sadie_airr
from sadie.reference import Reference


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


def test_airr_single_sequence(fixture_setup):
    """Test we can run a single sequence"""
    pg9_heavy_seq = fixture_setup.get_pg9_heavy_sequence().seq.__str__()
    pg9_light_seq = fixture_setup.get_pg9_light_sequence().seq.__str__()
    air_api = Airr("human", adaptable=False)
    airr_table = air_api.run_single("PG9", pg9_heavy_seq)
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
        airr_table = air_api.run_single(9, pg9_heavy_seq)

    # super edge case, two heavy or two light scfv
    air_api.run_single("two_heavy", pg9_heavy_seq + pg9_heavy_seq, scfv=True)
    air_api.run_single("two_light", pg9_light_seq + pg9_light_seq, scfv=True)

    # run with 1 cpu
    air_api = Airr("human", adaptable=False, num_cpus=1)
    airr_table = air_api.run_single("PG9", pg9_heavy_seq)


def test_run_multiple(fixture_setup):
    # Run with seqrecords
    airr = Airr("human", adaptable=False)

    airr.run_records([fixture_setup.get_pg9_heavy_sequence(), fixture_setup.get_pg9_light_sequence()])
    with pytest.raises(TypeError):
        airr.run_records(fixture_setup.get_pg9_heavy_sequence())


def test_scfv(fixture_setup):
    airr = Airr("human", adaptable=False)
    linked = airr.run_records(fixture_setup.get_scfv_sequences(), scfv=True)
    assert isinstance(linked, LinkedAirrTable)


def test_run_multiple_scfv(fixture_setup):
    """Run on multiple long scfv sequences"""
    scfv = fixture_setup.get_long_scfv_fastq()
    airr = Airr("human", adaptable=False)
    results = airr.run_records(scfv, scfv=True)
    assert isinstance(results, LinkedAirrTable)


def test_vdj_overlap(fixture_setup):
    air_api = Airr("human", allow_vdj_overlap=True, adaptable=False)
    air_api.run_single("PG9", fixture_setup.get_pg9_heavy_sequence().seq.__str__())
    assert air_api.igblast.allow_vdj_overlap.value is True


def test_airr_from_dataframe(fixture_setup):
    """Test we can pass a dataframe to runtime"""
    dog_df = pd.read_csv(fixture_setup.get_dog_airrtable(), index_col=0)
    airr_api = Airr("dog")
    unjoined_df = airr_api.run_dataframe(dog_df, "sequence_id", "sequence")
    assert isinstance(unjoined_df, AirrTable)
    joined_df = airr_api.run_dataframe(dog_df, "sequence_id", "sequence", return_join=True)
    assert isinstance(joined_df, pd.DataFrame)


def test_airr_from_file(fixture_setup):
    """Test we can pass a dataframe to runtime"""
    f = str(fixture_setup.get_pg9_heavy_multiple_fasta())
    airr_api = Airr("human")
    result = airr_api.run_fasta(f)
    assert isinstance(result, AirrTable)

    # see if we can take in paths
    result = airr_api.run_fasta(f)
    assert isinstance(result, AirrTable)

    f = fixture_setup.get_card()
    with pytest.raises(BadRequstedFileType) as execinfo:
        airr_api.run_fasta(f)
    assert execinfo.value.__str__()

    f = fixture_setup.get_scfv_fasta()
    airr_api = Airr("human", temp_directory="a_non_existant_directory")
    airr_api.run_fasta(f, scfv=True)

    # cleanup just in case test fails
    if Path("a_non_existant_directory").exists():
        Path("a_non_existant_directory").rmdir()


def test_adaptable_penalty(fixture_setup):
    test_sequence = fixture_setup.get_adaptable_pentalty_test_seq()

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

    air_api = Airr("human", adaptable=True)
    airr_table = air_api.run_single("scfv", fixture_setup.get_adaptable_pentalty_test_seq_scfv(), scfv=True)
    assert airr_table.fwr1_aa_heavy.iloc[0] == "EVQLLESGGGLVQPGGSLRLSCAAS"
    assert airr_table.fwr2_aa_heavy.iloc[0] == "MSWVRQAPGXGAGV"
    assert isinstance(airr_table, pd.DataFrame)
    assert isinstance(airr_table, AirrTable)
    assert isinstance(airr_table, LinkedAirrTable)


def test_adaptable_catnap(fixture_setup):
    air_api = Airr("human", adaptable=True)
    liable = air_api.run_fasta(fixture_setup.get_catnap_heavy_nt())
    assert not liable["liable"].any()
    liable = air_api.run_fasta(fixture_setup.get_catnap_light_nt())
    assert not liable["liable"].any()


def test_mutational_analysis(heavy_catnap_airrtable, light_catnap_airrtable):
    # run on heavy from fixture of airr
    airrtable_with_analysis = airr_methods.run_mutational_analysis(heavy_catnap_airrtable, "kabat", run_multiproc=False)
    assert "mutations" in airrtable_with_analysis.columns

    # run on light
    airrtable_with_analysis = airr_methods.run_mutational_analysis(light_catnap_airrtable, "kabat")
    assert "mutations" in airrtable_with_analysis.columns
    link_table = heavy_catnap_airrtable.merge(light_catnap_airrtable, on="cellid", suffixes=["_heavy", "_light"])

    # make sure we can run it on linked airr table
    joined_airr_table = LinkedAirrTable(link_table, key_column="cellid")
    joined_airr_table_with_analysis = airr_methods.run_mutational_analysis(joined_airr_table, "kabat")
    assert "mutations_heavy" in joined_airr_table_with_analysis.columns
    assert "mutations_light" in joined_airr_table_with_analysis.columns


def test_igl_assignment(heavy_catnap_airrtable, light_catnap_airrtable):
    # test if we can assign igl
    airrtable_with_igl = airr_methods.run_igl_assignment(heavy_catnap_airrtable)
    assert "iGL" in airrtable_with_igl.columns

    airrtable_with_igl = airr_methods.run_igl_assignment(light_catnap_airrtable)
    assert "iGL" in airrtable_with_igl.columns

    link_table = heavy_catnap_airrtable.merge(light_catnap_airrtable, on="cellid", suffixes=["_heavy", "_light"])
    joined_airr_table = LinkedAirrTable(link_table, key_column="cellid")
    joined_airr_table_with_analysis = airr_methods.run_igl_assignment(joined_airr_table)
    assert "iGL_heavy" in joined_airr_table_with_analysis.columns
    assert "iGL_light" in joined_airr_table_with_analysis.columns


def test_runtime_referecne(fixture_setup):
    reference = Reference()
    reference.add_gene({"species": "custom", "sub_species": "human", "gene": "IGHV1-2*01", "database": "imgt"})
    reference.add_gene({"species": "custom", "sub_species": "human", "gene": "IGHV3-15*01", "database": "imgt"})
    reference.add_gene({"species": "custom", "sub_species": "human", "gene": "IGHJ6*01", "database": "imgt"})
    reference.add_gene({"species": "custom", "sub_species": "human", "gene": "IGHD3-3*01", "database": "imgt"})
    airr_api = Airr(reference, adaptable=True)
    airr_table = airr_api.run_single("PG9", fixture_setup.get_pg9_heavy_sequence().seq.__str__())
    assert airr_table.species.iloc[0] == "custom"
    assert airr_table.v_call_top.iloc[0] == "IGHV3-15*01"
    assert airr_table.species.iloc[0] == "custom"

    reference = Reference()
    reference.add_gene({"species": "custom", "sub_species": "human", "gene": "IGHV1-2*01", "database": "imgt"})
    reference.add_gene({"species": "custom", "sub_species": "dog", "gene": "IGHV1-15*01", "database": "imgt"})
    reference.add_gene({"species": "custom", "sub_species": "human", "gene": "IGHJ6*01", "database": "imgt"})
    reference.add_gene({"species": "custom", "sub_species": "human", "gene": "IGHD3-3*01", "database": "imgt"})
    airr_api = Airr(reference, adaptable=True)
    airr_table = airr_api.run_single("PG9", fixture_setup.get_pg9_heavy_sequence().seq.__str__())
    assert airr_table.species.iloc[0] == "custom"
    assert airr_table.v_call_top.iloc[0] == "human|IGHV1-2*01"


def _run_cli(args, tmpfile):
    runner = CliRunner(echo_stdin=True)
    result = runner.invoke(sadie_airr, args)
    if result.exit_code != 0:
        print(f"Exception - {result.exception.__str__()}")
        print(f"{args} produces error")
    assert result.exit_code == 0
    assert os.path.exists(tmpfile.name)
    return True


def test_cli(caplog, fixture_setup):
    """Confirm the CLI works as expecte"""

    # we need this fixture because... well i don't know
    # https://github.com/pallets/click/issues/824
    caplog.set_level(200000)
    # IMGT DB
    quereies = fixture_setup.get_fasta_files()
    species = ["dog", "rat", "human", "mouse", "macaque", "se09"]
    products = product(species, ["imgt"], quereies)

    with tempfile.NamedTemporaryFile(suffix=".csv") as tmpfile:
        # run 1 with mutational analysis
        cli_input = [
            "-v",
            "--species",
            "human",
            "--db-type",
            "imgt",
            str(fixture_setup.get_pg9_heavy_fasta()),
            tmpfile.name,
        ]
        print(f"CLI input {' '.join(cli_input)}")
        test_success = _run_cli(cli_input, tmpfile)
        assert test_success

        for p_tuple in products:
            cli_input = [
                "-v",
                "--species",
                p_tuple[0],
                "--db-type",
                p_tuple[1],
                str(p_tuple[2]),
                tmpfile.name,
                "--skip-mutation",
            ]
            print(f"CLI input {' '.join(cli_input)}")
            test_success = _run_cli(cli_input, tmpfile)
            assert test_success


def test_cli_custom(caplog, fixture_setup):
    queries = fixture_setup.get_fasta_files()
    species = ["cat", "macaque", "dog"]
    products = product(species, ["custom"], queries)
    caplog.set_level(200000)

    with tempfile.NamedTemporaryFile(suffix=".csv") as tmpfile:
        for p_tuple in products:
            cli_input = [
                "-v",
                "--species",
                p_tuple[0],
                "--db-type",
                p_tuple[1],
                str(p_tuple[2]),
                tmpfile.name,
                "--skip-mutation",
            ]
            print(f"CLI input {' '.join(cli_input)}")
            test_success = _run_cli(cli_input, tmpfile)
            assert test_success


def test_cli_scfv(caplog, fixture_setup):
    queries = str(fixture_setup.get_scfv_fasta())
    species = ["human"]
    caplog.set_level(200000)
    products = product(species, ["imgt"], [queries])
    with tempfile.NamedTemporaryFile(suffix=".csv") as tmpfile:
        for p_tuple in products:
            cli_input = [
                "-v",
                "--species",
                p_tuple[0],
                "--db-type",
                p_tuple[1],
                str(p_tuple[2]),
                tmpfile.name,
                "--skip-mutation",
            ]
            print(f"CLI input {' '.join(cli_input)}")
            test_success = _run_cli(cli_input, tmpfile)
            assert test_success


def test_edge_cases(fixture_setup):
    airr_api = Airr("macaque", database="custom", adaptable=False)
    airr_api.run_single("bad_seq", fixture_setup.get_monkey_edge_seq())
