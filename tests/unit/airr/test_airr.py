"""Unit tests for analysis interface."""
import os
import platform
import tempfile
from itertools import product
from pathlib import Path

import pandas as pd
import pytest
import semantic_version
from Bio.SeqRecord import SeqRecord
from click.testing import CliRunner
from numpy import isnan
from sadie.airr import Airr, AirrTable, BadDataSet, BadRequstedFileType, GermlineData, LinkedAirrTable
from sadie.airr import __file__ as sadie_airr_file
from sadie.airr import exceptions as airr_exceptions
from sadie.airr import igblast
from sadie.airr import methods as airr_methods
from sadie.app import airr as sadie_airr
from sadie.reference import Reference


def test_antibody_igblast_setup():
    """
    testing if we can manually setup IgBlast databases
    """

    system = platform.system().lower()
    executable = os.path.join(os.path.dirname(sadie_airr_file), f"bin/{system}/igblastn")

    ig_blast = igblast.IgBLASTN(executable)
    assert ig_blast.version == semantic_version.Version("1.17.1")
    germline_ref = os.path.join(os.path.dirname(os.path.abspath(sadie_airr_file)), "data/germlines")
    assert os.path.exists(germline_ref)

    # Set data
    for species in ["human", "mouse", "rat", "dog"]:
        db_ref = os.path.join(germline_ref, "imgt/Ig/blastdb/")
        internal_ref = os.path.join(germline_ref, "imgt/Ig/")
        aux_ref = os.path.join(germline_ref, "imgt/aux_db/")
        assert os.path.exists(db_ref)
        ig_blast.germline_db_v = os.path.join(db_ref, "{}_V".format(species))
        ig_blast.germline_db_d = os.path.join(db_ref, "{}_D".format(species))
        ig_blast.germline_db_j = os.path.join(db_ref, "{}_J".format(species))

        aux_ref = os.path.join(aux_ref, f"{species}_gl.aux")
        ig_blast.aux_path = aux_ref
        ig_blast.organism = species
        ig_blast.igdata = internal_ref
        ig_blast._pre_check()
    assert os.path.exists(germline_ref)

    # bad germeline ref
    germline_ref = "reference/germlines/"

    # Set data
    with pytest.raises(airr_exceptions.BadIgDATA):
        ig_blast.igdata = germline_ref
        ig_blast._pre_check()
    for species in ["human", "mouse", "rat", "dog"]:
        db_ref = os.path.join(germline_ref, "imgt/Ig/blastdb/")
        internal_ref = os.path.join(germline_ref, "imgt/Ig/")
        aux_ref = os.path.join(germline_ref, "imgt/aux_db/")
        with pytest.raises(airr_exceptions.BadIgBLASTArgument):
            ig_blast.germline_db_v = os.path.join(db_ref, "{}_V".format(species))
        with pytest.warns(UserWarning):
            ig_blast.germline_db_d = os.path.join(db_ref, "{}_D".format(species))
        with pytest.raises(airr_exceptions.BadIgBLASTArgument):
            ig_blast.germline_db_j = os.path.join(db_ref, "{}_J".format(species))
        with pytest.raises(airr_exceptions.BadIgBLASTArgument):
            ig_blast.aux_path = aux_ref

    with pytest.raises(airr_exceptions.BadIgBLASTExe):
        # pass an executable that is not igblast
        igblast.IgBLASTN("ls")


def test_antibody_igblast_file_run(fixture_setup):
    """
    testing if we can manually run igblast
    """

    # this should be fixture
    system = platform.system().lower()
    executable = os.path.join(os.path.dirname(sadie_airr_file), f"bin/{system}/igblastn")

    ig_blast = igblast.IgBLASTN(executable)
    germline_ref = os.path.join(os.path.dirname(os.path.abspath(sadie_airr_file)), "data/germlines")
    db_ref = os.path.join(germline_ref, "imgt/Ig/blastdb/")
    aux_path = os.path.join(germline_ref, "imgt/aux_db/")

    # Set data
    ig_blast.igdata = os.path.join(germline_ref, "imgt/Ig")

    # Grab from fixtures
    query = fixture_setup.get_pg9_heavy_fasta()
    ig_blast.germline_db_v = os.path.join(db_ref, "human_V")
    ig_blast.germline_db_d = os.path.join(db_ref, "human_D")
    ig_blast.germline_db_j = os.path.join(db_ref, "human_J")
    ig_blast.aux_path = os.path.join(aux_path, "human_gl.aux")
    ig_blast.organism = "human"
    ig_blast._pre_check()
    # have to make this -1 to get a more specifc j gene match
    ig_blast.j_penalty = -1
    csv_dataframe = ig_blast.run_file(query).fillna("")

    query = fixture_setup.get_pg9_heavy_fasta()
    csv_dataframe = ig_blast.run_file(query)
    csv_dataframe["d_call"] = csv_dataframe["d_call"].fillna("")


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
    """Test if we can init airr"""
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


def test_airrtable_init(fixture_setup):
    test_csv = fixture_setup.get_dog_airrtable()
    airr_table = AirrTable.read_csv(test_csv)
    assert not airr_table.empty
    non_airr_columns = [
        "d_call_top",
        "d_mutation",
        "d_mutation_aa",
        "j_call_top",
        "j_mutation",
        "j_mutation_aa",
        "liable",
        "note",
        "species",
        "v_call_top",
        "v_mutation",
        "v_mutation_aa",
        "vdj_aa",
        "vdj_nt",
    ]
    assert sorted(airr_table.non_airr_columns) == non_airr_columns
    # Tes twe can init with a direct pandas call
    airr_table = AirrTable(pd.read_csv(test_csv))
    assert not airr_table.empty
    assert isinstance(airr_table, AirrTable)

    # test we can read in json too
    test_json = fixture_setup.get_json_as_dataframe()
    airr_table = AirrTable.read_json(test_json)
    assert not airr_table.empty

    # gen bank
    genbanks = airr_table.get_genbank()
    assert all([isinstance(i, SeqRecord) for i in genbanks])

    # I will not accept a busted table sam I am
    busted_table = fixture_setup.get_busted_airrtable()
    with pytest.raises(airr_exceptions.MissingAirrColumns) as e:
        AirrTable.read_csv(busted_table)
    assert e.value.__str__()


def test_indel_correction(fixture_setup):
    test_csv = fixture_setup.get_dog_airrtable()
    # Test we can initllize with staic meathod
    airr_table = AirrTable.read_csv(test_csv)
    airr_table_indel = AirrTable.correct_indel(airr_table)
    assert "germline_alignment_aa_corrected" in airr_table_indel.columns
    assert "v_germline_alignment_aa_corrected" in airr_table_indel.columns
    assert isinstance(airr_table_indel, AirrTable)


def test_scfv_airrtable(fixture_setup):
    file_path = fixture_setup.get_linked_airrtable()
    dummy_scfv_table = pd.read_csv(file_path, index_col=0)
    linked_table = LinkedAirrTable(dummy_scfv_table)
    # test if we can split
    heavy_table, light_table = linked_table.get_split_table()
    assert isinstance(heavy_table, AirrTable)
    assert isinstance(light_table, AirrTable)
    rebuild_table = LinkedAirrTable(heavy_table.merge(light_table, on="sequence_id", suffixes=["_heavy", "_light"]))
    assert rebuild_table == rebuild_table
    assert rebuild_table == linked_table

    heavy_table["cell_id"] = heavy_table["sequence_id"]
    light_table["cell_id"] = light_table["sequence_id"]
    with pytest.raises(airr_exceptions.MissingAirrColumns):
        LinkedAirrTable(heavy_table.merge(light_table, on="cell_id", suffixes=["_heavy", "_light"]))
    rebuild_data = LinkedAirrTable(
        heavy_table.merge(light_table, on="cell_id", suffixes=["_heavy", "_light"]), key_column="cell_id"
    )
    heavy_table_split, light_table_split = rebuild_data.get_split_table()
    assert heavy_table_split.columns.difference(heavy_table.columns).empty
    assert light_table_split.columns.difference(light_table.columns).empty
    assert heavy_table == heavy_table_split[heavy_table.columns]
    assert light_table == light_table_split[light_table.columns]
    assert rebuild_data.key_column == "cell_id"
    assert rebuild_data.suffixes == ["_heavy", "_light"]

    rebuild_data = LinkedAirrTable(
        heavy_table.merge(light_table, on="cell_id", suffixes=["_h", "_l"]), suffixes=["_h", "_l"], key_column="cell_id"
    )
    assert rebuild_data.suffixes == ["_h", "_l"]


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
