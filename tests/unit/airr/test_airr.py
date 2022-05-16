"""Unit tests for analysis interface."""
import os
import platform
import shutil
from pathlib import Path

import pandas as pd
import pytest
import semantic_version
from Bio.SeqRecord import SeqRecord
from numpy import isnan
from sadie.airr import Airr, AirrSeries, AirrTable, GermlineData, LinkedAirrTable
from sadie.airr.exceptions import BadDataSet, BadRequstedFileType
from sadie.airr import __file__ as sadie_airr_file
from sadie.airr import exceptions as airr_exceptions
from sadie.airr import igblast
from sadie.airr import methods as airr_methods
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
        igblast.IgBLASTN("find")


def test_antibody_igblast_file_run(fixture_setup):
    """
    testing if we can manually run igblast withou airr abstraction
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


def test_airr_init(tmpdir, monkeypatch):
    """Test if we can init airr"""
    monkeypatch.setenv("TMPDIR", str(tmpdir / "monkeyairr"))
    for species in ["human", "mouse", "rat", "dog"]:
        air_api = Airr(species)
        air_api.get_available_datasets()
        air_api.get_available_species()
        assert isinstance(air_api, Airr)
    # show we can catch bad species inputs
    with pytest.raises(BadDataSet) as execinfo:
        air_api = Airr("robot")
    assert execinfo.value.__str__()

    # have an already created directory
    temporary_directory = tmpdir / "temp_output"
    air_api = Airr("human", temp_directory=str(temporary_directory))
    shutil.rmtree(temporary_directory)
    assert air_api.__repr__() == air_api.__str__()

    # odd settings produce coverage
    monkeypatch.delenv("TMPDIR")
    air_api = Airr("human", allow_vdj_overlap=True, d_gene_penalty=-1, j_gene_penalty=-1)
    with pytest.raises(TypeError):
        air_api.adapt_penalty = "robot"

    # test a temp file that is not a directory
    with pytest.raises(TypeError):
        air_api = Airr("human", temp_directory=30)

    # test a non writeable file
    with pytest.raises(IOError):
        air_api = Airr("human", temp_directory="/usr/bin")

    # test air init with bad igblastn executable.
    # this mocks that igblastn was not shipped in the repo and to search along path
    monkeypatch.setenv("IGBLASTN_MONEKY_PATCH", "blastp")
    with pytest.raises(airr_exceptions.BadIgBLASTExe):
        air_api = Airr("human")
    with pytest.raises(airr_exceptions.BadIgBLASTExe):
        air_api = Airr("human", igblast_exe="robot")

    monkeypatch.delenv("IGBLASTN_MONEKY_PATCH")
    air_api = Airr("human")

    # coverage for string and executable not being an executable file
    with pytest.raises(airr_exceptions.BadIgBLASTExe):
        air_api.executable = __file__

    # has to be IgBlastN class
    with pytest.raises(TypeError):
        air_api.igblast = str


def test_custom_mice_init():
    """Test we can initialize custom Mice"""
    Airr("se09")


def test_airr_single_sequence(fixture_setup):
    """Test we can run a single sequence."""
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
    assert round(v_mutation, 2) == 13.99
    assert round(d_mutation, 2) == 17.86
    assert round(j_mutation, 2) == 11.32
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
    """Test we can run with seq records biopython instances"""
    airr = Airr("human", adaptable=False)

    airr.run_records([fixture_setup.get_pg9_heavy_sequence(), fixture_setup.get_pg9_light_sequence()])
    with pytest.raises(TypeError):
        airr.run_records(fixture_setup.get_pg9_heavy_sequence())


# def test_scfv(fixture_setup):
# """Run on short scfv sequence"""
# airr = Airr("human", adaptable=False)
# linked = airr.run_records(fixture_setup.get_scfv_sequences(), scfv=True)
# assert isinstance(linked, LinkedAirrTable)

# """Run on multiple long scfv sequences"""
# scfv = fixture_setup.get_long_scfv_fastq()
# airr = Airr("human", adaptable=False)
# results = airr.run_records(scfv, scfv=True)
# assert isinstance(results, LinkedAirrTable)


def test_airr_from_dataframe(fixture_setup):
    """Test we can pass a dataframe to runtime"""
    dog_df = pd.read_csv(fixture_setup.get_dog_airrtable(), sep="\t")
    airr_api = Airr("dog")
    unjoined_df = airr_api.run_dataframe(dog_df, "sequence_id", "sequence")
    assert isinstance(unjoined_df, AirrTable)
    joined_df = airr_api.run_dataframe(dog_df, "sequence_id", "sequence", return_join=True)
    assert isinstance(joined_df, pd.DataFrame)

    with pytest.raises(TypeError):
        airr_api.run_dataframe(fixture_setup.get_oas_fasta(), "sequence_id", "sequence")

    double_dog = pd.concat([dog_df, dog_df])
    with pytest.raises(ValueError):
        airr_api.run_dataframe(double_dog, "sequence_id", "sequence")

    with pytest.raises(ValueError):
        airr_api.run_dataframe(fixture_setup.get_dog_airrtable_with_missing_sequences(), "sequence_id", "sequence")


def test_airr_from_file(fixture_setup):
    """Test we can pass a fasta file to runtime"""
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

    # more coverage
    airr_api = Airr("human")
    with pytest.raises(FileNotFoundError):
        airr_api.run_fasta("a_non_existant_directory/a_non_existant_file.fasta")


def test_vdj_overlap(fixture_setup):
    """Test if we can allow vdj_overlap"""
    air_api = Airr("human", allow_vdj_overlap=True, adaptable=False)
    air_api.run_single("PG9", fixture_setup.get_pg9_heavy_sequence().seq.__str__())
    assert air_api.igblast.allow_vdj_overlap.value is True

    #  make sure a warning is thrown if we have both adapt penalty and vdj overlap
    with pytest.warns(UserWarning):
        air_api = Airr("human", allow_vdj_overlap=True, adaptable=True)


def test_adaptable_penalty(fixture_setup):
    """Test we can run adapatable pentalty methods"""
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


def test_adaptable_correction(fixture_setup):
    air_api = Airr("human", adaptable=False)

    # these sequences are liable.
    liable = air_api.run_fasta(fixture_setup.get_OAS_correctable_pentalty_file())
    assert (liable["v_penalty"] == -1).all()
    assert (liable["d_penalty"] == -1).all()
    assert (liable["j_penalty"] == -2).all()

    air_api = Airr("human", adaptable=True)
    liable_corrected = air_api.run_dataframe(liable, "sequence_id", "sequence")
    # correct 30 of them
    assert not (liable_corrected["liable"]).any()

    liable_mix = air_api.run_fasta(fixture_setup.get_OAS_liable_file())
    assert liable_mix["liable"].value_counts().to_dict() == {True: 30, False: 24}


def test_mutational_analysis(heavy_catnap_airrtable, light_catnap_airrtable):
    """Test we can run mutational analysis post hoc method"""

    # run on heavy from fixture of airr
    heavy_airrtable_with_analysis = airr_methods.run_mutational_analysis(
        heavy_catnap_airrtable, "kabat", run_multiproc=True
    )
    assert "mutations" in heavy_airrtable_with_analysis.columns

    # run on light
    light_airrtable_with_analysis = airr_methods.run_mutational_analysis(
        light_catnap_airrtable, "kabat", run_multiproc=False
    )
    assert "mutations" in light_airrtable_with_analysis.columns
    link_table = heavy_catnap_airrtable.merge(light_catnap_airrtable, on="cellid", suffixes=["_heavy", "_light"])

    # make sure we can run it on linked airr table
    joined_airr_table = LinkedAirrTable(link_table, key_column="cellid")
    joined_airr_table_with_analysis = airr_methods.run_mutational_analysis(
        joined_airr_table, "kabat", run_multiproc=False
    )
    assert "mutations_heavy" in joined_airr_table_with_analysis.columns
    assert "mutations_light" in joined_airr_table_with_analysis.columns
    return heavy_airrtable_with_analysis, light_airrtable_with_analysis, joined_airr_table_with_analysis


def test_igl_assignment(heavy_catnap_airrtable, light_catnap_airrtable):
    """test if we can assign igl through method analysis"""
    airrtable_with_igl = airr_methods.run_igl_assignment(heavy_catnap_airrtable)
    assert "iGL" in airrtable_with_igl.columns

    airrtable_with_igl = airr_methods.run_igl_assignment(light_catnap_airrtable)
    assert "iGL" in airrtable_with_igl.columns

    link_table = heavy_catnap_airrtable.merge(light_catnap_airrtable, on="cellid", suffixes=["_heavy", "_light"])
    joined_airr_table = LinkedAirrTable(link_table, key_column="cellid")
    joined_airr_table_with_analysis = airr_methods.run_igl_assignment(joined_airr_table)
    assert "iGL_heavy" in joined_airr_table_with_analysis.columns
    assert "iGL_light" in joined_airr_table_with_analysis.columns


def test_runtime_reference(fixture_setup):
    """Test JIT reference generation"""
    reference = Reference()
    reference.add_gene({"species": "custom", "sub_species": "human", "gene": "IGHV1-2*01", "database": "imgt"})
    reference.add_gene({"species": "custom", "sub_species": "macaque", "gene": "IGHV1-105*01", "database": "custom"})

    # make sure we add a J gene
    with pytest.raises(ValueError):
        airr_api = Airr(reference, adaptable=True)
    reference.add_gene({"species": "custom", "sub_species": "human", "gene": "IGHJ6*01", "database": "imgt"})

    # make sure that we don't accept heterogenous database types
    with pytest.raises(ValueError):
        airr_api = Airr(reference, adaptable=True)

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
    test_tsv = fixture_setup.get_dog_airrtable()
    airr_table = AirrTable.read_airr(test_tsv)
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
    airr_table = AirrTable(pd.read_csv(test_tsv, sep="\t"))
    assert not airr_table.empty
    assert isinstance(airr_table, AirrTable)

    # test we can read in json too
    test_json = fixture_setup.get_json_as_dataframe()
    airr_table = AirrTable(pd.read_json(test_json, orient="records"))
    assert not airr_table.empty

    # gen bank
    genbanks = airr_table.get_genbank()
    assert all([isinstance(i, SeqRecord) for i in genbanks])

    # I will not accept a busted table sam I am
    busted_table = fixture_setup.get_busted_airrtable()
    with pytest.raises(airr_exceptions.MissingAirrColumns) as e:
        AirrTable.read_airr(busted_table)
    assert e.value.__str__()


def test_indel_correction(fixture_setup):
    test_csv = fixture_setup.get_dog_airrtable()
    # Test we can initllize with staic meathod
    airr_table = AirrTable.read_airr(test_csv)
    airr_table_indel = AirrTable.correct_indel(airr_table)
    assert "germline_alignment_aa_corrected" in airr_table_indel.columns
    assert "v_germline_alignment_aa_corrected" in airr_table_indel.columns
    assert isinstance(airr_table_indel, AirrTable)


# def test_scfv_airrtable(fixture_setup):
#     file_path = fixture_setup.get_linked_airrtable()
#     dummy_scfv_table = pd.read_csv(file_path, index_col=0)
#     linked_table = LinkedAirrTable(dummy_scfv_table)
#     # test if we can split
#     heavy_table, light_table = linked_table.get_split_table()
#     assert isinstance(heavy_table, AirrTable)
#     assert isinstance(light_table, AirrTable)
#     rebuild_table = LinkedAirrTable(heavy_table.merge(light_table, on="sequence_id", suffixes=["_heavy", "_light"]))
#     assert rebuild_table == rebuild_table
#     assert rebuild_table == linked_table

#     heavy_table["cell_id"] = heavy_table["sequence_id"]
#     light_table["cell_id"] = light_table["sequence_id"]
#     with pytest.raises(airr_exceptions.MissingAirrColumns):
#         LinkedAirrTable(heavy_table.merge(light_table, on="cell_id", suffixes=["_heavy", "_light"]))
#     rebuild_data = LinkedAirrTable(
#         heavy_table.merge(light_table, on="cell_id", suffixes=["_heavy", "_light"]), key_column="cell_id"
#     )
#     heavy_table_split, light_table_split = rebuild_data.get_split_table()
#     assert heavy_table_split.columns.difference(heavy_table.columns).empty
#     assert light_table_split.columns.difference(light_table.columns).empty
#     assert heavy_table == heavy_table_split[heavy_table.columns]
#     assert light_table == light_table_split[light_table.columns]
#     assert rebuild_data.key_column == "cell_id"
#     assert rebuild_data.suffixes == ["_heavy", "_light"]

#     rebuild_data = LinkedAirrTable(
#         heavy_table.merge(light_table, on="cell_id", suffixes=["_h", "_l"]), suffixes=["_h", "_l"], key_column="cell_id"
#     )
#     assert rebuild_data.suffixes == ["_h", "_l"]


def test_write_and_check_airr(fixture_setup):
    """Check that the offical airr can validate our airr tables"""
    catnap_heavy = fixture_setup.get_catnap_heavy_nt()
    output_file = fixture_setup.tmp_path / Path("test_write_and_check_airr.tsv")
    airr_api = Airr("human")
    airr_table = airr_api.run_fasta(catnap_heavy)
    assert isinstance(airr_table, AirrTable)
    airr_table.to_airr(output_file)
    import airr

    # the official airr package will validate
    d = airr.load_rearrangement(output_file, debug=True, validate=True)
    assert isinstance(d, pd.DataFrame)

    # if we drop a requried field, the official airr module should cause a validation error
    d.drop("sequence_id", axis=1).to_csv(output_file, sep="\t")
    with pytest.raises(airr.schema.ValidationError):
        airr.load_rearrangement(output_file, debug=True, validate=True)


def test_airr_series(fixture_setup):
    """Test that we can create an airr series"""
    df = pd.read_feather(fixture_setup.get_catnap_heavy_airrtable())
    series = AirrSeries(df.iloc[0])
    assert isinstance(series, AirrSeries)
    table = AirrTable([series, series])
    assert isinstance(table, AirrTable)
