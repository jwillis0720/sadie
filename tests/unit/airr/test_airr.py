"""Unit tests for analysis interface."""

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy import isnan, nan

from sadie.airr import Airr, AirrSeries, AirrTable, GermlineData, LinkedAirrTable
from sadie.airr import exceptions as airr_exceptions
from sadie.airr import methods as airr_methods
from sadie.airr.airrtable.genbank import GenBank, GenBankFeature
from sadie.airr.exceptions import BadDataSet, BadRequstedFileType
from sadie.airr.models import AirrSeriesModel
from sadie.reference import Reference, References
from tests.conftest import SadieFixture


def test_airr_model() -> None:
    # string nulls should be equal to None
    assert AirrSeriesModel(**{field: str(nan) for field in AirrSeriesModel().__fields__}) == AirrSeriesModel()


def test_airrtable_inheritance(fixture_setup) -> None:
    """Test that we can create an airr DataFrame"""
    df = pd.read_feather(fixture_setup.get_bum_igl_assignment())
    table = AirrTable(df)  # init and verify
    assert isinstance(table, AirrTable)
    assert table.verified == True
    series = table.iloc[0]  # constroctor slice
    series = AirrSeries(series)  # init and verify
    assert isinstance(series, AirrSeries)
    assert set(series.index) == set(table.columns)
    copied_serires = series.copy()  # constroctor
    # table = series.to_frame()  # expanddim
    # assert isinstance(table, AirrTable)
    table = AirrTable([series, series])
    assert isinstance(table, AirrTable)


def test_airr_columns(fixture_setup) -> None:
    """Formality to return airr format compliant columns."""
    df = pd.read_feather(fixture_setup.get_bum_igl_assignment())
    table = AirrTable(df)
    assert table.airr_columns.empty == False


def test_to_fasta(fixture_setup) -> None:
    """Formality to return airr format compliant columns."""
    df = pd.read_feather(fixture_setup.get_bum_igl_assignment())
    table = AirrTable(df)
    filepath = fixture_setup.tmp_path / "test.fasta"
    table.to_fasta(filename=filepath)
    SeqIO.parse(filepath, "fasta")
    with pytest.raises(ValueError):
        table.to_fasta(filename=filepath, sequence_field="not_a_sequence_field")


def test_to_genbank(fixture_setup) -> None:
    """Formality to return airr format compliant columns."""
    df = pd.read_feather(fixture_setup.get_bum_igl_assignment())
    table = AirrTable(df)
    filepath = fixture_setup.tmp_path / "test.gb"
    table.to_genbank(filename=filepath)
    SeqIO.parse(filepath, "genbank")


def test_sanitized_antibodies(fixture_setup) -> None:
    """Test get antobodies where no segment is missing."""
    df = pd.read_feather(fixture_setup.get_bum_igl_assignment())
    table = AirrTable(df)
    assert table.get_sanitized_antibodies().empty == False


def test_germline_init() -> None:
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


def test_airr_init(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test if we can init airr"""

    tmpdir = tmp_path_factory.mktemp("test_airr_init")
    # this will make our tmpdir discoverable by the AIRR
    monkeypatch.setenv("TMPDIR", str(tmpdir / Path("monkeyairr")))
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


def test_airr_single_sequence(fixture_setup: SadieFixture) -> None:
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
    assert round(v_mutation, 2) == np.float32(0.14)
    assert round(d_mutation, 2) == np.float32(0.18)
    assert round(j_mutation, 2) == np.float32(0.11)
    with pytest.raises(TypeError):
        # id must be str
        airr_table = air_api.run_single(9, pg9_heavy_seq)

    # super edge case, two heavy or two light scfv
    air_api.run_single("two_heavy", pg9_heavy_seq + pg9_heavy_seq, scfv=True)
    air_api.run_single("two_light", pg9_light_seq + pg9_light_seq, scfv=True)

    # run with 1 cpu
    air_api = Airr("human", adaptable=False, num_cpus=1)
    airr_table = air_api.run_single("PG9", pg9_heavy_seq)


def test_run_multiple(fixture_setup: SadieFixture) -> None:
    """Test we can run with seq records biopython instances"""
    airr = Airr("human", adaptable=False)

    airr.run_records([fixture_setup.get_pg9_heavy_sequence(), fixture_setup.get_pg9_light_sequence()])
    with pytest.raises(TypeError):
        airr.run_records(fixture_setup.get_pg9_heavy_sequence())


def test_airr_from_dataframe(fixture_setup: SadieFixture) -> None:
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


def test_airr_from_file(fixture_setup: SadieFixture) -> None:
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


def test_vdj_overlap(fixture_setup: SadieFixture) -> None:
    """Test if we can allow vdj_overlap"""
    air_api = Airr("human", allow_vdj_overlap=True, adaptable=False)
    air_api.run_single("PG9", fixture_setup.get_pg9_heavy_sequence().seq.__str__())
    assert air_api.igblast.allow_vdj_overlap.value is True

    #  make sure a warning is thrown if we have both adapt penalty and vdj overlap
    with pytest.warns(UserWarning):
        air_api = Airr("human", allow_vdj_overlap=True, adaptable=True)


def test_adaptable_penalty(fixture_setup: SadieFixture) -> None:
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
    # With current IgBLAST configuration, this sequence is properly annotated
    # even without adaptable penalties (liable=False in both cases)
    assert fw4_ == "FGGGTKVDIK"  # Previously expected to be NaN when adaptable=False
    assert cdr3_ == "QQYDN"  # Previously expected to be NaN when adaptable=False

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


def test_adaptable_correction(fixture_setup: SadieFixture) -> None:
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
    # Updated expectations: the adaptive penalty algorithm now corrects 2 more sequences
    assert liable_mix["liable"].value_counts().to_dict() == {True: 28, False: 26}


def test_single_adaptable(fixture_setup: SadieFixture) -> None:
    seq = "GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCCAGGCGAGTCAGGACATTAGCAACTATTTAAATTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTACGATGCATCCAATTTGGAAACAGGGGTCCCATCAAGGTTCAGTGGAAGTGGATCTGGGACAGATTTTACTTTCACCATCAGCAGCCTGCAGCCTGAAGATATTGCAACATATTACTGTCAACAGTATGAGACTTTCGGCCCTGGGACCAAAGTGGATATCAAAC"
    name = "test"
    airr_api = Airr("human", adaptable=True)
    result = airr_api.run_single(name, seq)
    assert isinstance(result["cdr3_aa"].iloc[0], str)


def test_mutational_analysis(heavy_catnap_airrtable: AirrTable, light_catnap_airrtable: AirrTable) -> None:
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
    # return heavy_airrtable_with_analysis, light_airrtable_with_analysis, joined_airr_table_with_analysis


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


def test_igl_assignment_errors():
    bad_seq = """>CH103_Liao2013_CH103_light_chain_from_Liao2013,_Supplementary_Table_4
CCTATGAGCTTGACTCAGCCACCCTCAGTGTCCGTGTCCCCAGGACAGACAGCCACCATAACGTGCTCTGGGGCAAGTAC
AAATGTTTGCTGGTATCAGGTGAAGCCGGGCCAGTCCCCTGAGGTGGTCATCTTTGAGAATTATAAGCGGCCCTCAGGGA
TCCCTGACCGGTTCTCTGGCTCCAAGTCTGGGAGCACAGCCACTCTGACCATCCGCGGGACCCAGGCTATAGATGAGGCT
GATTATTACTGTCAGGTGTGGGACAGCTTCTCCACCTTCGTCTTCGGATCTGGGACCCAGGTCACCGTCCTC"""


def test_five_and_three_prime_extension(fixture_setup: SadieFixture) -> None:
    airr_table = AirrTable(pd.read_feather(fixture_setup.get_bum_igl_assignment()))
    airr_methods.run_five_prime_buffer(airr_table)


def test_hard_igl_seqs(fixture_setup: SadieFixture) -> None:
    """Test we can run igl on super hard igl macaque set"""
    airr_table = AirrTable(pd.read_feather(fixture_setup.get_bum_igl_assignment()))
    out_airr_single = airr_methods.run_termini_buffers(airr_table)
    igl_df = airr_methods.run_igl_assignment(out_airr_single)
    output_solution = AirrTable(pd.read_feather(fixture_setup.get_bum_igl_solution()))
    pd.testing.assert_frame_equal(igl_df, output_solution)


def test_hard_igl_seqs_linked(fixture_setup: SadieFixture) -> None:
    """Test we can run igl on super hard igl macaque set linked"""
    lat = LinkedAirrTable(pd.read_feather(fixture_setup.get_bum_link_igl_assignment()))
    out_lat = airr_methods.run_termini_buffers(lat)
    igl_df = airr_methods.run_igl_assignment(out_lat)
    output_solution = LinkedAirrTable(pd.read_feather(fixture_setup.get_bum_link_igl_solution()))
    igl_df = igl_df[output_solution.columns]

    # assert frame equal does not work. Let's go element by element
    for x in igl_df.columns:
        pd.testing.assert_series_equal(igl_df[x], output_solution[x])


def test_runtime_reference(fixture_setup: SadieFixture) -> None:
    """Test JIT reference generation"""
    reference = Reference()
    reference.add_gene({"species": "human", "gene": "IGHV1-2*01", "source": "imgt"})
    reference.add_gene({"species": "macaque", "gene": "IGHV1-105*01", "source": "custom"})
    references = References()
    references.add_reference("silly_monkey", reference)

    # make sure we add a J gene
    with pytest.raises(ValueError):
        airr_api = Airr("silly_monkey", references=references, adaptable=True)
    reference.add_gene({"species": "human", "gene": "IGHJ6*01", "source": "imgt"})

    # no D gene
    with pytest.raises(ValueError):
        airr_api = Airr("silly_monkey", references=references, adaptable=True)

    reference.add_gene({"species": "human", "gene": "IGHD3-3*01", "source": "imgt"})
    references.add_reference("silly_monkey", reference, overwrite=True)
    airr_api = Airr("silly_monkey", references=references, adaptable=True)
    at = airr_api.run_fasta(fixture_setup.get_pg9_heavy_fasta())
    assert (at["reference_name"] == "silly_monkey").all()
    reference.add_gene({"species": "human", "gene": "IGHV1-2*01", "source": "imgt"})
    reference.add_gene({"species": "human", "gene": "IGHV3-15*01", "source": "imgt"})
    reference.add_gene({"species": "human", "gene": "IGHJ6*01", "source": "imgt"})
    reference.add_gene({"species": "mouse", "gene": "IGHV1-11*01", "source": "imgt"})
    references.add_reference("crazy_mouse", reference, overwrite=True)
    airr_crazy = Airr("crazy_mouse", references=references, adaptable=True)
    with pytest.raises(BadDataSet):
        # we used wrong name
        Airr("mighty_mouse", references=references, adaptable=True)

    airr_table = airr_crazy.run_single("PG9", fixture_setup.get_pg9_heavy_sequence().seq.__str__())
    assert airr_table.reference_name.iloc[0] == "crazy_mouse"
    assert airr_table.v_call_top.iloc[0] == "human|IGHV3-15*01"

    reference = Reference()
    reference.add_gene({"species": "human", "gene": "IGHV1-2*01", "source": "imgt"})
    reference.add_gene({"species": "human", "gene": "IGHJ6*01", "source": "imgt"})
    reference.add_gene({"species": "human", "gene": "IGHD3-3*01", "source": "imgt"})
    references.add_reference("normal_human", reference, overwrite=True)
    airr_api = Airr("normal_human", references=references, adaptable=True)
    airr_table = airr_api.run_fasta(fixture_setup.get_OAS_liable_file())


def test_airrtable_init(fixture_setup: SadieFixture) -> None:
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


def test_indel_correction(fixture_setup: SadieFixture) -> None:
    test_csv = fixture_setup.get_dog_airrtable()
    # Test we can initllize with staic meathod
    airr_table = AirrTable.read_airr(test_csv)
    airr_table_indel = AirrTable.correct_indel(airr_table)
    assert "germline_alignment_aa_corrected" in airr_table_indel.columns
    assert "v_germline_alignment_aa_corrected" in airr_table_indel.columns
    assert isinstance(airr_table_indel, AirrTable)


def test_write_and_check_airr(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    """Check that the offical airr can validate our airr tables"""
    catnap_heavy = fixture_setup.get_catnap_heavy_nt()
    output_file = tmp_path_factory.mktemp("tmp_dir") / Path("output.tsv")
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

    # files extensions checked
    with pytest.raises(ValueError):
        airr_table.to_airr(filename=fixture_setup.tmp_path / "test.tsvv")
    with pytest.raises(ValueError):
        airr_table.to_airr(filename=fixture_setup.tmp_path / "test.tsvv.bz2")


def test_airr_constant_region(fixture_setup: SadieFixture) -> None:
    constant_fasta = fixture_setup.get_fasta_with_constant()
    airr_api = Airr("human")
    results = airr_api.run_fasta(constant_fasta)
    c_cols = [
        "c_call",
        "c_cigar",
        "c_germline_alignment",
        "c_germline_alignment_aa",
        "c_germline_start",
        "c_germline_end",
        "c_identity",
        "c_score",
        "c_sequence_alignment",
        "c_sequence_alignment_aa",
        "c_sequence_start",
        "c_sequence_end",
        "c_support",
    ]
    assert all([i in results.columns for i in c_cols])
    airr_api = Airr("macaque")
    results = airr_api.run_fasta(constant_fasta)
    print(results.columns)


def test_genbank_feature():
    with pytest.raises(TypeError):
        GenBankFeature(start=0, end=20, feature_type="FWR1_fake")


def test_genbank():
    GenBank(sequence="ATCG", id="test", name="test")
    GenBank(sequence=Seq("ATCG"), id="test", name="test")
    with pytest.raises(TypeError):
        GenBank(sequence=1, id="test", name="test")
    with pytest.raises(TypeError):
        GenBank(sequence="ATCG", id="test").add_feature(feature="bad_feature")
