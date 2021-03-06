"""Unit tests for antibody."""
import logging
import os
import platform

# third party
import pytest
import semantic_version

# package level
from sadie import airr
from sadie.airr.exceptions import BadIgBLASTExe

loger = logging.getLogger()


def get_file(file):
    """Helper method for test execution."""
    _file = os.path.join(os.path.abspath(os.path.dirname(__file__)), f"fixtures/{file}")
    if not os.path.exists(_file):
        raise FileNotFoundError(_file)
    return _file


def assert_df(expected, actual):
    assert len(expected) == len(actual)
    for i, j in zip(expected, actual):
        if isinstance(i, float):
            assert pytest.approx(i) == pytest.approx(j)
        else:
            assert i == j


def test_antibody_igblast_setup():
    """
    testing if we can manually setup IgBlast databases
    """
    system = platform.system().lower()
    executable = os.path.join(os.path.dirname(airr.__file__), f"bin/{system}/igblastn")

    ig_blast = airr.igblast.IgBLASTN(executable)
    assert ig_blast.version == semantic_version.Version("1.17.1")
    germline_ref = os.path.join(os.path.dirname(os.path.abspath(airr.__file__)), "data/germlines")
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
    with pytest.raises(airr.exceptions.BadIgDATA):
        ig_blast.igdata = germline_ref
        ig_blast._pre_check()
    for species in ["human", "mouse", "rat", "dog"]:
        db_ref = os.path.join(germline_ref, "imgt/Ig/blastdb/")
        internal_ref = os.path.join(germline_ref, "imgt/Ig/")
        aux_ref = os.path.join(germline_ref, "imgt/aux_db/")
        with pytest.raises(airr.exceptions.BadIgBLASTArgument):
            ig_blast.germline_db_v = os.path.join(db_ref, "{}_V".format(species))
        with pytest.warns(UserWarning):
            ig_blast.germline_db_d = os.path.join(db_ref, "{}_D".format(species))
        with pytest.raises(airr.exceptions.BadIgBLASTArgument):
            ig_blast.germline_db_j = os.path.join(db_ref, "{}_J".format(species))
        with pytest.raises(airr.exceptions.BadIgBLASTArgument):
            ig_blast.aux_path = aux_ref

    with pytest.raises(BadIgBLASTExe):
        # pass an executable that is not igblast
        airr.igblast.IgBLASTN("ls")


def test_antibody_igblast_file_run():
    """
    testing if we can manually run igblast
    """

    # this should be fixture
    system = platform.system().lower()
    executable = os.path.join(os.path.dirname(airr.__file__), f"bin/{system}/igblastn")

    ig_blast = airr.igblast.IgBLASTN(executable)
    germline_ref = os.path.join(os.path.dirname(os.path.abspath(airr.__file__)), "data/germlines")
    db_ref = os.path.join(germline_ref, "imgt/Ig/blastdb/")
    aux_path = os.path.join(germline_ref, "imgt/aux_db/")

    # Set data
    ig_blast.igdata = os.path.join(germline_ref, "imgt/Ig")

    # Grab from fixtures
    query = get_file("fasta_inputs/PG9_H.fasta")
    ig_blast.germline_db_v = os.path.join(db_ref, "human_V")
    ig_blast.germline_db_d = os.path.join(db_ref, "human_D")
    ig_blast.germline_db_j = os.path.join(db_ref, "human_J")
    ig_blast.aux_path = os.path.join(aux_path, "human_gl.aux")
    ig_blast.organism = "human"
    ig_blast._pre_check()
    # have to make this -1 to get a more specifc j gene match
    ig_blast.j_penalty = -1
    csv_dataframe = ig_blast.run_file(query).fillna("")

    query = get_file("fasta_inputs/PG9_L.fasta")
    csv_dataframe = ig_blast.run_file(query)
    csv_dataframe["d_call"] = csv_dataframe["d_call"].fillna("")
