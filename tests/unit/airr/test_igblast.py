import os
import platform
import re
from pathlib import Path

import pytest
import semantic_version

from sadie.airr import __file__ as sadie_airr_file
from sadie.airr import exceptions as airr_exceptions
from sadie.airr import igblast
from sadie.airr.exceptions import (
    BadIgBLASTArgument,
    BadIgBLASTExe,
    BadIgDATA,
    IgBLASTRunTimeError,
    MissingIgBLASTArgument,
)
from sadie.airr.igblast.igblast import IgBLASTArgument, IgBLASTN
from tests.conftest import SadieFixture


def test_antibody_igblast_setup() -> None:
    """
    testing if we can manually setup IgBlast databases
    """

    system = platform.system().lower()
    executable = os.path.join(os.path.dirname(sadie_airr_file), f"bin/{system}/igblastn")

    ig_blast = igblast.IgBLASTN(executable)
    # Check that version is a valid semantic version (e.g., 1.19.0, 1.22.0, etc.)
    assert isinstance(ig_blast.version, semantic_version.Version)
    # Ensure it's at least version 1.22.0 or higher
    assert ig_blast.version >= semantic_version.Version("1.22.0")
    germline_ref = os.path.join(os.path.dirname(os.path.abspath(sadie_airr_file)), "data/germlines")
    assert os.path.exists(germline_ref)

    # Set data
    for name in ["human", "mouse", "rat", "dog"]:
        db_ref = os.path.join(germline_ref, f"Ig/blastdb/{name}")
        internal_ref = os.path.join(germline_ref, f"Ig/internal_data/{name}")
        aux_ref = os.path.join(germline_ref, "aux_db/imgt/")
        assert os.path.exists(db_ref)
        ig_blast.germline_db_v = os.path.join(db_ref, "{}_V".format(name))
        ig_blast.germline_db_d = os.path.join(db_ref, "{}_D".format(name))
        ig_blast.germline_db_j = os.path.join(db_ref, "{}_J".format(name))
        ig_blast.germline_db_c = os.path.join(db_ref, "{}_C".format(name))

        aux_ref = os.path.join(aux_ref, f"{name}_gl.aux")
        ig_blast.aux_path = aux_ref
        ig_blast.organism = name
        ig_blast.igdata = internal_ref
        ig_blast.pre_check()
    assert os.path.exists(germline_ref)

    # bad germeline ref
    germline_ref = "reference/germlines/"

    # Set data
    with pytest.raises(airr_exceptions.BadIgDATA):
        ig_blast.igdata = germline_ref
        ig_blast.pre_check()
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


def test_antibody_igblast_file_run(fixture_setup: SadieFixture) -> None:
    """
    testing if we can manually run igblast withou airr abstraction
    """

    # this should be fixture
    system = platform.system().lower()
    executable = os.path.join(os.path.dirname(sadie_airr_file), f"bin/{system}/igblastn")

    ig_blast = igblast.IgBLASTN(executable)
    germline_ref = os.path.join(os.path.dirname(os.path.abspath(sadie_airr_file)), "data/germlines")
    db_ref = os.path.join(germline_ref, "Ig/blastdb/human/")
    aux_path = os.path.join(germline_ref, "aux_db/imgt/")

    # Set data
    # internal_data/human/human_V must be discoverd by IgBLAST
    ig_blast.igdata = os.path.join(germline_ref, "Ig/")

    # Grab from fixtures
    query = fixture_setup.get_pg9_heavy_fasta()
    ig_blast.germline_db_v = os.path.join(db_ref, "human_V")
    ig_blast.germline_db_d = os.path.join(db_ref, "human_D")
    ig_blast.germline_db_j = os.path.join(db_ref, "human_J")
    ig_blast.germline_db_c = os.path.join(db_ref, "human_C")
    ig_blast.aux_path = os.path.join(aux_path, "human_gl.aux")
    ig_blast.organism = "human"
    ig_blast.pre_check()

    # have to make this -1 to get a more specifc j gene match
    ig_blast.j_penalty = -1
    csv_dataframe = ig_blast.run_file(query)

    query = fixture_setup.get_pg9_heavy_fasta()
    csv_dataframe = ig_blast.run_file(query)
    csv_dataframe["d_call"] = csv_dataframe["d_call"].fillna("")


def test_IgBLASTArgument_exceptions() -> None:
    igblast_argument = IgBLASTArgument(name="name", arg_key="key", arg_value="value", required=True)
    assert igblast_argument.name == "name"
    assert igblast_argument.key == "key"
    assert igblast_argument.value == "value"
    assert igblast_argument.required is True
    assert str(igblast_argument) == "name-key"


class TestIgBLASTN:
    igblastn_exe = os.path.join(os.path.dirname(sadie_airr_file), f"bin/{platform.system().lower()}/igblastn")
    germline_ref = os.path.join(os.path.dirname(os.path.abspath(sadie_airr_file)), "data/germlines")
    db_ref = os.path.join(germline_ref, "Ig/blastdb/human/")
    aux_path = os.path.join(germline_ref, "aux_db/imgt/")

    def setup_class(self) -> None:
        self.igblastn = IgBLASTN(self.igblastn_exe)

    def get_populated_igblastn(self) -> None:
        """Assumes properies have already been tested"""
        igblastn = IgBLASTN(self.igblastn_exe)
        igblastn.igdata = os.path.join(self.germline_ref, "Ig/")
        igblastn.germline_db_v = os.path.join(self.db_ref, "human_V")
        igblastn.germline_db_d = os.path.join(self.db_ref, "human_D")
        igblastn.germline_db_j = os.path.join(self.db_ref, "human_J")
        igblastn.germline_db_c = os.path.join(self.db_ref, "human_C")
        igblastn.aux_path = os.path.join(self.aux_path, "human_gl.aux")
        igblastn.organism = "human"
        return igblastn

    def test_init(self) -> None:
        # executable poxpath
        IgBLASTN(Path(self.igblastn_exe))
        # executable str path
        IgBLASTN(str(self.igblastn_exe))
        # str tmp output path
        IgBLASTN(self.igblastn_exe, tmp_dir=".")
        # poxpath tmp output path
        IgBLASTN(self.igblastn_exe, tmp_dir=Path("."))

    def test_init_exceptions(self) -> None:
        # Path does not exist or is not executable
        with pytest.raises(BadIgBLASTExe, match="Cant find IgBLAST"):
            igblastn = IgBLASTN("igblastn-fake")
            igblastn._get_version()  # type: ignore
        # Path is executable but not igblast
        with pytest.raises(BadIgBLASTExe, match="Cant find IgBLAST"):
            igblastn = IgBLASTN("robot")
            igblastn._get_version()  # type: ignore
        # Bad tmp output path
        with pytest.raises(IOError):
            IgBLASTN(self.igblastn_exe, tmp_dir=Path("fake"))

    def test_str_repr(self) -> None:
        assert str(self.igblastn).startswith("IgBLAST")

    def test_version(self) -> None:
        assert re.match("[0-9]{1,}.[0-9]{1,}.[0-9]{1,}", str(self.igblastn.version)) is not None

    def test_min_d_match(self) -> None:
        self.igblastn.min_d_match = 5
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.min_d_match = 0

    def test_outfmt(self) -> None:
        self.igblastn.outfmt = 19
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.outfmt = 0

    def test_receptor(self) -> None:
        self.igblastn.receptor = "Ig"
        self.igblastn.receptor = "TCR"
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.receptor = "fake"

    def test_nomenclature(self) -> None:
        self.igblastn.nomenclature = "imgt"
        self.igblastn.nomenclature = "kabat"
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.nomenclature = "fake"

    def test_word_size(self) -> None:
        self.igblastn.word_size = 4
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.word_size = 0

    def test_gap_open(self) -> None:
        self.igblastn.gap_open = 0
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.gap_open = -1

    def test_gap_extend(self) -> None:
        self.igblastn.gap_extend = 0
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.gap_extend = -1

    def test_allow_vdj_overlap(self) -> None:
        self.igblastn.v_penalty = -1
        self.igblastn.d_penalty = -1
        with pytest.warns(UserWarning):
            self.igblastn.allow_vdj_overlap = True
        self.igblastn.allow_vdj_overlap = False

    def test_v_penalty(self) -> None:
        self.igblastn.v_penalty = -4
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.v_penalty = 1

    def test_d_penalty(self) -> None:
        self.igblastn.d_penalty = -3
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.d_penalty = 1

    def test_j_penalty(self) -> None:
        self.igblastn.j_penalty = -3
        with pytest.raises(BadIgBLASTArgument, match="Passed argument"):
            self.igblastn.j_penalty = 1

    def test_igdata(self) -> None:
        with pytest.raises(BadIgDATA, match="Bad IgDATA"):
            self.igblastn.igdata = "fake"

    def test_pre_check(self) -> None:
        with pytest.raises(MissingIgBLASTArgument, match="Missing"):
            self.igblastn.pre_check()

    def test_run_file_exceptions(self) -> None:
        igblastn = self.get_populated_igblastn()
        with pytest.raises(IgBLASTRunTimeError, match="Runtime Error"):
            igblastn.run_file(Path("fake"))
