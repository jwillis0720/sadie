"""Low level module for IgBLAST api calls

You probably only want to interact with this module through airr as the database files are extremely tricky to get right

"""
from __future__ import annotations

import glob
import logging
import os
import subprocess
import tempfile
import warnings
from multiprocessing import cpu_count
from pathlib import Path
from typing import Any, List, Union

# Third party
import pandas as pd
import semantic_version

from sadie.airr.airrtable.constants import IGBLAST_AIRR
from sadie.airr.exceptions import (
    BadIgBLASTArgument,
    BadIgBLASTExe,
    BadIgDATA,
    IgBLASTRunTimeError,
    MissingIgBLASTArgument,
)

# package/module level
from sadie.utility.util import is_tool

# get logger in global scope
logger = logging.getLogger("IgBLAST")


def ensure_prefix_to(path: Union[str, Path]) -> Union[Path, bool]:
    """Ensure that the blast db is actually a blast like database

    The problem is that a blast_db takes in a prefix file

    ex. /path/to/blast/human_V

    which is a file or path that does not actually exists, but blasts uses it as a file glob to match
    /path/to/blast/human_V.nod
    /path/to/blast/human_V.nsq
    /path/to/blast/human_V.fasta


    because of this, we don't have validated file path but an actual file path glob. This method validates that the file glob returns blast like files

    Parameters
    ----------
    path : str
        path glob

    Returns
    -------
    return : Union[str, bool]
       returns Path or False if not a file glob
    """

    # convert Path to str
    path = str(path)
    directory_path = os.path.dirname(path)
    if not os.path.exists(directory_path):
        return False

    # get abs path to directory
    directory_path = os.path.abspath(directory_path)

    # base name
    basename = os.path.basename(path)
    glob_path = os.path.join(directory_path, basename) + "*"
    # make sure that there are things that match this queyr
    if not glob.glob(glob_path):
        return False
    return Path(os.path.join(directory_path, basename))


class IgBLASTArgument:
    """
    A class for handling all IgBLAST Arguments
    """

    def __init__(self, name: str, arg_key: str, arg_value: Union[str, int, bool, Path], required: bool):
        """IgBLASTArgument Class constructor

        Parameters
        ----------
        name : str
            the internal name for the computer
        arg_key : str
           the argument key, ex. -germline_db_V
        arg_value : Union[str, int, bool]
           the value for the argument /path/germline/db/V
        required : bool
           is the argument required
        """
        self.name = name
        self.key = arg_key
        self.value = arg_value
        self.required = required

    @property
    def name(self) -> str:
        """
        An internal name for the argument
        """
        return self._name

    @name.setter
    def name(self, n: str) -> None:
        self._name = n

    @property
    def key(self) -> str:
        """The blast command key argument

        ex. '-germline_db_v

        Returns
        -------
        str
           blast command key
        """
        return self._key

    @key.setter
    def key(self, k: str) -> None:
        self._key = k

    @property
    def value(self) -> Union[str, int, bool, Path]:
        """Return the value of the argument

        Returns
        -------
        Union[str,int,bool]
           ex. /path/to/database
        """
        return self._value

    @value.setter
    def value(self, v: Union[str, int, bool, Path]) -> None:
        self._value = v

    @property
    def required(self) -> bool:
        """Returns if the argument is required

        Returns
        -------
        bool
            if argument is required
        """
        return self._required

    @required.setter
    def required(self, r: bool) -> None:
        self._required = r

    def get_formatted_blast_arg(self) -> List[str]:
        """Return the blast formatted argument as a list

        Returns
        -------
        Union[List[str],List[str,str]]
           Either returns a single argument ['-arg'] for bool args or key value arguments ['-arg', 'value']
        """
        # If value is a bool, we return the key
        if isinstance(self.value, bool) and self.value:
            return ["-" + self.key]
        else:
            # If its not a bool, we check if it has been set
            if self.value:
                return ["-" + str(self.key), str(self.value)]
        return []  # return an empty str if its been set

    def __str__(self) -> str:
        return "{}-{}".format(self.name, self.key)


ArgumentType = Union[IgBLASTArgument, int, str, Path]


class IgBLASTN:
    """IgBLASTN

    IgBLASTN class.  A tool for immunoglobulin (IG) and T cell receptor (TR) V domain sequences from nucletodies.

    This is a lower level class and you should probably use sadie.airr to interact

    Examples
    --------
    >>> ig_blast = igblast.IgBLASTN('igblastn')
    >>> germline_ref = "reference/germlines/"
    >>> db_ref = "reference/germlines/blastdb/Ig"
    >>> aux_path = "reference/germlines/aux_data)

    # Set data
    >>> ig_blast.igdata = germline_ref

    >>> query = "fasta_inputs/PG9_H.fasta"
    >>> ig_blast.germline_db_v = os.path.join(db_ref, "human/human_V")
    >>> ig_blast.germline_db_d = os.path.join(db_ref, "human/human_D")
    >>> ig_blast.germline_db_j = os.path.join(db_ref, "human/human_J")
    >>> ig_blast.aux_path = os.path.join(aux_path, "human_gl.aux")
    >>> ig_blast.organism = "human"
    >>> csv_dataframe = ig_blast.run(query)
    """

    # Only allow these attributes
    __slots__ = [
        "_executable",
        "_version",
        "_min_d_match",
        "_num_v",
        "_num_d",
        "_num_j",
        "_outfmt",
        "_receptor",
        "_word_size",
        "_nomenclature",
        "_gap_open",
        "_gap_extend",
        "_num_threads",
        "_extend_5",
        "_extend_3",
        "_j_penalty",
        "_v_penalty",
        "_d_penalty",
        "_organism",
        "_germline_db_v",
        "_germline_db_d",
        "_germline_db_j",
        "_germline_db_c",
        "_aux_path",
        "_igdata",
        "_temp_dir",
        "_allow_vdj_overlap",
        "debug",
    ]

    def __init__(self, executable: Union[Path, str] = "igblastn", tmp_dir: Union[Path, str] = None):
        """IgBLASTN with a query. Set everything up with a setter

        Parameters
        ----------
        excutable : Union[Path,str], optional
            Path to the igblastn executable, by default "igblastn"
        tmp_dir : Union[Path,str], optional
            Path to the temporary directory, by default None
        """

        # set the executable dynamically
        self.executable = Path(executable)
        self._version = self._get_version()

        # setup all the default values if we don't add them
        self.min_d_match = 5
        self.num_v = 3
        self.num_d = 3
        self.num_j = 3
        self.outfmt = 19
        self.receptor = "Ig"
        self.word_size = 5
        self.nomenclature = "imgt"
        self.gap_open = 5
        self.gap_extend = 2
        self.num_threads = cpu_count()
        self.extend_5 = True
        self.extend_3 = True
        self.j_penalty = -2
        self.v_penalty = -1
        self.d_penalty = -1
        self.allow_vdj_overlap = False

        # Make these blank, if they are not set by the caller, then we will complain during runtime. They must be set dynamically
        self._organism = IgBLASTArgument("organism", "organism", "", True)
        self._germline_db_v = IgBLASTArgument("germline_db_v", "germline_db_V", "", True)
        self._germline_db_d = IgBLASTArgument("germline_db_d", "germline_db_D", "", True)
        self._germline_db_j = IgBLASTArgument("germline_db_j", "germline_db_J", "", True)
        self._germline_db_c = IgBLASTArgument("c_region_db", "c_region_db", "", True)
        self._aux_path = IgBLASTArgument("aux_path", "auxiliary_data", "", True)

        # Igdata is not an official blast argument, it is an enviroment
        self._igdata = Path(".")
        self.temp_dir = tmp_dir or Path(".")

        # Debug mode - prints commands before execution
        self.debug = False

    def _get_version(self) -> semantic_version.Version:
        """Private method to parse igblast -version and get semantic_version

        Returns
        -------
        semantic_version.Version
            the igblast version
        """
        try:
            process = subprocess.run([self.executable, "-version"], capture_output=True)
            stdout = process.stdout.decode("utf-8")
        except (FileNotFoundError, subprocess.CalledProcessError) as e:
            raise BadIgBLASTExe(self.executable, f"igblastn -version failed: {e}")
        version = stdout.split("\n")[0].split(":")[-1].strip()
        try:
            version = semantic_version.Version(version)
        except ValueError:
            raise BadIgBLASTExe(self.executable, f"semantic version can't parse {stdout}.")
        return version

    @property
    def version(self) -> semantic_version.Version:
        return self._version

    @property
    def executable(self) -> Path:
        """The igblastn executable

        Returns
        -------
        Path
            igblastn path
        """
        return self._executable

    @executable.setter
    def executable(self, exe: Path) -> None:
        self._executable = exe

    @property
    def temp_dir(self) -> Path:
        """The path to the tempdata directory for spliting igblast

        Returns
        -------
        Path
           A valid temporary directory path
        """
        return self._temp_dir

    @temp_dir.setter
    def temp_dir(self, data: Union[Path, str]) -> None:
        if isinstance(data, str):
            data = Path(data)
        self._temp_dir = data.absolute()
        if not os.access(self._temp_dir, os.W_OK):
            raise IOError(self._temp_dir, "Unable to write to temp dir")

    @property
    def igdata(self) -> Path:
        """The path to IGDATA which contains the internal_data needed to make a recomdination

        Returns
        -------
        Path
           A valid IGDATA path
        """
        return Path(self._igdata)

    @igdata.setter
    def igdata(self, data: Union[Path, str]) -> None:
        if isinstance(data, str):
            data = Path(data)
        if not data.exists() or not data.is_dir():
            raise BadIgDATA(data)
        self._igdata = Path(data.absolute())

    @property
    def min_d_match(self) -> ArgumentType:
        """Required minimal consecutive nucleotide base matches for D genes

        Returns
        -------
        IgBLASTArgument
        """
        return self._min_d_match

    @min_d_match.setter
    def min_d_match(self, d: int) -> None:
        """Required minimal consecutive nucleotide base matches for D genes

        Parameters
        ----------
        d : int
            The number of matches

        Raises
        ------
        BadIgBLASTArgument
            Must be >= 5
        """
        if d < 5:
            raise BadIgBLASTArgument(d, ">= 5")
        self._min_d_match = IgBLASTArgument("min_d_match", "min_D_match", d, False)

    @property
    def num_v(self) -> ArgumentType:
        """
           Number of Germline sequences to show alignments for

        Returns
        -------
        IgBLASTArgument
        """
        return self._num_v

    @num_v.setter
    def num_v(self, v: int) -> None:
        self._num_v = IgBLASTArgument("num_v", "num_alignments_V", v, False)

    @property
    def num_d(self) -> ArgumentType:
        """
           Number of Germline sequences to show alignments for D gene

        Returns
        -------
        IgBLASTArgument
        """
        return self._num_d

    @num_d.setter
    def num_d(self, d: int) -> None:
        self._num_d = IgBLASTArgument("num_d", "num_alignments_D", d, False)

    @property
    def num_j(self) -> ArgumentType:
        """
           Number of Germline sequences to show alignments for J gene

        Returns
        -------
        IgBLASTArgument
        """
        return self._num_j

    @num_j.setter
    def num_j(self, j: int) -> None:
        self._num_j = IgBLASTArgument("num_j", "num_alignments_J", j, False)

    @property
    def organism(self) -> ArgumentType:
        """The organism for your query sequence.

        Returns
        -------
        IgBLASTArgument
        """
        return self._organism

    @organism.setter
    def organism(self, o: str) -> None:
        """Organism

        Parameters
        ----------
        o : str
            an organism string
        """
        # I don't want to hardcode in the organisms here.
        # I will handle that logic at a higher level,
        # this is because blast has no preset organims and it's all about the v,d,j blast paths which are set dynamically
        self._organism = IgBLASTArgument("organism", "organism", o, True)

    @property
    def outfmt(self) -> ArgumentType:
        """alignment view options:
            3 = Flat query-anchored, show identities,
            4 = Flat query-anchored, no identities,
            7 = Tabular with comment lines
            19 = Rearrangement summary report (AIRR format)

        Returns
        -------
        IgBLASTArgument
        """
        return self._outfmt

    @outfmt.setter
    def outfmt(self, fmt: int) -> None:
        """alignment view options:
            3 = Flat query-anchored, show identities,
            4 = Flat query-anchored, no identities,
            7 = Tabular with comment lines
            19 = Rearrangement summary report (AIRR format)

        Exceptions
        ----------
        BadIgBLASTArgument
            Only 19 is supported for now
        """
        # only accept 19 for now
        if fmt != 19:
            raise BadIgBLASTArgument(fmt, 19)
        self._outfmt = IgBLASTArgument("outfmt", "outfmt", fmt, True)

    @property
    def receptor(self) -> ArgumentType:
        """Specify Ig or T cell receptor sequence

        Returns
        -------
        IgBLASTArgument
        """
        return self._receptor

    @receptor.setter
    def receptor(self, r: str) -> None:
        """Specify Ig or T cell receptor sequence

        Exceptions
        ----------
        BadIgBLASTArgument
            Only 'Ig' or 'TCR' are supported
        """
        if r not in ["Ig", "TCR"]:
            raise BadIgBLASTArgument(r, ["Ig", "TCR"])
        self._receptor = IgBLASTArgument("receptor", "ig_seqtype", r, True)

    @property
    def nomenclature(self) -> ArgumentType:
        """Domain system to be used for segment annotation

        Returns
        -------
        IgBLASTArgument
        """
        return self._nomenclature

    @nomenclature.setter
    def nomenclature(self, system: str) -> None:
        """Domain system to be used for segment annotation

        Exceptions
        ----------
        BadIgBLASTArgument
            Only 'imgt' or 'kabat' are supported
        """
        if system.lower() not in ["imgt", "kabat"]:
            raise BadIgBLASTArgument(system, "['imgt','kabat']")
        self._nomenclature = IgBLASTArgument("nomenclature", "domain_system", system, True)

    @property
    def aux_path(self) -> ArgumentType:
        """Auxilary data path. This is needed to lookup the J genes and tell them when the CDR3 stops.

        Returns
        -------
        IgBLASTArgument
        """
        return self._aux_path

    @aux_path.setter
    def aux_path(self, aux_path: Path | str) -> None:
        if isinstance(aux_path, str):
            aux_path = Path(aux_path)
        if not aux_path.exists():
            raise BadIgBLASTArgument(aux_path, "valid path to Auxilary database")
        self._aux_path = IgBLASTArgument("aux_path", "auxiliary_data", aux_path.absolute(), True)

    @property
    def germline_db_v(self) -> ArgumentType:
        """Path to V gene database prefix

        Returns
        -------
        IgBLASTArgument
        """
        return self._germline_db_v

    @germline_db_v.setter
    def germline_db_v(self, path: str | Path) -> None:
        """Path to V gene database prefix

        Exceptions
        ----------
        BadIgBLASTArgument
            Only valid path to V Database is supported
        """
        abs_path = ensure_prefix_to(path)
        if not abs_path:
            raise BadIgBLASTArgument(path, "Valid path to V Database")
        self._germline_db_v = IgBLASTArgument("germline_db_v", "germline_db_V", path, True)

    @property
    def germline_db_d(self) -> ArgumentType:
        """Path to D gene database prefix

        Returns
        -------
        IgBLASTArgument
        """
        return self._germline_db_d

    @germline_db_d.setter
    def germline_db_d(self, path: str | Path) -> None:
        """Path to D gene database prefix"""
        abs_path = ensure_prefix_to(path)
        if not abs_path:
            warnings.warn(f"{path} is not found, No D gene segment", UserWarning)
            # raise BadIgBLASTArgument(path, "Valid path to D Database")
            self._germline_db_d = IgBLASTArgument("germline_db_d", "germline_db_D", "", False)
        else:
            self._germline_db_d = IgBLASTArgument("germline_db_d", "germline_db_D", path, True)

    @property
    def germline_db_c(self) -> ArgumentType:
        """Path to C gene database prefix

        Returns
        -------
        IgBLASTArgument
        """
        return self._germline_db_c

    @germline_db_c.setter
    def germline_db_c(self, path: str | Path) -> None:
        """Path to C gene database prefix"""
        abs_path = ensure_prefix_to(path)
        if not abs_path:
            warnings.warn(f"{path} is not found, No C gene segment", UserWarning)
            self._germline_db_c = IgBLASTArgument("c_region_db", "c_region_db", "", False)
        else:
            self._germline_db_c = IgBLASTArgument("c_region_db", "c_region_db", path, True)

    @property
    def germline_db_j(self) -> ArgumentType:
        """Path to J gene database prefix

        Returns
        -------
        IgBLASTArgument
        """
        return self._germline_db_j

    @germline_db_j.setter
    def germline_db_j(self, path: str | Path) -> None:
        """Path to J gene database prefix

        Exceptions
        ----------
        BadIgBLASTArgument
            Only valid path to J Database is supported
        """
        abs_path = ensure_prefix_to(path)
        if not abs_path:
            raise BadIgBLASTArgument(path, "Valid path to J Database")
        self._germline_db_j = IgBLASTArgument("germline_db_j", "germline_db_J", path, True)

    @property
    def word_size(self) -> ArgumentType:
        """Word size for wordfinder algorithm (length of best perfect match)

        Returns
        -------
        IgBLASTArugment
        """
        return self._word_size

    @word_size.setter
    def word_size(self, word_size: int) -> None:
        """Word size for wordfinder algorithm (length of best perfect match)

        Exceptions
        ----------
        BadIgBLASTArgument
            Only integer greater than 4 is supported
        """
        if word_size < 4:
            raise BadIgBLASTArgument(word_size, ">=4")
        self._word_size = IgBLASTArgument("word_size", "word_size", word_size, False)

    @property
    def gap_open(self) -> ArgumentType:
        """Cost to open a gap

        Returns
        -------
        IgBLASTArgument
        """
        return self._gap_open

    @gap_open.setter
    def gap_open(self, go: int) -> None:
        """Cost to open a gap

        Exceptions
        ----------
        BadIgBLASTArgument
            Only integer greater than 0 is supported
        """
        if go < 0:
            raise BadIgBLASTArgument(go, ">=0")
        self._gap_open = IgBLASTArgument("gap_open", "gapopen", go, False)

    @property
    def gap_extend(self) -> ArgumentType:
        """Cost to extend a gap

        Returns
        -------
        IgBLASTArgument
        """
        return self._gap_extend

    @gap_extend.setter
    def gap_extend(self, ge: int) -> None:
        """Cost to extend a gap

        Exceptions
        ----------
        BadIgBLASTArgument
            Only integer greater than 0 is supported
        """
        if ge < 0:
            raise BadIgBLASTArgument(ge, ">=0")
        self._gap_extend = IgBLASTArgument("gap_open", "gapextend", ge, False)

    @property
    def num_threads(self) -> ArgumentType:
        """
            Number of threads (CPUs) to use in the BLAST search

        Returns
        -------
        IgBLASTArgument
        """
        return self._num_threads

    @num_threads.setter
    def num_threads(self, num_threads: int) -> None:
        # Some architectures have more threads than number of cores. Letting the users decide the number of threads.
        self._num_threads = IgBLASTArgument("number_threds", "num_threads", num_threads, False)

    @property
    def extend_5(self) -> ArgumentType:
        """Extend V gene alignment at 5' end

        Returns
        -------
        IgBLASTArgument
        """
        return self._extend_5

    @extend_5.setter
    def extend_5(self, extend_5: bool) -> None:
        self._extend_5 = IgBLASTArgument("extend_5", "extend_align5end", extend_5, False)

    @property
    def extend_3(self) -> ArgumentType:
        """Extend V gene alignment at 3' end

        Returns
        -------
        IgBLASTArgument
        """
        return self._extend_3

    @extend_3.setter
    def extend_3(self, extend_3: bool) -> None:
        self._extend_3 = IgBLASTArgument("extend_3", "extend_align3end", extend_3, False)

    @property
    def allow_vdj_overlap(self) -> Any:
        """Allow the VDJ overlap

        This option is active only when D_penalty
        and J_penalty are set to -4 and -3, respectively

        Returns
        -------
        IgBLASTArgument
        """
        return self._allow_vdj_overlap  # type: ignore[has-type]

    @allow_vdj_overlap.setter
    def allow_vdj_overlap(self, allow: bool) -> None:
        """Allow the VDJ overlap

        This option is active only when D_penalty
        and J_penalty are set to -4 and -3, respectively
        """
        j_penalty: IgBLASTArgument = self.j_penalty  # type: ignore[assignment]
        d_penalty: IgBLASTArgument = self.d_penalty  # type: ignore[assignment]
        if j_penalty.value != -3 and d_penalty.value != -4 and allow:
            warnings.warn(
                f"Allows vdj overlap set but j penalty and d penalty need to be -3 and -4, now are {self.j_penalty}, {self.d_penalty}",
                UserWarning,
            )
        self._allow_vdj_overlap = IgBLASTArgument("allow_vdj_overlap", "allow_vdj_overlap", allow, False)

    @property
    def d_penalty(self) -> ArgumentType:
        """What is the  D gene panalty

        Returns
        -------
        IgBLASTArgument
        """
        return self._d_penalty

    @d_penalty.setter
    def d_penalty(self, penalty: int) -> None:
        """What is the  D gene panalty

        Exceptions
        ----------
        BadIgBLASTArgument
            Only integer less than 0 and greater than -5 is supported
        """
        if not -5 < penalty < 1:
            raise BadIgBLASTArgument(penalty, "must be less than 0 and greater than -5")
        self._d_penalty = IgBLASTArgument("d_penalty", "D_penalty", penalty, True)

    @property
    def j_penalty(self) -> ArgumentType:
        """What is the  J gene panalty

        Returns
        -------
        IgBLASTArgument
        """
        return self._j_penalty

    @j_penalty.setter
    def j_penalty(self, penalty: int) -> None:
        """What is the  J gene panalty

        Exceptions
        ----------
        BadIgBLASTArgument
            Only integer less than 0 and greater than -4 is supported
        """
        if not -4 < penalty < 1:
            raise BadIgBLASTArgument(penalty, "must be less than 0 and greater than -4")
        self._j_penalty = IgBLASTArgument("j_penalty", "J_penalty", penalty, True)

    @property
    def v_penalty(self) -> ArgumentType:
        """What is the  v gene panalty

        Returns
        -------
        IgBLASTArgument
        """
        return self._v_penalty

    @v_penalty.setter
    def v_penalty(self, penalty: int) -> None:
        if not -5 < penalty < 1:
            raise BadIgBLASTArgument(penalty, "must be less than 0 and greater than -5")
        self._v_penalty = IgBLASTArgument("v_penalty", "V_penalty", penalty, True)

    @property
    def arguments(self) -> List[IgBLASTArgument]:
        """return a list of IgBLASTArugments

        Returns
        -------
        List[IgBLASTArguments]
        """
        # lots of type ignores since these are IgBLASTArguments set in the setter, but are read from the property
        return [
            self.min_d_match,  # type: ignore
            self.num_v,  # type: ignore
            self.num_j,  # type: ignore
            self.num_d,  # type: ignore
            self.organism,  # type: ignore
            self.receptor,  # type: ignore
            self.germline_db_v,  # type: ignore
            self.germline_db_d,  # type: ignore
            self.germline_db_j,  # type: ignore
            self.germline_db_c,  # type: ignore
            self.aux_path,  # type: ignore
            self.outfmt,  # type: ignore
            self.nomenclature,  # type: ignore
            self.word_size,  # type: ignore
            self.gap_open,  # type: ignore
            self.gap_extend,  # type: ignore
            self.j_penalty,  # type: ignore
            self.v_penalty,  # type: ignore
            self.d_penalty,  # type: ignore
            self.num_threads,  # type: ignore
            self.extend_5,  # type: ignore
            self.extend_3,  # type: ignore
            self.allow_vdj_overlap,
        ]

    @property
    def cmd(self) -> List[str]:
        """Return the blast cmd that will be run by subprocess"""
        _cmd = [str(self.executable)]
        for blast_arg in self.arguments:
            kv = blast_arg.get_formatted_blast_arg()  # can return non if we already set it twice
            if kv:  # only set on boolean if they are true
                _cmd += kv
        return _cmd

    def pre_check(self) -> None:
        """Ensures we have set everything right

        Raises
        ------
        MissingIgBLASTArg
            We have set the IGDATA field
        BadIgBLASTExe
            Correct IGblast executable
        BadIgDATA
            If any of the fields are not set properly
        """
        # TODO: move to pydanic model
        # Ensure required arguments werer set
        for blast_arg in self.arguments:
            if blast_arg.required and not (blast_arg.value):
                raise MissingIgBLASTArgument(f"Missing Blast argument. Need to set IgBLASTN.{blast_arg.name}")

    # Run methods
    def run_file(self, file: Union[Path, str]) -> pd.DataFrame:
        """Run IgBlast on a file

        Parameters
        ----------
        file : Path
            the fasta file path

        Returns
        -------
        pd.DataFrame
            A dataframe with the IgBLAST results

        Raises
        ------
        IgBLASTRunTimeError
           for any given runtime error for igblastn
        """

        # because igblast uses IGDATA as the internal file structure, we should pass the enviroment to the subprocess
        local_env = os.environ.copy()
        local_env["IGDATA"] = str(self.igdata)

        # take the cmd and finally add the query file
        cmd = self.cmd
        cmd += ["-query", str(file)]

        # run a precheck to make sure everything passed was working
        self.pre_check()

        # Print debug info if debug mode is enabled
        if self.debug:
            # Clean up the command for display by resolving paths
            display_cmd = []
            i = 0
            while i < len(cmd):
                if i < len(cmd) - 1 and cmd[i].startswith("-") and not cmd[i + 1].startswith("-"):
                    # This is a parameter with a value
                    display_cmd.append(cmd[i])
                    # Check if the value looks like a path
                    if "/" in cmd[i + 1] or "\\" in cmd[i + 1]:
                        try:
                            # Try to resolve the path to remove ".."
                            resolved_path = str(Path(cmd[i + 1]).resolve())
                            display_cmd.append(resolved_path)
                        except:
                            # If path resolution fails, use original
                            display_cmd.append(cmd[i + 1])
                    else:
                        display_cmd.append(cmd[i + 1])
                    i += 2
                else:
                    display_cmd.append(cmd[i])
                    i += 1

            logger.info(f"IgBLAST command: {' '.join(display_cmd)}")
            logger.info(f"IGDATA environment variable: {self.igdata.resolve()}")

        # while we can certainly do this as an output stream on stdout,
        # It's probably best to take advantage of IGblast output and tempfile
        with tempfile.NamedTemporaryFile(dir=self.temp_dir, suffix="_igblast.tsv") as tmpfile:
            cmd += ["-out", tmpfile.name]
            process = subprocess.run(cmd, env=local_env, capture_output=True)
            if process.stderr:
                raise IgBLASTRunTimeError(process.stderr)
            # we read the dataframe from the tempfile, it should always be in .TSV.
            # We can also cast it to IGBLAST_AIRR dtypes to save memory
            df = pd.read_csv(tmpfile.name, sep="\t", dtype=IGBLAST_AIRR)  # type: ignore

        df["v_identity"] = df["v_identity"] / 100
        df["d_identity"] = df["d_identity"] / 100
        df["j_identity"] = df["j_identity"] / 100
        return df

    def __repr__(self) -> str:
        return "IgBLAST: env IGDATA={} {}".format(str(self.igdata), " ".join(self.cmd))

    def __str__(self) -> str:
        return self.__repr__()
