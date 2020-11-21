"""Low level module for IgBLAST api calls

You probably only want to interact with this module through airr as the database files are extremely tricky to get right

"""
import glob
import logging
import os
import subprocess
import sys
import tempfile
from io import StringIO
from pathlib import Path
from typing import List, Union
import platform
import shutil

import pandas as pd

from .util import is_tool

logger = logging.getLogger("IgBLAST")


def ensure_prefix_to(path: str):
    """Ensure that the blast db is actually a balst like database

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
    return Path
       returns Path or none if not a file glob
    """

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
    return os.path.join(directory_path, basename)


class Error(Exception):
    """Base class for exceptions in this module."""


class BadIgBLASTExe(Error):
    """Exception raised for not finiding the igblast module

    Attributes:
    """

    def __init__(self, passed_executable, msg):
        super().__init__()
        self.passed_arguments = passed_executable
        self.msg = msg

    def __str__(self):
        _env = os.environ["PATH"]
        return f"Cant find IgBLAST {self.passed_arguments}. Check {_env}\n {self.msg}"


class EmtpyFileError(Error):
    """Exception raised for a file passed to igblast being empty

    this is needed because blasts accepts empty files but we will not

    """

    def __init__(self, file):
        super().__init__()
        self.passed_arguments = file

    def __str__(self):
        return "{} is empty".format(self.passed_arguments)


class MissingIgBLASTArg(Error):
    """Missing a required IgBLAST command argument

    If a required command is missing

    """

    def __init__(self, msg):
        super().__init__()
        self.msg = msg

    def __str__(self):
        return self.msg


class BadIgBLASTArgument(Error):
    """Exception raised for passing incorrect params to an igblast arguments"""

    def __init__(self, passed_arguments, accepted_argumetns):
        super().__init__()
        self.passed_arguments = passed_arguments
        self.accepted_arguments = accepted_argumetns

    def __str__(self):
        return "Passed argument {}. Only accepts {}".format(self.passed_arguments, self.accepted_arguments)


class BadIgDATA(Error):
    """Exception raised for IgData path (which is crucial) not being found

    Attributes:
    """

    def __init__(self, passed_arguments):
        super().__init__()
        self.passed_arguments = passed_arguments

    def __str__(self):
        return f"Bad IgDAta path {self.passed_arguments} - please provide where IgDATA is located"


class IgBLASTRunTimeError(Error):
    """Exception raised for Igblast runtime error"""

    def __init__(self, stderr):
        super().__init__()
        self.stderr = stderr

    def __str__(self):
        return "Runtime Error with Blast {}".format(self.stderr.decode("utf-8"))


class IgBLASTArgument:
    """
    A class for handling all IgBLAST Arguments
    """

    def __init__(self, name: str, arg_key: str, arg_value: Union[str, int, bool], required: bool):
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
    def name(self, n: str):
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
    def key(self, k: str):
        self._key = k

    @property
    def value(self) -> Union[str, int, bool]:
        """Return the value of the argument

        Returns
        -------
        Union[str,int,bool]
           ex. /path/to/database
        """
        return self._value

    @value.setter
    def value(self, v):
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
    def required(self, r):
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
        return False

    def __str__(self):
        return "{}-{}".format(self.name, self.key)


class IgBLASTN:
    """IgBLASTN

    IgBLASTN class.  A tool for immunoglobulin (IG) and T cell receptor (TR) V domain sequences from nucletodies.

    This is a lower level class and probably should use airr to interact

    Examples
    --------
    >>> ig_blast = igblast.IgBLASTN()
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
        "_show_translation",
        "_extend_5",
        "_extend_3",
        "_j_penalty",
        "_organism",
        "_germline_db_v",
        "_germline_db_d",
        "_germline_db_j",
        "_aux_path",
        "_igdata",
        "_temp_dir",
    ]

    def __init__(self):
        """IgBLASTN with a query. Set everything up with a setter"""

        ###setup all the default values
        self.executable = ""
        self.min_d_match = 5
        self.num_v = 3
        self.num_d = 3
        self.num_j = 3
        self.outfmt = 19
        self.receptor = "Ig"
        self.word_size = 11
        self.nomenclature = "imgt"
        self.gap_open = 5
        self.gap_extend = 2
        self.num_threads = os.cpu_count()
        self.show_translation = True
        self.extend_5 = True
        self.extend_3 = True
        self.j_penalty = -1

        ##Make these blank, if they are not set by the caller, then we will complain during runtime.
        self._organism = IgBLASTArgument("organism", "organism", "", True)
        self._germline_db_v = IgBLASTArgument("germline_db_v", "germline_db_V", "", True)
        self._germline_db_d = IgBLASTArgument("germline_db_d", "germline_db_D", "", True)
        self._germline_db_j = IgBLASTArgument("germline_db_j", "germline_db_J", "", True)
        self._aux_path = IgBLASTArgument("aux_path", "auxiliary_data", "", True)

        ##Igdata is not an official blast argument, it is an enviroment
        self._igdata = ""
        self.temp_dir = "."

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
    def executable(self, path: Path):
        _executable = "igblastn"
        if not path:  # try and use package hmmscan
            system = platform.system().lower()
            igblastn_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"bin/{system}/{_executable}")
            # check if its
            if os.path.exists(igblastn_path):
                # check if it's executable
                if shutil.which(igblastn_path):
                    self._executable = shutil.which(igblastn_path)
                else:
                    # If it's not check the access
                    _access = os.access(igblastn_path, os.X_OK)
                    raise BadIgBLASTExe(igblastn_path, f"is, this executable? Executable-{_access}")
            else:  # The package igblastn is not working
                logger.warning(
                    f"Can't find igblast executable in {igblastn_path}, with system {system} within package {__package__}. Trying to find system installed hmmer"
                )
                igblastn_path = shutil.which(_executable)
                if igblastn_path:
                    self._executable = igblastn_path
                else:
                    raise BadIgBLASTExe(
                        igblastn_path, f"Can't find igblastn in package {__package__} or in path {os.env['PATH']}"
                    )
        else:  # User specifed custome path
            logger.debug(f"User passed custom igblastn {path}")
            igblastn_path = shutil.which(path)
            if igblastn_path:
                self._executable = igblastn_path
            else:
                _access = os.access(igblastn_path, os.X_OK)
                raise BadIgBLASTExe(
                    igblastn_path, f"Custom igblastn path is not executable {igblastn_path}, {_access} "
                )

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
    def temp_dir(self, data: Union[Path, str]):
        if isinstance(data, str):
            data = Path(data)
        self._temp_dir = data.absolute()

    @property
    def igdata(self) -> Path:
        """The path to IGDATA which contains the internal_data needed to make a recomdination

        Returns
        -------
        Path
           A valid IGDATA path
        """
        return self._igdata

    @igdata.setter
    def igdata(self, data: Path):
        if isinstance(data, str):
            data = Path(data)
        if not data.exists() or not data.is_dir():
            raise BadIgDATA(data)
        self._igdata = data.absolute()

    @property
    def min_d_match(self) -> IgBLASTArgument:
        """Required minimal consecutive nucleotide base matches for D genes

        Returns
        -------
        IgBLASTArgument
        """
        return self._min_d_match

    @min_d_match.setter
    def min_d_match(self, d: int):
        if not isinstance(d, int) and d < 5:
            raise BadIgBLASTArgument(d, ">5")
        self._min_d_match = IgBLASTArgument("min_d_match", "min_D_match", d, False)

    @property
    def num_v(self) -> IgBLASTArgument:
        """
           Number of Germline sequences to show alignments for

        Returns
        -------
        IgBLASTArgument
        """
        return self._num_v

    @num_v.setter
    def num_v(self, v):
        if not isinstance(v, int):
            raise BadIgBLASTArgument(v, int)
        self._num_v = IgBLASTArgument("num_v", "num_alignments_V", v, False)

    @property
    def num_d(self) -> IgBLASTArgument:
        """
           Number of Germline sequences to show alignments for D gene

        Returns
        -------
        IgBLASTArgument
        """
        return self._num_d

    @num_d.setter
    def num_d(self, d: int):
        if not isinstance(d, int):
            raise BadIgBLASTArgument(d, int)
        self._num_d = IgBLASTArgument("num_d", "num_alignments_D", d, False)

    @property
    def num_j(self) -> IgBLASTArgument:
        """
           Number of Germline sequences to show alignments for J gene

        Returns
        -------
        IgBLASTArgument
        """
        return self._num_j

    @num_j.setter
    def num_j(self, j: int):
        if not isinstance(j, int):
            raise BadIgBLASTArgument(j, int)
        self._num_j = IgBLASTArgument("num_j", "num_alignments_J", j, False)

    @property
    def organism(self) -> IgBLASTArgument:
        """The organism for your query sequence.

        Returns
        -------
        IgBLASTArgument
        """
        return self._organism

    @organism.setter
    def organism(self, o: str):
        # I don't want to hardcode in the organisms here.
        # I will handle that logic at a higher level,
        # this is because blast has no preset organims and it's all about the v,d,j blast paths which are set dynamically
        if not isinstance(o, str):
            raise BadIgBLASTArgument(o, str)
        self._organism = IgBLASTArgument("organism", "organism", o, True)

    @property
    def outfmt(self) -> IgBLASTArgument:
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
    def outfmt(self, fmt: int):
        ##only accept 19 for now
        if fmt != 19:
            raise BadIgBLASTArgument(fmt, 19)
        self._outfmt = IgBLASTArgument("outfmt", "outfmt", fmt, True)

    @property
    def receptor(self) -> IgBLASTArgument:
        """
        Specify Ig or T cell receptor sequence

        Returns
        -------
        IgBLASTArgument
        """
        return self._receptor

    @receptor.setter
    def receptor(self, r):
        if not isinstance(r, str) and not (r in ["Ig", "TCR"]):
            raise BadIgBLASTArgument(r, ["Ig", "TCR"])
        self._receptor = IgBLASTArgument("receptor", "ig_seqtype", r, True)

    @property
    def nomenclature(self) -> IgBLASTArgument:
        """Domain system to be used for segment annotation

        Returns
        -------
        IgBLASTArgument
        """
        return self._nomenclature

    @nomenclature.setter
    def nomenclature(self, system: str):
        if system.lower() not in ["imgt", "kabat"]:
            raise BadIgBLASTArgument(system, "['imgt','kaba']")
        self._nomenclature = IgBLASTArgument("nomenclature", "domain_system", system, True)

    @property
    def aux_path(self) -> IgBLASTArgument:
        """Auxilary data path. This is needed to lookup the J genes and tell them when the CDR3 stops.

        Returns
        -------
        IgBLASTArgument
        """
        return self._aux_path

    @aux_path.setter
    def aux_path(self, aux_path: Path):
        if isinstance(aux_path, str):
            aux_path = Path(aux_path)
        if not aux_path.exists():
            raise BadIgBLASTArgument(aux_path, "valid path to Auxilary database")
        self._aux_path = IgBLASTArgument("aux_path", "auxiliary_data", aux_path.absolute(), True)

    @property
    def germline_db_v(self) -> IgBLASTArgument:
        """Path to V gene database prefix

        Returns
        -------
        IgBLASTArgument
        """
        return self._germline_db_v

    @germline_db_v.setter
    def germline_db_v(self, path: str):
        abs_path = ensure_prefix_to(path)
        if not abs_path:
            raise BadIgBLASTArgument(path, "Valid path to V Database")
        self._germline_db_v = IgBLASTArgument("germline_db_v", "germline_db_V", path, True)

    @property
    def germline_db_d(self) -> IgBLASTArgument:
        """Path to D gene database prefix

        Returns
        -------
        IgBLASTArgument
        """
        return self._germline_db_d

    @germline_db_d.setter
    def germline_db_d(self, path: str):
        abs_path = ensure_prefix_to(path)
        if not abs_path:
            raise BadIgBLASTArgument(path, "Valid path to D Database")
        self._germline_db_d = IgBLASTArgument("germline_db_d", "germline_db_D", path, True)

    @property
    def germline_db_j(self: str) -> IgBLASTArgument:
        """Path to J gene database prefix

        Returns
        -------
        IgBLASTArgument
        """
        return self._germline_db_j

    @germline_db_j.setter
    def germline_db_j(self, path):

        abs_path = ensure_prefix_to(path)
        if not abs_path:
            raise BadIgBLASTArgument(path, "Valid path to J Database")
        self._germline_db_j = IgBLASTArgument("germline_db_j", "germline_db_J", path, True)

    @property
    def word_size(self) -> IgBLASTArgument:
        """Word size for wordfinder algorithm (length of best perfect match)

        Returns
        -------
        IgBLASTArugment
        """
        return self._word_size

    @word_size.setter
    def word_size(self, word_size: int):
        if not isinstance(word_size, int) and word_size < 4:
            raise BadIgBLASTArgument(word_size, ">4")
        self._word_size = IgBLASTArgument("word_size", "word_size", word_size, False)

    @property
    def gap_open(self) -> IgBLASTArgument:
        """Cost to open a gap

        Returns
        -------
        IgBLASTArgument
        """
        return self._gap_open

    @gap_open.setter
    def gap_open(self, go: int):
        if not isinstance(go, int) and go > 0:
            raise BadIgBLASTArgument(go, ">0")
        self._gap_open = IgBLASTArgument("gap_open", "gapopen", go, False)

    @property
    def gap_extend(self) -> IgBLASTArgument:
        """Cost to extend a gap

        Returns
        -------
        IgBLASTArgument
        """
        return self._gap_extend

    @gap_extend.setter
    def gap_extend(self, ge: int):
        if not isinstance(ge, int) and ge > 0:
            raise BadIgBLASTArgument(ge, ">0")
        self._gap_extend = IgBLASTArgument("gap_open", "gapextend", ge, False)

    @property
    def num_threads(self) -> IgBLASTArgument:
        """
            Number of threads (CPUs) to use in the BLAST search

        Returns
        -------
        IgBLASTArgument
        """
        return self._num_threads

    @num_threads.setter
    def num_threads(self, num_threads: int):
        if num_threads > os.cpu_count():
            raise BadIgBLASTArgument(num_threads, "<" + str(os.cpu_count()))
        self._num_threads = IgBLASTArgument("number_threads", "num_threads", num_threads, False)

    @property
    def show_translation(self) -> IgBLASTArgument:
        """show translated_alignments

        Returns
        -------
        IgBLASTArgument
        """
        return self._show_translation

    @show_translation.setter
    def show_translation(self, b: bool):
        self._show_translation = IgBLASTArgument("show_translation", "show_translation", b, False)

    @property
    def extend_5(self) -> IgBLASTArgument:
        """Extend V gene alignment at 5' end

        Returns
        -------
        IgBLASTArgument
        """
        return self._extend_5

    @extend_5.setter
    def extend_5(self, extend_5: bool):
        self._extend_5 = IgBLASTArgument("extend_5", "extend_align5end", extend_5, False)

    @property
    def extend_3(self) -> IgBLASTArgument:
        """Extend V gene alignment at 3' end

        Returns
        -------
        IgBLASTArgument
        """
        return self._extend_3

    @extend_3.setter
    def extend_3(self, extend_3: bool):
        self._extend_3 = IgBLASTArgument("extend_3", "extend_align3end", True, False)

    @property
    def j_penalty(self) -> IgBLASTArgument:
        """What is the  J gene panalty

        Returns
        -------
        IgBLASTArgument
        """
        return self._j_penalty

    @j_penalty.setter
    def j_penalty(self, penalty: int):
        if not -5 < penalty < 1:
            raise BadIgBLASTArgument(penalty, "must be less than 0 and greater than -5")
        self._j_penalty = IgBLASTArgument("j_penalty", "J_penalty", penalty, True)

    @property
    def arguments(self) -> List[IgBLASTArgument]:
        """return a list of IgBLASTArugments

        Returns
        -------
        List[IgBLASTArguments]
        """
        return [
            self.min_d_match,
            self.num_v,
            self.num_j,
            self.num_d,
            self.organism,
            self.receptor,
            self.germline_db_v,
            self.germline_db_d,
            self.germline_db_j,
            self.aux_path,
            self.outfmt,
            self.nomenclature,
            self.word_size,
            self.gap_open,
            self.gap_extend,
            self.j_penalty,
            self.num_threads,
            self.show_translation,
            self.extend_5,
            self.extend_3,
        ]

    @property
    def cmd(self) -> list:
        """Return the blast cmd that will be run by subprocess"""
        _cmd = [self.executable]
        for blast_arg in self.arguments:
            kv = blast_arg.get_formatted_blast_arg()
            if kv:
                _cmd += kv
        return _cmd

    def _pre_check(self):
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
        # Ensure required arguments werer set
        for blast_arg in self.arguments:
            if blast_arg.required and not (blast_arg.value):
                raise MissingIgBLASTArg(f"Missing Blast argument. Need to set IgBLASTN.{blast_arg.name}")

            # # Check the executable
            if not is_tool(self.executable):
                raise BadIgBLASTExe(self.executable)

            if not self.igdata:
                raise BadIgDATA("No IGDATA set, set with IgBLASTN.igdata")

            else:
                if not os.path.exists(self.igdata):
                    raise BadIgDATA(self.igdata)

    def run_file(self, file: Path) -> pd.DataFrame:
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
        EmtpyFileError
           if the fasta file is empty
        IgBLASTRunTimeError
           for any given runtime error for igblastN
        """
        local_env = os.environ.copy()
        if os.path.getsize(file) == 0:
            raise EmtpyFileError(file)
        local_env["IGDATA"] = self.igdata
        cmd = self.cmd
        cmd += ["-query", file]
        self._pre_check()
        # while we can certainly do this as an output stream on stdout,
        # It's probably best to take advantage of IGblast output
        with tempfile.NamedTemporaryFile(dir=self.temp_dir, suffix="_igblast.tsv") as tmpfile:
            cmd += ["-out", tmpfile.name]
            if sys.version_info.minor == 6:
                process = subprocess.run(cmd, env=local_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if process.stderr:
                    raise IgBLASTRunTimeError(process.stderr)
            else:
                process = subprocess.run(cmd, env=local_env, capture_output=True)
                if process.stderr:
                    raise IgBLASTRunTimeError(process.stderr)
            return pd.read_csv(tmpfile.name, sep="\t")

    def run_single(self, q: str) -> pd.DataFrame:
        """Run Igblast on a single fasta string

        Parameters
        ----------
        q : str
            A string fasta, ex ">my_file\nATCACA..."

        Returns
        -------
        pd.DataFrame
            A dataframe with the IgBLAST results

        Raises
        ------
        BadIgBLASTArgument
            During precheck if any blast arguments are incorrect
        IgBLASTRunTimeError
            Any runtime error for igblast
        """
        if not isinstance(q, str):
            raise BadIgBLASTArgument(type(q), "needs to be instance str")
        if not q:
            raise BadIgBLASTArgument(q, "Input query is null, please provide sequence")
        # Local Env
        local_env = os.environ.copy()
        local_env["IGDATA"] = self.igdata
        self._pre_check()

        # run process
        if sys.version_info.minor == 6:
            process = subprocess.run(
                self.cmd, env=local_env, input=q.encode("utf-8"), stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        else:
            process = subprocess.run(self.cmd, env=local_env, input=q.encode("utf-8"), capture_output=True)
        stdout = process.stdout.decode("utf-8")
        string_io = StringIO(stdout)
        stderr = process.stderr
        if process.stderr:
            raise IgBLASTRunTimeError(stderr)
        else:
            return pd.read_csv(string_io, sep="\t")

    def __repr__(self):
        return "IgBLAST: env IGDATA={} {}".format(self.igdata, " ".join(self.cmd))

    def __str__(self):
        return self.__repr__()


if __name__ == "__main__":
    ig_blast = IgBLASTN()