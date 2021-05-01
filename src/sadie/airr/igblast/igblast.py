"""Low level module for IgBLAST api calls

You probably only want to interact with this module through airr as the database files are extremely tricky to get right

"""
import glob
import logging
import os
import subprocess
import tempfile
import warnings
import semantic_version
from shutil import which

from pathlib import Path
from typing import List, Union

# Third party
import pandas as pd

# package/module level
from sadie.utility.util import is_tool
from sadie.airr.airrtable.constants import IGBLAST_AIRR
from sadie.airr.exceptions import (
    BadIgBLASTArgument,
    BadIgBLASTExe,
    BadIgDATA,
    MissingIgBLASTArgument,
    EmtpyFileError,
    IgBLASTRunTimeError,
)


# get logger in global scope
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
        "_show_translation",
        "_extend_5",
        "_extend_3",
        "_j_penalty",
        "_v_penalty",
        "_d_penalty",
        "_organism",
        "_germline_db_v",
        "_germline_db_d",
        "_germline_db_j",
        "_aux_path",
        "_igdata",
        "_temp_dir",
        "_allow_vdj_overlap",
    ]

    def __init__(self, executable="igblastn"):
        """IgBLASTN with a query. Set everything up with a setter"""

        # set the executable dynamically
        self.executable = executable
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
        self.num_threads = os.cpu_count()
        self.show_translation = True
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
        self._aux_path = IgBLASTArgument("aux_path", "auxiliary_data", "", True)

        # Igdata is not an official blast argument, it is an enviroment
        self._igdata = ""
        self.temp_dir = "."

    def _get_version(self) -> semantic_version.Version:
        """Private method to parse igblast -version and get semantic_version

        Returns
        -------
        semantic_version.Version
            the igblast version
        """
        process = subprocess.run([self.executable, "-version"], capture_output=True)
        stdout = process.stdout.decode("utf-8")
        if process.stderr:
            logger.error(
                f"{self.executable}, has no returned and error when checking version,. Tried igblastn -version: {process.stderr.decode('utf-8')}"
            )
            raise BadIgBLASTExe(self.executable, process.stderr.decode("utf-8"))
        version = stdout.split("\n")[0].split(":")[-1].strip()
        try:
            version = semantic_version.Version(version)
        except ValueError:
            raise BadIgBLASTExe(self.executable, f"semantic version can't parse {version}")
        return version

    @property
    def executable(self) -> Path:
        """The igblastn executable

        Returns
        -------
        Path
            igblastn path
        """
        return self._executable

    @property
    def version(self) -> semantic_version.Version:
        return self._version

    @executable.setter
    def executable(self, exe: Path):
        if isinstance(exe, str):
            exe = Path(exe)
        full_exe_path = which(exe)
        if not full_exe_path:
            raise BadIgBLASTExe(exe, f"{exe} must exist")
        if not os.access(full_exe_path, os.X_OK):
            raise BadIgBLASTExe(exe, f"{full_exe_path} must be executable")
        self._executable = Path(exe)

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
        """Organism

        Parameters
        ----------
        o : str
            an organism string

        Raises
        ------
        BadIgBLASTArgument
            if igblast is not a str
        """
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
        # only accept 19 for now
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
            warnings.warn(f"{path} is not found, No D gene segment", UserWarning)
            # raise BadIgBLASTArgument(path, "Valid path to D Database")
            self._germline_db_d = IgBLASTArgument("germline_db_d", "germline_db_D", "", False)
        else:
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
        self._extend_3 = IgBLASTArgument("extend_3", "extend_align3end", extend_3, False)

    @property
    def allow_vdj_overlap(self) -> IgBLASTArgument:
        """Allow the VDJ overlap
        This option is active only when D_penalty
        and J_penalty are set to -4 and -3, respectively
        Returns
        -------
        IgBLASTArgument
        """
        return self._allow_vdj_overlap

    @allow_vdj_overlap.setter
    def allow_vdj_overlap(self, allow: bool):
        if self.j_penalty.value != -3 and self.d_penalty.value != -4 and allow:
            warnings.warn(
                f"Allows vdj overlap set but j penalty and d penalty need to be -3 and -4, now are {self.j_penalty}, {self.d_penalty}",
                UserWarning,
            )
        self._allow_vdj_overlap = IgBLASTArgument("allow_vdj_overlap", "allow_vdj_overlap", allow, False)

    @property
    def d_penalty(self) -> IgBLASTArgument:
        """What is the  D gene panalty

        Returns
        -------
        IgBLASTArgument
        """
        return self._d_penalty

    @d_penalty.setter
    def d_penalty(self, penalty: int):
        if not -5 < penalty < 1:
            raise BadIgBLASTArgument(penalty, "must be less than 0 and greater than -5")
        self._d_penalty = IgBLASTArgument("d_penalty", "D_penalty", penalty, True)

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
        if not -4 < penalty < 1:
            raise BadIgBLASTArgument(penalty, "must be less than 0 and greater than -4")
        self._j_penalty = IgBLASTArgument("j_penalty", "J_penalty", penalty, True)

    @property
    def v_penalty(self) -> IgBLASTArgument:
        """What is the  v gene panalty

        Returns
        -------
        IgBLASTArgument
        """
        return self._v_penalty

    @v_penalty.setter
    def v_penalty(self, penalty: int):
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
            self.v_penalty,
            self.d_penalty,
            self.num_threads,
            self.show_translation,
            self.extend_5,
            self.extend_3,
            self.allow_vdj_overlap,
        ]

    @property
    def cmd(self) -> list:
        """Return the blast cmd that will be run by subprocess"""
        _cmd = [str(self.executable)]
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
                raise MissingIgBLASTArgument(f"Missing Blast argument. Need to set IgBLASTN.{blast_arg.name}")

            #  Check the executable
            if not is_tool(self.executable):
                raise BadIgBLASTExe(self.executable, "Is not an executable tool")

            if not self.igdata:
                raise BadIgDATA("No IGDATA set, set with IgBLASTN.igdata")

            else:
                if not os.path.exists(self.igdata):
                    raise BadIgDATA(self.igdata)

    # Run methods
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
           for any given runtime error for igblastn
        """

        # because igblast uses IGDATA as the internal file structure, we should pass the enviroment to the subprocess
        local_env = os.environ.copy()
        local_env["IGDATA"] = str(self.igdata)

        # we want to ensure they actually passed a file with stuff in it
        if os.path.getsize(file) == 0:
            raise EmtpyFileError(file)

        # take the cmd and finally add the query file
        cmd = self.cmd
        cmd += ["-query", file]

        # run a precheck to make sure everything passed was working
        self._pre_check()

        # while we can certainly do this as an output stream on stdout,
        # It's probably best to take advantage of IGblast output and tempfile
        with tempfile.NamedTemporaryFile(dir=self.temp_dir, suffix="_igblast.tsv") as tmpfile:
            cmd += ["-out", tmpfile.name]

            process = subprocess.run(cmd, env=local_env, capture_output=True)
            if process.stderr:
                raise IgBLASTRunTimeError(process.stderr)
            # we read the dataframe from the tempfile, it should always be in .TSV.
            # We can also cast it to IGBLAST_AIRR dtypes to save memory
            df = pd.read_csv(tmpfile.name, sep="\t", dtype=IGBLAST_AIRR)
        if Path(tmpfile.name).exists():
            logger.debug(f"{tmpfile.name} was not deleted after it exited scope")
            Path(tmpfile.name).unlink()

        return df

    def __repr__(self):
        return "IgBLAST: env IGDATA={} {}".format(str(self.igdata), " ".join(self.cmd))

    def __str__(self):
        return self.__repr__()


if __name__ == "__main__":
    ig_blast = IgBLASTN()
