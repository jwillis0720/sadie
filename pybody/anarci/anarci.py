"""The ANARCI Abstraction to make it usuable"""

import bz2
import glob
import gzip
import json
import logging
import multiprocessing
import os
import pickle
import shutil
import tempfile
import platform

# Std library
import warnings
from multiprocessing import cpu_count
from pathlib import Path
from typing import List, Tuple, Union

import filetype
import pandas as pd

# third party
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .aa._anarci import check_for_j, dataframe_output, number_sequences_from_alignment, run_hmmer
from .numbering import scheme_numbering
from .exception import AnarciDuplicateIdError, BadAnarciArgument, BadRequstedFileType, HmmerExecutionError
from .result import AnarciResult, AnarciResults
from ..utility import split_fasta

logger = logging.getLogger("ANARCI")

# Get out of here with your partial codon warnigns
warnings.filterwarnings("ignore", "Partial codon")


class Anarci:
    def __init__(
        self,
        scheme="imgt",
        region_assign="imgt",
        allowed_chain=["H", "K", "L"],
        assign_germline=True,
        allowed_species=["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca", "dog", "cat"],
        tempdir="",
        hmmerpath="",
    ):
        """Anarci constructor"""
        self.scheme = scheme
        self.region_definition = region_assign
        self.allowed_chains = allowed_chain
        self.assign_germline = assign_germline
        self.allowed_species = allowed_species
        self.tempdir = tempdir
        self.hmmerpath = hmmerpath
        self.num_cpus = cpu_count()

    @property
    def hmmerpath(self) -> Path:
        """Path to the hmmrscan executable

        Returns
        -------
        Path
            The path of the hmmrscan
        """
        return self._hmmerpath

    @hmmerpath.setter
    def hmmerpath(self, path: Path):
        _executable = "hmmscan"
        if not path:  # try and use package hmmscan
            system = platform.system().lower()
            hmmer_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"bin/{system}/{_executable}")
            # check if its
            if os.path.exists(hmmer_path):
                # check if it's executable
                if shutil.which(hmmer_path):
                    self._hmmerpath = hmmer_path
                else:
                    # If it's not check the access
                    _access = os.access(hmmer_path, os.X_OK)
                    raise HmmerExecutionError(hmmer_path, f"is, this executable? Executable-{_access}")
            else:  # The package hmmscan is not working
                logger.warning(
                    f"Can't find hmmscan executable in {hmmer_path}, with system {system} within package {__package__}. Trying to find system installed hmmer"
                )
                hmmer_path = shutil.which(_executable)
                if hmmer_path:
                    self._hmmerpath = hmmer_path
                else:
                    raise HmmerExecutionError(
                        hmmer_path, f"Can't find hmmrscan in package {__package__} or in path {os.env['PATH']}"
                    )
        else:  # User specifed custome path
            logger.debug(f"User passed custom hmmer path {path}")
            hmmer_path = shutil.which(path)
            if hmmer_path:
                self._hmmerpath = hmmer_path
            else:
                _access = os.access(hmmer_path, os.X_OK)
                raise HmmerExecutionError(hmmer_path, f"Custom hmmer path is not executable {hmmer_path}, {_access} ")

    @property
    def region_definition(self) -> str:
        """Region defiition, eg. imgt, chotia"""
        return self._region_definition

    @region_definition.setter
    def region_definition(self, definition: str):
        """The region defitinions that should be applied scheme that should be applied

        accepted: imgt, kabat, chotia, martin, abm

        """
        _accepted_defs = ["imgt", "kabat", "chothia", "abm", "contact", "scdr"]
        if definition.lower() not in _accepted_defs:
            raise BadAnarciArgument(definition, _accepted_defs)
        self._region_definition = definition

    @property
    def scheme(self) -> str:
        """The numbering scheme that should be applied"""
        return self._scheme

    @scheme.setter
    def scheme(self, scheme: str):
        """The numbering scheme that should be applied

        accepted: imgt, kabat, chotia, martin

        """
        _accepted_schemenes = ["imgt", "kabat", "chothia"]  #
        __future_schemes = ["martin"]
        if scheme.lower() not in _accepted_schemenes:
            print(f"need support for {__future_schemes} numbering schemes. See abysis")
            raise BadAnarciArgument(scheme, _accepted_schemenes)
        self._scheme = scheme.lower()

    @property
    def tempdir(self) -> Path:
        return self._path

    @tempdir.setter
    def tempdir(self, path: str):
        if not path:
            if "TMPFILE" in os.environ:
                _path = os.environ["TMPFILE"]
            else:
                _path = "/tmp"
        else:
            _path = path
        _path = Path(_path)
        # check path exists
        if _path.exists() and _path.is_dir():
            self._path = _path
        elif not _path.exists():
            _path.mkdir()
            self._path = _path
        else:
            raise TypeError(f"{path} is not found or can't be created")

    @property
    def allowed_chains(self) -> List[str]:
        """The chain types to consider in the alignment,

        Returns
        -------
        List[str]
            [description]
        """
        return self._allowed_chains

    @allowed_chains.setter
    def allowed_chains(self, allowed_chains: List[str]):
        """A list of single character chains

        H - Heavy
        K - Kappa
        L - Lambda
        A - Alpha T cell receptor
        B - Beta T cell receptor
        G - Gamma T cell receptor
        D - Delta T cell receptor


        Parameters
        ----------
        allowed_chains : list,
            e.g, ['H','K'] will only search heavy and kappa chains
        """
        _allowed_chain = ["H", "K", "L", "A", "B", "G", "D"]
        _diff = list(set(map(lambda x: x.upper(), allowed_chains)).difference(_allowed_chain))
        if _diff:
            raise BadAnarciArgument(_diff, _allowed_chain)
        self._allowed_chains = allowed_chains

    @property
    def allowed_species(self) -> list:
        return self._allowed_species

    @allowed_species.setter
    def allowed_species(self, allowed_species: List[str]):
        """If assign_germline is true, limit the species that can be assigned to a limited set.
        Useful when the animal species is known or when performing closest germline experiments


        Parameters
        ----------
        allowed_species: list,
            ["human", "mouse", "rat", "rabbit", "rhesus ", "pig", "alpaca"],
        """
        _allowed_species = ["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca", "dog", "cat"]
        _diff = list(set(map(lambda x: x.lower(), allowed_species)).difference(_allowed_species))
        if _diff:
            raise BadAnarciArgument(_diff, _allowed_species)
        self._allowed_species = allowed_species

    @property
    def assign_germline(self) -> bool:
        """Should Anarci try to assign germline

        Returns
        -------
        bool
        """
        return self._assign_germline

    @assign_germline.setter
    def assign_germline(self, assign: bool):
        """Should anarci try to assign germline

        Parameters
        ----------
        assign : bool

        Raises
        ------
        BadAnarciArgument
            if not a bool
        """
        if not isinstance(assign, bool):
            raise BadAnarciArgument(assign, bool)
        self._assign_germline = assign

    def _run(self, sequences: List[Tuple]):
        """
        private method to run Anarci

        Parameters
        ----------
        sequences : List[Tuple]
            list or tuple of (Id, Sequence) pairs
                              e.g. [ ("seq1","EVQLQQSGAEVVRSG ..."),
                                     ("seq2","DIVMTQSQKFMSTSV ...")
        """

        # Perform the alignments of the sequences to the hmm database
        _alignments = run_hmmer(
            sequences,
            hmm_database="ALL",
            hmmerpath=self.hmmerpath,
            ncpu=self.num_cpus,
            bit_score_threshold=80,
            tempdir=self.tempdir,
        )

        # Check the numbering for likely very long CDR3s that will have been missed by the first pass.
        # Modify alignments in-place
        check_for_j(sequences, _alignments, self.scheme)

        # # Apply the desired numbering scheme to all sequences
        _numbered, _alignment_details, _hit_tables = number_sequences_from_alignment(
            sequences,
            _alignments,
            scheme=self.scheme,
            allow=self.allowed_chains,
            assign_germline=self.assign_germline,
            allowed_species=self.allowed_species,
        )

        _summary, _alignment = dataframe_output(sequences, _numbered, _alignment_details)
        _summary["scheme"] = self.scheme
        _alignment["scheme"] = self.scheme
        _summary["allowed_species"] = ",".join(self.allowed_species)
        _summary["allowed_chains"] = ",".join(self.allowed_chains)
        _summary["region_definition"] = self.region_definition

        _results = []
        if _summary.empty:
            # If we have nothing
            return [None]
        for group_id, summary_df in _summary.groupby("Id"):
            if len(summary_df) > 1:
                dup_ids = list(summary_df["Id"].unique())
                found_number = len(summary_df)
                raise AnarciDuplicateIdError(dup_ids, found_number)
            alignment_df = _alignment.loc[_alignment["Id"] == group_id]
            summary_series = summary_df.iloc[0]
            _results.append(AnarciResult(summary_series, alignment_df))
        return _results

    def run_single(self, seq_id: str, seq: str, scfv=False) -> "AnarciResult":
        """Run a single string sequence on an amino acid

        Parameters
        ----------
        seq_id : str
           the sequence_id of the string object, ex. "my_sequence"
        seq : str
            The string nucletodide sequence, ex. "EVQLQQSGAEVVRSG ..."

        Returns
        -------
            AnarchiResult Object
        """

        if not scfv:
            sequences = [(seq_id, seq)]

        return self._run(sequences)[0]

    def run_multiple(self, seqrecords: List[SeqRecord], scfv=False) -> "AnarciResults":
        """Run multiple seq records

        Parameters
        ----------
        seqrecords : List[SeqRecord]
            A list of sequence records of amino acids.

        Returns
        -------
            ANARCIResults - Holds many results

        Raises
        ------
        TypeError
            if you don't pass a list of SeqRecords
        """
        if not isinstance(seqrecords, (list, type(SeqIO.FastaIO.FastaIterator))):
            raise TypeError(f"seqrecords must be of type {list} pased {type(seqrecords)}")

        if isinstance(seqrecords, list) and not all([type(i) is SeqRecord for i in seqrecords]):
            raise TypeError("seqrecords argument must be of a list of Seqrecords")

        _sequences = []
        _seen = set()
        for seq in seqrecords:
            if seq.id in _seen:
                raise AnarciDuplicateIdError(seq.id, 1)
            _sequences.append((seq.id, str(seq.seq)))
            _seen.add(seq.id)

        _results = self._run(_sequences)
        return AnarciResults(results=_results)

    def run_file(self, file: Path, multi=False) -> "AnarciResults":
        """Run anarci annotator on a fasta file

        Parameters
        ----------
        file: Path
            The fasta file to run
        multi: Bool, defaults=True
            split and run file as multiprocess

        Returns
        -------
        AnarciResults
            Returns AnarciResults object

        Raises
        ------
        FileExistsError
            if file does not exist
        BadRequstedFileType
            if file is not fasta

        """
        _filetype = filetype.guess(file)
        if not os.path.exists(file):
            raise FileExistsError(f"{file} not found")
        if _filetype:
            logger.info("Guess File Type is %s ", _filetype.extension)
            with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
                if _filetype.extension == "gz":
                    logger.info("File type is compressed gzip")
                    file_buffer = gzip.open(file)
                elif _filetype.extension == "bz2":
                    logger.info("File type is compressed bzip2")
                    file_buffer = bz2.open(file)
                else:
                    raise BadRequstedFileType(_filetype, ["bzip2", "gzip"])
                shutil.copyfileobj(file_buffer, tmpfile)
                file = tmpfile.name

        if not multi:
            # run on fasta
            _results = self.run_multiple(list(SeqIO.parse(file, "fasta")))
            # if we had a file delete it
            if _filetype:
                os.unlink(tmpfile.name)
        else:
            _results = self._run_mp(file)
        return _results

    def _run_mp(self, table: Path) -> "AnarciResults":
        results = []
        if isinstance(table, (str, Path)):
            with tempfile.TemporaryDirectory(dir=".") as d:
                split_fasta(table, 1000, d)
                files = list(glob.glob(d + "/*.fasta*"))
                p = multiprocessing.Pool()
                results = AnarciResults.concat(p.map(self.run_file, files))
        return results


if __name__ == "__main__":

    from pprint import pprint

    anarci_api = Anarci(scheme="chothia", region_assign="scdr")
    result = anarci_api.run_single(
        "MySweetAntibody",
        "EVQLLESGGGLVQPGGSLRLSCAASGFTFPVYNMAWVRQAPGKGLEWVSGIAHNGRNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKGHEISRFSRWSSFDYWGQGTLVTVSS",
    )

    pprint(result.segment_table_no_gaps)
    anarci_api = Anarci(scheme="kabat", region_assign="scdr")
    result = anarci_api.run_single(
        "MySweetAntibody",
        "DIQMTQSPSSLSASVGDRVTITCRPNQNIATYINWYQQKPGKAPKLLIYAASGLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQHSWEIPYTFGQGTKVEIK",
    )
    pprint(result.segment_table_no_gaps)
