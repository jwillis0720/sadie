"""The ANARCI Abstraction to make it usuable"""

# Std library
import bz2
import gzip
import logging
import multiprocessing
import filetype
import os
import shutil
import tempfile
import platform

import warnings
from multiprocessing import cpu_count
from pathlib import Path
from typing import List, Tuple, Union, Generator

# third party
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

from .aa._anarci import (
    check_for_j,
    parsed_output,
    number_sequences_from_alignment,
    run_hmmer,
)
from .exception import (
    AnarciDuplicateIdError,
    BadAnarciArgument,
    BadRequstedFileType,
    HmmerExecutionError,
)
from .result import AnarciResults
from .constants import ANARCI_RESULTS
from .numbering import scheme_numbering

logger = logging.getLogger("ANARCI")

# Get out of here with your partial codon warnigns
warnings.filterwarnings("ignore", "Partial codon")


class Error(Exception):
    """Base class for exceptions in this module."""


class Anarci:
    def __init__(
        self,
        scheme="imgt",
        region_assign="imgt",
        allowed_chain=["H", "K", "L"],
        assign_germline=True,
        allowed_species=[
            "human",
            "mouse",
            "rat",
            "rabbit",
            "rhesus",
            "pig",
            "alpaca",
            "dog",
            "cat",
        ],
        tempdir="",
        hmmerpath="",
        threshold=80,
        run_multiproc=True,
    ):
        """[summary]

        Parameters
        ----------
        scheme : str, optional
            [description], by default "imgt"
        region_assign : str, optional
            [description], by default "imgt"
        allowed_chain : list, optional
            [description], by default ["H", "K", "L"]
        assign_germline : bool, optional
            [description], by default True
        allowed_species : list, optional
            [description], by default [ "human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca", "dog", "cat", ]
        tempdir : str, optional
            [description], by default ""
        hmmerpath : str, optional
            [description], by default ""
        threshold : int, optional
            [description], by default 80
        run_multiproc : bool, optional
            [description], by default True

        Raises
        ------
        NotImplementedError
            [description]
        """
        self.scheme = scheme
        self.region_definition = region_assign
        self.allowed_chains = allowed_chain
        self.assign_germline = assign_germline
        self.allowed_species = allowed_species
        self.tempdir = tempdir
        self.hmmerpath = hmmerpath
        self.num_cpus = cpu_count()
        self.run_multiproc = run_multiproc
        self.threshold_bit = threshold
        if not self.check_combination(self.scheme, self.region_definition):
            raise NotImplementedError(f"{self.scheme} with {self.region_definition} has not been implemented yet")

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
            hmmer_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                f"bin/{system}/{_executable}",
            )
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
                        hmmer_path,
                        f"Can't find hmmrscan in package {__package__} or in path {os.environ['PATH']}",
                    )
        else:  # User specifed custome path
            logger.debug(f"User passed custom hmmer path {path}")
            hmmer_path = shutil.which(path)
            if hmmer_path:
                self._hmmerpath = hmmer_path
            else:
                _access = os.access(hmmer_path, os.X_OK)
                raise HmmerExecutionError(
                    hmmer_path,
                    f"Custom hmmer path is not executable {hmmer_path}, {_access} ",
                )

    @property
    def region_definition(self) -> str:
        """Region defiition, eg. imgt, chotia"""
        return self._region_definition

    @region_definition.setter
    def region_definition(self, definition: str):
        """The region defitinions that should be applied scheme that should be applied

        accepted: imgt, kabat, chotia, martin, abm

        """
        if definition.lower() not in self.get_available_region_definitions():
            raise BadAnarciArgument(definition, self.get_available_region_definitions())
        self._region_definition = definition

    @staticmethod
    def get_available_region_definitions() -> List:
        """Get currently available antibody region definitions

        Returns
        -------
        List
            a list of region defitions, ex. ["imgt", "kabat", "chothia", "abm", "contact", "scdr"]

        """
        _accepted_defs = ["imgt", "kabat", "chothia", "abm", "contact", "scdr"]
        return _accepted_defs

    @property
    def scheme(self) -> str:
        """The numbering scheme that should be applied"""
        return self._scheme

    @scheme.setter
    def scheme(self, scheme: str):
        """The numbering scheme that should be applied

        accepted: imgt, kabat, chotia, martin

        """
        __future_schemes = ["martin"]
        if scheme.lower() not in self.get_available_region_definitions():
            logger.warning(f"need support for {__future_schemes} numbering schemes. See abysis")
            raise BadAnarciArgument(scheme, self.get_available_region_definitions())
        self._scheme = scheme.lower()

    @staticmethod
    def get_available_numbering_schemes() -> List:
        """Get currently available antibody numbering schemes

        Returns
        -------
        List
            a list of region defitions, ex. ["imgt", "kabat", "chothia"]

        """
        _accepted_schemes = ["imgt", "kabat", "chothia"]  #
        return _accepted_schemes

    @staticmethod
    def check_combination(scheme: str, region: str) -> bool:
        scheme_keys = scheme_numbering[scheme]
        try:
            scheme_keys["heavy"][region]
            scheme_keys["light"][region]
        except KeyError:
            return False
        return True

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
        _allowed_chain = self.get_allowed_chains()
        _diff = list(set(map(lambda x: x.upper(), allowed_chains)).difference(_allowed_chain))
        if _diff:
            raise BadAnarciArgument(_diff, _allowed_chain)
        self._allowed_chains = allowed_chains

    @staticmethod
    def get_allowed_chains() -> List:
        """Get the allowed chains options. Which chains can you align against

        Returns
        -------
        List
            A list of one letter codes that correspond to chains
        """
        _allowed_chain = ["H", "K", "L", "A", "B", "G", "D"]
        return _allowed_chain

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
        _allowed_species = self.get_allowed_species()
        _diff = list(set(map(lambda x: x.lower(), allowed_species)).difference(_allowed_species))
        if _diff:
            raise BadAnarciArgument(_diff, _allowed_species)
        self._allowed_species = allowed_species

    @staticmethod
    def get_allowed_species() -> List:
        """Get allowed species that we should align against.

        Returns
        -------
        List
            A list of currently implmented allowed species
        """
        _allowed_species = [
            "human",
            "mouse",
            "rat",
            "rabbit",
            "rhesus",
            "pig",
            "alpaca",
            "dog",
            "cat",
        ]
        return _allowed_species

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
            bit_score_threshold=self.threshold_bit,
            tempdir=self.tempdir,
        )

        # Check the numbering for likely very long CDR3s that will have been missed by the first pass.
        # Modify alignments in-place
        check_for_j(sequences, _alignments, self.scheme, self.hmmerpath)

        # # Apply the desired numbering scheme to all sequences
        _numbered, _alignment_details, _hit_tables = number_sequences_from_alignment(
            sequences,
            _alignments,
            scheme=self.scheme,
            allow=self.allowed_chains,
            assign_germline=self.assign_germline,
            allowed_species=self.allowed_species,
        )

        _summary = parsed_output(sequences, _numbered, _alignment_details)
        anarci_results = pd.DataFrame(_summary)
        if anarci_results.empty:
            return AnarciResults()

        # I really want to set the scheme and region in the constructor
        # https://stackoverflow.com/questions/66647680/subclassing-pandas-dataframe-and-setting-field-in-constuctor
        anarci_results = AnarciResults(
            anarci_results.astype(ANARCI_RESULTS),
        )

        # Must set these schemes before we set the segments
        anarci_results["scheme"] = self.scheme
        anarci_results["region_definition"] = self.region_definition
        anarci_results["allowed_species"] = ",".join(self.allowed_species)
        anarci_results["allowed_chains"] = ",".join(self.allowed_chains)
        anarci_results = anarci_results._add_segment_regions()

        if len(anarci_results["Id"].unique()) != len(anarci_results):
            logger.warning(f"multiple results for {anarci_results[anarci_results['Id'].duplicated()]} is duplicated")
            anarci_results = anarci_results.sort_values("score").groupby("Id").head(1)

        # segment the region
        # anarci_results = anarci_results.add_segment_regions()
        return anarci_results

    def run_single(self, seq_id: str, seq: str, scfv=False) -> AnarciResults:
        """Run a single string sequence on an amino acid

        Parameters
        ----------
        seq_id : str
           the sequence_id of the string object, ex. "my_sequence"
        seq : str
            The string nucletodide sequence, ex. "EVQLQQSGAEVVRSG ..."

        Returns
        -------
            AnarchiResults Object
        """

        if not scfv:
            sequences = [(seq_id, seq)]

        return self._run(sequences)

    def run_multiple(self, seqrecords: List[SeqRecord], scfv=False) -> AnarciResults:
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
        if not isinstance(seqrecords, (list, type(SeqIO.FastaIO.FastaIterator), Generator)):
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

        if self.run_multiproc:
            # split a list into evenly sized chunks
            def chunks(list_to_split, n):
                return [list_to_split[i : i + n] for i in range(0, len(list_to_split), n)]

            # split sequences into chunks
            _sequences = chunks(_sequences, min(multiprocessing.cpu_count(), len(_sequences)))
            multiproc = multiprocessing.Pool()
            _results = pd.concat(multiproc.map(self._run, _sequences))
        else:
            _results = self._run(_sequences)
        return _results

    def run_dataframe(
        self,
        dataframe: pd.DataFrame,
        seq_id_field: Union[str, int],
        seq_field: Union[str, int],
        return_join=False,
    ) -> AnarciResults:
        """Pass dataframe and field and run airr.

        Parameters
        ----------
        dataframe : pd.DataFrame
            The input dataframe to run airr on

        seq_field: Union[str,int]
           The field in the dataframe to run airr on

        seq_id_field: Union[str,int]:
            The field that you want the "Sequence ID" in the airr table to correspond to.

        Returns
        -------
        AnarciResults
            AnarciResults object

        ToDo
        -------
        Default seq_id to be index. But have to account for it being a multi index
        """

        def _get_seq_generator():
            for seq_id, seq in zip(
                dataframe.reset_index()[seq_id_field],
                dataframe.reset_index()[seq_field],
            ):
                yield SeqRecord(id=str(seq_id), name=str(seq_id), description="", seq=Seq(str(seq)))

        if return_join:
            dataframe[seq_id_field] = dataframe[seq_id_field].astype(str)
            _df = self.run_multiple(_get_seq_generator())
            # convert seq id field to stry stince sequence_id is cast to string
            return dataframe.merge(
                _df,
                left_on=seq_id_field,
                right_on="sequence_id",
            )
        else:
            return self.run_multiple(_get_seq_generator())

    def run_file(self, file: Path) -> "AnarciResults":
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

            # run on fasta
        _results = self.run_multiple(list(SeqIO.parse(file, "fasta")))
        # if we had a file delete it
        if _filetype:
            os.unlink(tmpfile.name)
        return _results


if __name__ == "__main__":
    anarci_api = Anarci(scheme="chothia", region_assign="scdr")
    anarci_api.run_file("tests/integration/airr/fixtures/catnap_aa_heavy.fasta.gz")
    #     "MySweetAntibody",
    #     "EVQLLESGGGLVQPGGSLRLSCAASGFTFPVYNMAWVRQAPGKGLEWVSGIAHNGRNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKGHEISRFSRWSSFDYWGQGTLVTVSS",
    # )

    # pprint(result.segment_table_no_gaps)
    # anarci_api = Anarci(scheme="kabat", region_assign="scdr")
    # result = anarci_api.run_single(
    #     "MySweetAntibody",
    #     "DIQMTQSPSSLSASVGDRVTITCRPNQNIATYINWYQQKPGKAPKLLIYAASGLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQHSWEIPYTFGQGTKVEIK",
    # )
    # pprint(result.segment_table_no_gaps)
