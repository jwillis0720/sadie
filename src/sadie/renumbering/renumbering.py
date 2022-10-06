"""The NUMBERING Abstraction to make it usuable"""

# Std library
import gzip
import logging
import multiprocessing

import warnings
from multiprocessing import cpu_count
from pathlib import Path
from typing import Any, List, Tuple, Union, Generator

# third party
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

from sadie.renumbering.aligners import HMMER
from sadie.numbering import Numbering

from .exception import (
    NumberingDuplicateIdError,
    BadNumberingArgument,
)
from .result import NumberingResults
from .constants import NUMBERING_RESULTS
from sadie.numbering.scheme_numbering import scheme_numbering

logger = logging.getLogger("RENUMBERING")

# Get out of here with your partial codon warnigns
warnings.filterwarnings("ignore", "Partial codon")


class Error(Exception):
    """Base class for exceptions in this module."""


class Renumbering:

    hmmer = HMMER()
    numbering = Numbering()

    def __init__(
        self,
        aligner: str = "hmmer",
        scheme: str = "imgt",
        region_assign: str = "imgt",
        allowed_chain: List[str] = ["H", "K", "L"],
        assign_germline: bool = True,
        allowed_species: List[str] = ["human"],
        threshold: int = 80,
        run_multiproc: bool = True,
        num_cpus: int = cpu_count(),
        use_numbering_hmms: bool = False,
        *args,
        **kwargs,
    ):
        """HMMER v3 wrapper that runs hmmersearch using G3 built HMMs and numbering schema from Numbering.

        Parameters
        ----------
        aligner: str
            The aligner to use. Currently only hmmer is supported
        scheme : str, optional
            scheme of alignment, by default imgt,
            options: Chothia, Kabat, Martin (Extended Chothia), Aho
        region_assign : str, optional
            assigning frw1-cdr1-fwr2-cdr2-fwr3-cdr3-fwr4, by default "imgt"
            options: imgt, kabat, chothia
        allowed_chain : list, optional
            antibody chains, by default ["H", "K", "L"]
            all options: ["L", "H", "K", "A", "B", "G", "D"]
            options with numbering scheme + region: ["H", "K", "L"]
        assign_germline : bool, optional
            assign germline; falls back on hardcoded dict in germlines.py, by default True
        allowed_species : list, optional
            query sequnce only hits HMMs build with the species requested, by default ["human"]
            options: human, mouse, rat, rabbit, rhesus, pig, alpaca, dog, cat
        threshold : int, optional
            HMMER specific bitscore threshold determined by the best domain hit found, by default 80
            notes: anything over 160 is a likely position and anything below 80 is a likely false positive hit. Anything inbetween is open to interpretation.
        run_multiproc : bool, optional
            Runs each sequence as a input concurrenty, by default True
            notes: Each sequence in a fasta file give is run in parallel in the order given.
        num_cpus : int, optional
            number of cpus to use if multiple inputs are given, by default all cpus.
        use_numbering_hmms : bool, optional
            if True, will use only backup hmms, by default False
            note: these backup hmms are legacy from the ANARCI team and are not updated.
        *args, **kwargs  # for backwards compatibility options

        Raises
        ------
        NotImplementedError
            If the scheme + region assign combo is not implemented
        """
        self.scheme = scheme
        self.region_definition = region_assign
        self.allowed_chains = allowed_chain
        self.assign_germline = assign_germline
        self.allowed_species = allowed_species
        self.num_cpus = num_cpus
        self.run_multiproc = run_multiproc
        self.threshold_bit = threshold
        self.hmmer.use_numbering_hmms = use_numbering_hmms
        if not self.check_combination(self.scheme, self.region_definition):
            raise NotImplementedError(f"{self.scheme} with {self.region_definition} has not been implemented yet")

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
            raise BadNumberingArgument(definition, self.get_available_region_definitions())
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
        __future_schemes = ["martin", "aho"]
        if scheme.lower() not in ["imgt", "chothia", "kabat"]:
            logger.warning(f"need support for {__future_schemes} numbering schemes. See abysis")
            raise BadNumberingArgument(scheme, self.get_available_region_definitions())
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
            raise BadNumberingArgument(_diff, _allowed_chain)
        self._allowed_chains = allowed_chains

    @staticmethod
    def get_allowed_chains() -> List[str]:
        """Get the allowed chains options. Which chains can you align against

        Returns
        -------
        List
            A list of one letter codes that correspond to chains
        """
        _allowed_chain = ["H", "K", "L", "A", "B", "G", "D"]
        return _allowed_chain

    @property
    def allowed_species(self) -> List[str]:
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
            raise BadNumberingArgument(_diff, _allowed_species)
        self._allowed_species = allowed_species

    @staticmethod
    def get_allowed_species() -> List[str]:
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

    def _run(self, sequences: List[Tuple[str, str]]):
        """
        private method to run Numbering

        Parameters
        ----------
        sequences : List[Tuple]
            list or tuple of (Id, Sequence) pairs
                              e.g. [ ("seq1","EVQLQQSGAEVVRSG ..."),
                                     ("seq2","DIVMTQSQKFMSTSV ...")
        """
        # Perform the alignments of the sequences to the hmm database
        _alignments = self.hmmer.hmmsearch(
            sequences=sequences,
            species=self.allowed_species,
            chains=self.allowed_chains,
            bit_score_threshold=self.threshold_bit,
            limit=1,
            for_numbering=True,
        )

        # Check the numbering for likely very long CDR3s that will have been missed by the first pass.
        # Modify alignments in-place
        self.hmmer.check_for_j(
            sequences=sequences,
            alignments=_alignments,
            species=self.allowed_species,
            chains=self.allowed_chains,
        )

        # Apply the desired numbering scheme to all sequences
        _numbered, _alignment_details, _hit_tables = self.numbering.number_sequences_from_alignment(
            sequences,
            _alignments,
            scheme=self.scheme,
            allow=self.allowed_chains,
            assign_germline=self.assign_germline,
            allowed_species=self.allowed_species,
        )

        _summary = self.numbering.parsed_output(sequences, _numbered, _alignment_details)
        numbering_results = pd.DataFrame(_summary)

        if numbering_results.empty:
            return NumberingResults()

        # I really want to set the scheme and region in the constructor
        # https://stackoverflow.com/questions/66647680/subclassing-pandas-dataframe-and-setting-field-in-constuctor
        numbering_results = NumberingResults(
            numbering_results.astype(NUMBERING_RESULTS),
        )

        # Must set these schemes before we set the segments
        numbering_results["scheme"] = self.scheme
        numbering_results["region_definition"] = self.region_definition
        numbering_results["allowed_species"] = ",".join(self.allowed_species)
        numbering_results["allowed_chains"] = ",".join(self.allowed_chains)
        numbering_results = numbering_results._add_segment_regions()

        if len(numbering_results["Id"].unique()) != len(numbering_results):
            logger.warning(
                f"multiple results for {numbering_results[numbering_results['Id'].duplicated()]} is duplicated"
            )
            numbering_results = numbering_results.sort_values("score", ascending=False).groupby("Id").head(1)

        # segment the region
        # numbering_results = numbering_results.add_segment_regions()
        return numbering_results

    def run_single(self, seq_id: str, seq: str) -> NumberingResults:
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
        sequences = [(seq_id, seq)]

        return self._run(sequences)

    def run_multiple(self, seqrecords: List[SeqRecord], scfv=False) -> NumberingResults:
        """Run multiple seq records

        Parameters
        ----------
        seqrecords : List[SeqRecord]
            A list of sequence records of amino acids.

        Returns
        -------
            NUMBERINGResults - Holds many results

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
                raise NumberingDuplicateIdError(seq.id, 1)
            _sequences.append((seq.id, str(seq.seq)))
            _seen.add(seq.id)

        if self.run_multiproc:
            # split a list into evenly sized chunks
            def chunks(list_to_split: List[Any], n: int):
                return [list_to_split[i : i + n] for i in range(0, len(list_to_split), n)]

            try:
                # split sequences into chunks
                _sequences = chunks(_sequences, min(self.num_cpus, len(_sequences)))
                multiproc = multiprocessing.Pool()
                _results = pd.concat(multiproc.map(self._run, _sequences))

            # cleans up the pool
            finally:
                multiproc.close()
                multiproc.join()
        else:
            _results = self._run(_sequences)
        return _results

    def run_dataframe(
        self,
        dataframe: pd.DataFrame,
        seq_id_field: Union[str, int],
        seq_field: Union[str, int],
        return_join=False,
    ) -> NumberingResults:
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
        NumberingResults
            NumberingResults object

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

    def run_file(self, file: Path) -> "NumberingResults":
        """Run numbering annotator on a fasta file

        Parameters
        ----------
        file: Path
            The fasta file to run
        multi: Bool, defaults=True
            split and run file as multiprocess

        Returns
        -------
        NumberingResults
            Returns NumberingResults object

        Raises
        ------
        FileExistsError
            if file does not exist
        BadRequstedFileType
            if file is not fasta

        """
        file = Path(file)
        if file.is_file() is False:
            raise FileNotFoundError(f"{file} not found")

        if file.suffix == ".gz":
            with gzip.open(file, "rt") as handle:
                seqs = list(SeqIO.parse(handle, "fasta"))
        # Biopython natively handles bz2
        else:
            seqs = list(SeqIO.parse(file, "fasta"))

        return self.run_multiple(seqs)
