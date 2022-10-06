from functools import lru_cache
import gzip
import logging
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any, List, Tuple, Union
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy import nan
import pandas as pd

from sadie.antibody.exception import NumberingDecreasing
from sadie.renumbering.aligners import HMMER
from .exception import (
    NumberingDuplicateIdError,
    BadNumberingArgument,
)
from .constants import NUMBERING_RESULTS
from .germlines import all_germlines
from .result import NumberingResults
from .schemes import (
    number_kabat_heavy,
    number_kabat_light,
    number_chothia_heavy,
    number_chothia_light,
    # number_martin_heavy,
    # number_martin_light,
    number_imgt,
    # number_aho,
)

logger = logging.getLogger("RENUMBERING")

# Get out of here with your partial codon warnigns
warnings.filterwarnings("ignore", "Partial codon")


class Error(Exception):
    """Base class for exceptions in this module."""


class Renumbering:

    hmmer = HMMER()

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
        prioritize_cached_hmm: bool = True,
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
        prioritize_cached_hmm : bool, optional
            if True, will use cached hmms if they exist, by default True
            notes: A queried species and chain will be pulled from G3 if it exists, otherwise it will check a local backup before failing.
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
        self.allowed_species = allowed_species or ["human"]
        self.num_cpus = num_cpus
        self.run_multiproc = run_multiproc
        self.threshold_bit = threshold
        self.prioritize_cached_hmm = prioritize_cached_hmm
        self.hmmer.use_numbering_hmms = use_numbering_hmms

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
        # TODO: region numbering only supports H, K, L -- germline and aligner allow all chains.
        # _allowed_chain = ["H", "K", "L", "A", "B", "G", "D"]
        _allowed_chain = ["H", "K", "L"]
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

    @staticmethod
    def validate_numbering(xxx_todo_changeme, name_seq=[]):
        """
        Wrapper to do some basic validation of the numbering.

        Further validation could be done but at the moment we just check that the numbering indices are incremental (they should be)
        """
        (numbering, start, end) = xxx_todo_changeme
        name, seq = name_seq
        last = -1
        nseq = ""

        for (index, _), a in numbering:
            if index < last:
                raise NumberingDecreasing(name, "decreasing sequence count in the numbering")

                # , "Numbering was found to decrease along the sequence %s. Please report." % name
            last = index
            nseq += a.replace("-", "")

        assert nseq in seq.replace("-", ""), (
            "The algorithm did not number a contiguous segment for sequence %s. Please report" % name
        )

        return numbering, start, end

    @staticmethod
    def number_sequence_from_alignment(state_vector, sequence, scheme="imgt", chain_type=None):
        """
        Given you have an alignment. Give back the numbering

        @param state_vector: List of states from the hmm. Effectively these are imgt columns but CDR3 has not been redone.
        @param sequence: The original sequence string or list.
        @param scheme: The numbering scheme to apply
        @param chain_type: The type of chain to apply numbering for. Some schemes do not require this (IMGT). Others (e.g. Chothia/Wolfguy) do.

        @return: A list of numbering identifier / amino acids tuples over the domain that has been numbered. The indices of the start (inclusive) and end point (exclusive) in the sequence for the numbering
        """
        scheme = scheme.lower()
        if scheme == "imgt":
            return number_imgt(state_vector, sequence)
        elif scheme == "chothia":
            if chain_type == "H":
                return number_chothia_heavy(state_vector, sequence)
            elif chain_type in "KL":
                return number_chothia_light(state_vector, sequence)
            # else:
            #     raise AssertionError("Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
        elif scheme == "kabat":
            if chain_type == "H":
                return number_kabat_heavy(state_vector, sequence)
            elif chain_type in "KL":
                return number_kabat_light(state_vector, sequence)
            # else:
            #     raise AssertionError("Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
        # TODO: need scheme.py update for these formats
        # elif scheme == "martin":
        #     if chain_type == "H":
        #         return number_martin_heavy(state_vector, sequence)
        #     elif chain_type in "KL":
        #         return number_martin_light(state_vector, sequence)
        #     else:
        #         raise AssertionError("Unimplemented numbering scheme %s for chain %s" % (scheme, chain_type))
        # elif scheme == "aho":
        #     return number_aho(
        #         state_vector, sequence, chain_type
        #     )  # requires the chain type to heuristically put the CDR1 gap in position.

    @lru_cache(maxsize=None)
    def get_identity(self, state_sequence, germline_sequence):
        """
        Get the partially matched sequence identity between two aligned sequences.
        Partial in the sense that gaps can be in the state_sequence.
        """
        # Ensure that the sequences are the expected length
        assert len(state_sequence) == len(germline_sequence) == 128
        n, m = 0, 0
        for i in range(128):
            if germline_sequence[i] == "-":
                continue
            if state_sequence[i].upper() == germline_sequence[i]:
                m += 1
            n += 1

        return float(m) / n

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
            prioritize_cached_hmm=self.prioritize_cached_hmm,
            for_numbering=True,
        )

        # Check the numbering for likely very long CDR3s that will have been missed by the first pass.
        # Modify alignments in-place
        self.hmmer.check_for_j(
            sequences=sequences,
            alignments=_alignments,
            species=self.allowed_species,
            chains=self.allowed_chains,
            prioritize_cached_hmm=self.prioritize_cached_hmm,
        )

        # Apply the desired numbering scheme to all sequences
        _numbered, _alignment_details, _hit_tables = self.number_sequences_from_alignment(
            sequences,
            _alignments,
            scheme=self.scheme,
            allow=self.allowed_chains,
            assign_germline=self.assign_germline,
            allowed_species=self.allowed_species,
        )

        _summary = self.parsed_output(sequences, _numbered, _alignment_details)
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

        # TODO: Give multiple domain return options?
        # if len(numbering_results["Id"].unique()) != len(numbering_results):
        #     logger.warning(
        #         f"multiple results for {numbering_results[numbering_results['Id'].duplicated()]} is duplicated"
        #     )
        #     numbering_results = numbering_results.sort_values("score", ascending=False).groupby("Id").head(1)

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
                multiproc = Pool()
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

        # TODO: not used; seems to be more of a user end functionality that may not be common enough
        # if return_join:
        #     dataframe[seq_id_field] = dataframe[seq_id_field].astype(str)
        #     _df = self.run_multiple(_get_seq_generator())
        #     # convert seq id field to stry stince sequence_id is cast to string
        #     return dataframe.merge(
        #         _df,
        #         left_on=seq_id_field,
        #         right_on="sequence_id",
        #     )
        # else:
        #     return self.run_multiple(_get_seq_generator())

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

    def number_sequences_from_alignment(
        self,
        sequences,
        alignments,
        scheme="imgt",
        allow=set(["H", "K", "L", "A", "B", "G", "D"]),
        assign_germline=False,
        allowed_species=None,
    ):
        """
        Given a list of sequences and a corresponding list of alignments from run_hmmer apply a numbering scheme.
        """

        # Iteration over the sequence alignments performing the desired numbering
        numbered = []
        alignment_details = []
        hit_tables = []

        for i in range(len(sequences)):

            # Unpack
            hit_table, state_vectors, detailss = alignments[
                i
            ]  # We may have multiple domains per sequence (e.g. single chain fvs).

            # Iterate over all the domains in the sequence that have been recognised (typcially only 1 with the current hmms available)
            hit_numbered, hit_details = [], []
            for di in range(len(state_vectors)):
                state_vector = state_vectors[di]
                details = detailss[di]
                details["scheme"] = scheme
                details["query_name"] = sequences[i][0]

                # Only number things that are allowed. We still keep the alignment details and hit_table
                if state_vector and details["chain_type"] in allow:
                    # Do the numbering and validate (for development purposes)
                    try:
                        hit_numbered.append(
                            self.validate_numbering(
                                self.number_sequence_from_alignment(
                                    state_vector,
                                    sequences[i][1],
                                    scheme=scheme,
                                    chain_type=details["chain_type"],
                                ),
                                sequences[i],
                            )
                        )
                    except NumberingDecreasing:
                        warnings.warn(f"Skipping {details['query_name']}", UserWarning)
                        continue

                    if assign_germline:
                        details["germlines"] = self.run_germline_assignment(
                            state_vector,
                            sequences[i][1],
                            details["chain_type"],
                            allowed_species=allowed_species,
                        )
                    hit_details.append(details)

            if hit_numbered:
                numbered.append(hit_numbered)
                alignment_details.append(hit_details)
            else:
                numbered.append(None)
                alignment_details.append(None)
            hit_tables.append(hit_table)

        return numbered, alignment_details, hit_tables

    def parsed_output(self, sequences, numbered, details):
        """
        Write numbered sequences to dataframe
        Kappa and Lambda chains are written to the same file
        The sequences will written aligned to the numbering scheme. Gaps in the sequences with respect to the alignment are written
        as a '-'
        @param sequences: List of name, sequence tuples
        @param numbered: Numbered sequences in the same order as the sequences list.
        @param details: List of alignment details in the same order as the sequences list.
        """

        chain_types = {}
        pos_ranks = {}
        all_pos = {}
        _lc = {"K": "KL", "L": "KL"}

        # Divide the set into chain types and find how to order the numbering for each type.
        for i in range(len(sequences)):  # Iterate over entries
            if numbered[i] is None:
                continue

            for j in range(len(numbered[i])):  # Iterate over domains.
                # Record the chain type index
                c = details[i][j]["chain_type"]
                c = _lc.get(c, c)  # Consider lambda and kappa together.
                chain_types.setdefault(c, []).append((i, j))
                if c not in pos_ranks:
                    pos_ranks[c] = {}
                    all_pos[c] = set()

                # Update the insertion order for the scheme. i.e. is it A B C or C B A (e.g. imgt 111 and 112 repectively)
                tmp_p = -1
                r = 0
                for p, _ in numbered[i][j][0]:
                    if p[0] != tmp_p:
                        tmp_p = p[0]
                        r = 0
                    else:
                        r += 1
                    pos_ranks[c][p] = max(r, pos_ranks[c].get(p, r))
                    all_pos[c].add(p)

        summary_dataframes = []
        # Write a new file for each chain type. Kappa and lambda are written together as light chains.
        for cts in ["H", "KL", "A", "B", "G", "D"]:
            if cts in chain_types:

                # Sort the positions by index and insertion order
                positions = sorted(all_pos[cts], key=lambda p: (p[0], pos_ranks[cts][p]))

                # Header line
                fields = [
                    "Id",
                    "sequence",
                    "domain_no",
                    "hmm_species",
                    "chain_type",
                    "e-value",
                    "score",
                    "seqstart_index",
                    "seqend_index",
                    "identity_species",
                    "v_gene",
                    "v_identity",
                    "j_gene",
                    "j_identity",
                ]
                numbering_ = [("%d%s" % (p)).strip() for p in positions]
                # print(",".join(fields), file=out)

                # Iterate over the domains identified
                for i, j in chain_types[cts]:
                    j_gene = nan
                    j_gene_score = nan
                    if details[i][j].get("germlines", {}).get("j_gene", [0])[0]:
                        j_gene = details[i][j].get("germlines", {}).get("j_gene", [["", ""], 0])[0][1]
                        j_gene_score = "%.2f" % details[i][j].get("germlines", {}).get("j_gene", [["", ""], 0])[1]
                    line = [
                        sequences[i][0].replace(",", " "),
                        sequences[i][1],
                        str(j),
                        details[i][j].get("species", ""),
                        details[i][j].get("chain_type", ""),
                        str(details[i][j].get("evalue", "")),
                        str(details[i][j].get("bitscore", "")),
                        str(numbered[i][j][1]),
                        str(numbered[i][j][2]),
                        details[i][j].get("germlines", {}).get("v_gene", [["", ""], 0])[0][0],
                        details[i][j].get("germlines", {}).get("v_gene", [["", ""], 0])[0][1],
                        "%.2f" % details[i][j].get("germlines", {}).get("v_gene", [["", ""], 0])[1],
                        j_gene,
                        j_gene_score,
                    ]

                    # Hash the numbering. Insertion order has been preserved in the positions sort.
                    d = dict(numbered[i][j][0])
                    seq_ = [d.get(p, "-") for p in positions]
                    numbering_ = [[int(p[0]), p[1]] for p in positions]
                    assert len(line) == len(fields)
                    _df = pd.Series(line, index=fields)
                    _df["Chain"] = cts
                    _df["Numbering"] = list(map(lambda x: x[0], numbering_))
                    _df["Insertion"] = list(map(lambda x: x[1].strip(), numbering_))
                    _df["Numbered_Sequence"] = seq_
                    summary_dataframes.append(_df)

        return summary_dataframes

    def run_germline_assignment(self, state_vector, sequence, chain_type, allowed_species=None):
        """
        Find the closest sequence identity match.
        """
        genes = {"v_gene": [None, None], "j_gene": [None, None]}

        # Extract the positions that correspond to match (germline) states.
        state_dict = dict(((i, "m"), None) for i in range(1, 129))
        state_dict.update(dict(state_vector))
        state_sequence = "".join(
            [sequence[state_dict[(i, "m")]] if state_dict[(i, "m")] is not None else "-" for i in range(1, 129)]
        )

        # Iterate over the v-germline sequences of the chain type of interest.
        # The maximum sequence identity is used to assign the germline
        if chain_type in all_germlines["V"]:
            _allowed = []
            for sp in allowed_species:
                _allowed.append(sp)
            seq_ids = {}
            for species in allowed_species:
                for gene, germline_sequence in all_germlines["V"][chain_type][species].items():
                    seq_ids[(species, gene)] = self.get_identity(state_sequence, germline_sequence)
            genes["v_gene"][0] = max(seq_ids, key=lambda x: seq_ids[x])
            genes["v_gene"][1] = seq_ids[genes["v_gene"][0]]

            # Use the assigned species for the v-gene for the j-gene.
            # This assumption may affect exotically engineered abs but in general is fair.
            species = genes["v_gene"][0][0]
            if chain_type in all_germlines["J"]:
                if species in all_germlines["J"][chain_type]:
                    seq_ids = {}
                    for gene, germline_sequence in all_germlines["J"][chain_type][species].items():
                        seq_ids[(species, gene)] = self.get_identity(state_sequence, germline_sequence)
                    genes["j_gene"][0] = max(seq_ids, key=lambda x: seq_ids[x])
                    genes["j_gene"][1] = seq_ids[genes["j_gene"][0]]

        return genes
