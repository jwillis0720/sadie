import json
import logging
import warnings
from pathlib import Path
from typing import Tuple, Union

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

# from Bio.Align import substitution_matrices
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio.SubsMat import MatrixInfo as matlist
from numpy import nan
from Levenshtein._levenshtein import distance

from .constants import IGBLAST_AIRR
from .genbank import GenBank, GenBankFeature

logger = logging.getLogger("AIRRTable")
blosum_matrix = matlist.blosum62


class Error(Exception):
    """Base class for exceptions in this module."""


class MissingAirrColumns(Error):
    """Exception raised for not finiding the igblast module

    Attributes:
    """

    def __init__(self, missing_columns):
        super().__init__()
        self.missing_columns = missing_columns

    def __str__(self):
        return "Must have all AIRR columns defined, missing {}".format(self.missing_columns)


def _get_v_aa_distance(X) -> int:
    """Helper method for getting the amino acid distance

    Parameters
    ----------
    X : tuple
        (v_sequence_aligbment_aa,v_sequence_germline_alignment_aa)

    Returns
    -------
    float
        percentabe between germ and mature. i.e #(mutations+indels)/len(max(germ,mat))
    """
    # X = X.fillna("")
    if X.isna().any():
        return
    _mature = X["v_sequence_alignment_aa"]
    _germ = X["v_germline_alignment_aa"]
    return (distance(_mature, _germ) / max(len(_mature), len(_germ))) * 100


class AirrTable:
    """
        Object oriented Airr table
        compliance with - https://docs.airr-community.org/_/downloads/en/v1.2.1/pdf/

        Raises
        ------
    MissingAirrColumns
            If there is any missing columns from the airr standards


        Examples
        --------

        Pass in pandas dataframe
        >>> df = pd.read_csv("tests/fixtures/heavy_lev_sample.csv.gz")
        >>> at = AirrTable(df)

        The pandas table is held in the the table attribute
        >>> at.table.columns
        >>> ['sequence_id', 'cdr1', 'cdr1_aa', 'cdr1_end', 'cdr1_start', 'cdr2',
            'cdr2_aa', 'cdr2_end', 'cdr2_start', 'cdr3', 'cdr3_aa',
            'cdr3_aa_length', 'cdr3_end', 'cdr3_start', 'complete_vdj',
            'd_alignment_end', 'd_alignment_start', 'd_call', 'd_cigar', 'd_family',
            'd_germline_alignment', 'd_germline_alignment_aa', 'd_germline_end',
            'd_germline_start', 'd_identity', 'd_score', 'd_sequence_alignment',
            'd_sequence_alignment_aa', 'd_sequence_end', 'd_sequence_start',
            'd_support', 'fwr1', 'fwr1_aa', 'fwr1_end', 'fwr1_start', 'fwr2',
            'fwr2_aa', 'fwr2_end', 'fwr2_start', 'fwr3', 'fwr3_aa', 'fwr3_end',
            'fwr3_start', 'fwr4', 'fwr4_aa', 'fwr4_end', 'fwr4_start',
            'germline_alignment', 'germline_alignment_aa', 'j_alignment_end',
            'j_alignment_start', 'j_call', 'j_cigar', 'j_family',
            'j_germline_alignment', 'j_germline_alignment_aa', 'j_germline_end',
            'j_germline_start', 'j_identity', 'j_score', 'j_sequence_alignment',
            'j_sequence_alignment_aa', 'j_sequence_end', 'j_sequence_start',
            'j_support', 'junction', 'junction_aa', 'junction_aa_length',
            'junction_length', 'locus', 'np1', 'np1_length', 'np2', 'np2_length',
            'productive', 'rev_comp', 'sequence', 'sequence_alignment',
            'sequence_alignment_aa', 'stop_codon', 'v_alignment_end',
            'v_alignment_start', 'v_call', 'v_cigar', 'v_family',
            'v_germline_alignment', 'v_germline_alignment_aa', 'v_germline_end',
            'v_germline_start', 'v_identity', 'v_score', 'v_sequence_alignment',
            'v_sequence_alignment_aa', 'v_sequence_end', 'v_sequence_start',
            'v_support', 'vj_in_frame', 'chain', 'species'],

        Can also use static methods
        >>> at = AirrTable.read_csv("tests/fixtures/kappa_lev_sample.csv.gz")
    """

    def __init__(self, dataframe: pd.DataFrame, cast=True, infer_germline=True):
        """AirrTable consturctor

        Parameters
        ----------
        dataframe : pd.DataFrame
           pandas dataframe. Must have all Airr complient columns

        cast: bool
            cast dtypes to save memory
        infer_germline: bool
            infer germline sequences
        """
        self.cast = cast
        self.infer = infer_germline
        self.table = dataframe

    @property
    def table(self) -> pd.DataFrame:
        """Getter method for table

        Returns
        -------
        pd.DataFrame
           AirrTable
        """
        return self._table

    @table.setter
    def table(self, dataframe: pd.DataFrame):
        """Set Table

        Parameters
        ----------
        dataframe : pd.DataFrame
            set pandas table

        Raises
        ------
        MissingAirrColumns
           check complient columns
        """
        self._table = dataframe
        self._table.drop([i for i in self._table.columns if "Unnamed" in i], axis=1, inplace=True)
        self._table.loc[:, "note"] = ""

        # a check that allows us to see if the junction presentation is wrong
        liable_sequences = self._table[
            self._table[
                [
                    "fwr1_aa",
                    "cdr1_aa",
                    "fwr2_aa",
                    "cdr2_aa",
                    "fwr3_aa",
                    "cdr3_aa",
                    "fwr4_aa",
                ]
            ].apply(self._check_j_gene_liability, axis=1)
        ]

        # If these aren't productive, who cares
        if not liable_sequences.empty:
            logger.debug(f"Caution - sequences {list(liable_sequences['sequence_id'])} may need manual inspections")
            self._table.loc[liable_sequences.index, "note"] = "liable"

        missing_columns = set(IGBLAST_AIRR.keys()).difference(self._table.columns)
        if missing_columns:
            raise MissingAirrColumns(missing_columns)

        if self.cast:
            self._table = self._table.astype(IGBLAST_AIRR)

        # set those bool values
        self._table["productive"] = (
            self._table["productive"]
            .map(
                {"T": True, "F": False, True: True, False: False, nan: False},
                na_action="ignore",
            )
            .astype(bool)
        )
        self._table["stop_codon"] = (
            self._table["stop_codon"]
            .map(
                {"T": True, "F": False, True: True, False: False, nan: False},
                na_action="ignore",
            )
            .astype(bool)
        )
        self._table["vj_in_frame"] = (
            self._table["vj_in_frame"]
            .map(
                {"T": True, "F": False, True: True, False: False, nan: False},
                na_action="ignore",
            )
            .astype(bool)
        )
        self._table["rev_comp"] = (
            self._table["rev_comp"]
            .map(
                {"T": True, "F": False, True: True, False: False, nan: False},
                na_action="ignore",
            )
            .astype(bool)
        )
        self._table["v_frameshift"] = (
            self._table["v_frameshift"]
            .map(
                {"T": True, "F": False, True: True, False: False, nan: False},
                na_action="ignore",
            )
            .astype(bool)
        )
        self._table["complete_vdj"] = (
            self._table["complete_vdj"]
            .map(
                {"T": True, "F": False, True: True, False: False, nan: False},
                na_action="ignore",
            )
            .astype(bool)
        )
        self._table["locus"] = self._table["locus"].astype(
            pd.CategoricalDtype(categories=["IGH", "IGL", "IGK"], ordered=True)
        )

        self._table["vdj_nt"] = self._table[["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]].apply(
            lambda x: "".join([str(i) for i in x if not isinstance(i, float)]), axis=1
        )
        self._table["vdj_aa"] = self._table[
            [
                "fwr1_aa",
                "cdr1_aa",
                "fwr2_aa",
                "cdr2_aa",
                "fwr3_aa",
                "cdr3_aa",
                "fwr4_aa",
            ]
        ].apply(lambda x: "".join([str(i) for i in x if not isinstance(i, float)]), axis=1)

        # Insert the top call for each one
        for call in ["v_call", "d_call", "j_call"]:
            if call not in self._table.columns:
                continue
            if f"{call}_top" in self._table.columns:
                self._table.drop(f"{call}_top", inplace=True, axis=1)
            self._table.insert(
                self._table.columns.get_loc(call), f"{call}_top", self._table[call].str.split(",").str.get(0)
            )

        # indels are not handled in the amino acid sequence
        indel_indexes = self._table[
            (self._table["sequence_alignment_aa"].str.len() != self._table["germline_alignment_aa"].str.len())
            & (~self._table["sequence_alignment_aa"].isna())
            & (~self._table["germline_alignment_aa"].isna())
        ].index
        logger.info(f"Have {len(indel_indexes)} possible indels that are not in amino acid alignments")

        # indels are not handled in the amino acid sequence
        v_indel_indexes = self._table[
            (self._table["v_sequence_alignment_aa"].str.len() != self._table["v_germline_alignment_aa"].str.len())
            & (~self._table["v_sequence_alignment_aa"].isna())
            & (~self._table["v_germline_alignment_aa"].isna())
        ].index
        logger.info(f"Have {len(v_indel_indexes)} possible v-gene indels that are not in amino acid alignments")

        def _correct_alignment(X, field_1, field_2):
            alignment_aa_1 = X[field_1]
            alignment_aa_2 = X[field_2]
            if any([isinstance(alignment_aa_1, float), isinstance(alignment_aa_2, float)]):
                return pd.Series({field_1: alignment_aa_1, field_2: alignment_aa_2})

            # get alignments
            try:
                alignments = pairwise2.align.globalds(
                    alignment_aa_1,
                    alignment_aa_2,
                    blosum_matrix,
                    -12,
                    -4,
                    penalize_extend_when_opening=True,
                    penalize_end_gaps=False,
                )
            except SystemError:
                logger.debug(f"System error most likely due to * in germline alignment...falling back {alignment_aa_2}")
                alignments = pairwise2.align.globalms(
                    alignment_aa_1,
                    alignment_aa_2,
                    4,
                    -1,
                    -12,
                    -4,
                    penalize_extend_when_opening=True,
                    penalize_end_gaps=False,
                )
            alignment_1_aa_corrected = alignments[0][0]
            alignment_2_aa_corrected = alignments[0][1]
            return pd.Series(
                {
                    field_1: alignment_1_aa_corrected,
                    field_2: alignment_2_aa_corrected,
                }
            )

        if not indel_indexes.empty:
            correction_alignments = self._table.loc[indel_indexes, :].apply(
                lambda x: _correct_alignment(x, "sequence_alignment_aa", "germline_alignment_aa"), axis=1
            )
            # Correction in place
            self._table.update(correction_alignments)
            self._table.loc[indel_indexes, "alignment_correct"] = True

        if not v_indel_indexes.empty:
            correction_v_alignments = self._table.loc[v_indel_indexes, :].apply(
                lambda x: _correct_alignment(x, "v_sequence_alignment_aa", "v_germline_alignment_aa"), axis=1
            )
            # Correction in place
            self._table.update(correction_v_alignments)
            self._table.loc[v_indel_indexes, "alignment_correct"] = True

            logger.debug("Corrected gapped sequences")
        # get mutation frequency rather than identity
        # v identy are in percentage
        self._table.loc[:, "v_mutation"] = self._table["v_identity"].apply(lambda x: (100 - x))
        self._table.loc[:, "d_mutation"] = self._table["d_identity"].apply(lambda x: (100 - x))
        self._table.loc[:, "j_mutation"] = self._table["j_identity"].apply(lambda x: (100 - x))

        # mutation frequency in aa
        self._table.loc[:, "v_mutation_aa"] = self._table[["v_sequence_alignment_aa", "v_germline_alignment_aa"]].apply(
            lambda x: _get_v_aa_distance(x), axis=1
        )

        if self.infer:
            self._table["vdj_igl"] = self._table.apply(self._get_igl, axis=1)
            self._table.loc[self._table[self._table["vdj_igl"].isna()].index, "note"] = "liable"
            self._table["igl_mut_aa"] = self._table[["vdj_aa", "vdj_igl"]].apply(self._get_diff, axis=1)

        # finally add what is airr columns and what is not
        self._non_airr_columns = list(set(self._table.columns) - set(IGBLAST_AIRR.keys()))
        self._airr_columns = list(set(self._table.columns).intersection(IGBLAST_AIRR.keys()))

    def _get_diff(self, X):
        "get character levenshtrein distance"
        first = X[0]
        second = X[1]
        if not first or not second or isinstance(first, float) or isinstance(second, float):
            return
        return (distance(first, second) / max(len(first), len(first))) * 100

    def _get_igl(self, row: pd.Series) -> str:
        """Get infered germline sequenxe

        Parameters
        ----------
        row : pd.Series
            A row from the airr table
        Returns
        -------
        str
            the igl sequecne
        """

        # get germline components
        v_germline = row.v_germline_alignment_aa
        full_germline = row.germline_alignment_aa
        if isinstance(v_germline, float):
            return
        cdr3_j_germline = full_germline[len(v_germline) :]

        # get mature components
        v_mature = row.v_sequence_alignment_aa
        full_mature = row.sequence_alignment_aa
        cdr3_j_mature = full_mature[len(v_mature) :]

        # if the mature and cdr3 are not the same size
        # this will happen on non-productive
        if len(cdr3_j_mature) != len(cdr3_j_germline):
            logger.debug(f"{row.name} - strange iGL")
            return

            # # quick aligment
            # _alignments = align.globalxs(cdr3_j_mature, cdr3_j_germline, -10, -1)
            # cdr3_j_mature, cdr3_j_germline = _alignments[0][0], _alignments[0][1]

        iGL_cdr3 = ""
        for mature, germline in zip(cdr3_j_mature, cdr3_j_germline):
            if germline == "X" or germline == "-":
                iGL_cdr3 += mature
                continue
            iGL_cdr3 += germline

        full_igl = v_germline.replace("-", "") + iGL_cdr3.replace("-", "")
        return full_igl

    @property
    def non_airr_columns(self) -> list:
        """Getter for non airr columns that available in table

        Returns
        -------
        list
           columns that are not apart of the airr complicence
        """
        return sorted(self._non_airr_columns)

    @property
    def airr_columns(self) -> list:
        """Getter for airr complient columns

        Returns
        -------
        list
           columns that are apart of the airr complience
        """
        return sorted(self._airr_columns)

    @property
    def productive(self) -> pd.DataFrame:
        """Get productive airr dataframe

        Returns
        -------
        pd.DataFrame
            Dataframe of just productive sequences
        """
        return self._table[self._table["productive"]]

    @property
    def index(self) -> pd.Index:
        """

        Returns
        -------
        pd.Index
           returns pandas index
        """
        return self._table.index

    @property
    def iloc(self) -> pd.core.indexing._iLocIndexer:
        return self._table.iloc

    @property
    def loc(self) -> pd.core.indexing._LocIndexer:
        return self._table.loc

    @property
    def sanitized_antibodies(self) -> pd.DataFrame:
        """Getter method for airr entires that have  all antibody segments

        Returns
        -------
        pd.DataFrame
            dataframe of sanitized sequences
        """
        _sanitized = self[
            ~(
                (self["fwr1_aa"].isna())
                | (self["cdr1_aa"].isna())
                | (self["fwr2_aa"].isna())
                | (self["cdr2_aa"].isna())
                | (self["fwr3_aa"].isna())
                | (self["cdr3_aa"].isna())
                | (self["fwr4_aa"].isna())
            )
        ]
        # ensures that we start aligning at the very first codon
        return _sanitized[
            (_sanitized["v_germline_start"] == 1)
            & (_sanitized["productive"])
            & (_sanitized["vj_in_frame"])
            & (_sanitized["complete_vdj"])
            & (_sanitized["fwr4_aa"].str.len() >= 8)
        ]

    @staticmethod
    def parse_row_to_genbank(row: Tuple[int, pd.core.series.Series], suffix="") -> SeqRecord:
        single_seq = row[1]
        # these should be joined by sequence, so even if its heavy or light suffix, sequence should be in common
        sequence = single_seq["sequence"]

        # grab v,d and j calls
        v_call_string = "v_call" + suffix
        v_call = single_seq[v_call_string]
        d_call_string = "d_call" + suffix
        d_call = single_seq[d_call_string]
        j_call_string = "j_call" + suffix
        j_call = single_seq[j_call_string]

        # defaults
        top_v = ""
        top_d = ""
        top_j = ""
        other_v = ""
        other_d = ""
        other_j = ""
        if not isinstance(v_call, float):
            top_v = v_call.split(",")[0]
            other_v = v_call.split(",")[1:]

        if not isinstance(d_call, float):
            if d_call:
                top_d = d_call.split(",")[0]
                other_d = d_call.split(",")[1:]

        if not isinstance(j_call, float):
            top_j = j_call.split(",")[0]
            other_j = j_call.split(",")[1:]

        # Get start and end string
        _fw1_start_string = "fwr1_start" + suffix
        _fw2_start_string = "fwr2_start" + suffix
        _fw3_start_string = "fwr3_start" + suffix
        _fw4_start_string = "fwr4_start" + suffix
        _fw1_end_string = "fwr1_end" + suffix
        _fw2_end_string = "fwr2_end" + suffix
        _fw3_end_string = "fwr3_end" + suffix
        _fw4_end_string = "fwr4_end" + suffix

        # CDR start and end
        _cdr1_start_string = "cdr1_start" + suffix
        _cdr2_start_string = "cdr2_start" + suffix
        _cdr3_start_string = "cdr3_start" + suffix
        _cdr1_end_string = "cdr1_end" + suffix
        _cdr2_end_string = "cdr2_end" + suffix
        _cdr3_end_string = "cdr3_end" + suffix

        # Gene Segment
        _v_segment_start_string = "v_sequence_start" + suffix
        _v_segment_end_string = "v_sequence_end" + suffix
        _d_segment_start_string = "d_sequence_start" + suffix
        _d_segment_end_string = "d_sequence_end" + suffix
        _j_segment_start_string = "j_sequence_start" + suffix
        _j_segment_end_string = "j_sequence_end" + suffix

        # Get starts and ends _fw1_start = single_seq[_fw1_start_string]
        _fw1_start = single_seq[_fw1_start_string] - 1
        _fw2_start = single_seq[_fw2_start_string] - 1
        _fw3_start = single_seq[_fw3_start_string] - 1
        _fw4_start = single_seq[_fw4_start_string] - 1
        _fw1_end = single_seq[_fw1_end_string]
        _fw2_end = single_seq[_fw2_end_string]
        _fw3_end = single_seq[_fw3_end_string]
        _fw4_end = single_seq[_fw4_end_string] - 1

        # Get starts and ends of cdrs
        _cdr1_start = single_seq[_cdr1_start_string] - 1
        _cdr2_start = single_seq[_cdr2_start_string] - 1
        _cdr3_start = single_seq[_cdr3_start_string] - 1
        _cdr1_end = single_seq[_cdr1_end_string]
        _cdr2_end = single_seq[_cdr2_end_string]
        _cdr3_end = single_seq[_cdr3_end_string]

        # get v,d and j segments
        _v_segment_start = single_seq[_v_segment_start_string] - 1
        _v_segment_end = single_seq[_v_segment_end_string]
        _d_segment_start = single_seq[_d_segment_start_string] - 1
        _d_segment_end = single_seq[_d_segment_end_string]
        _j_segment_start = single_seq[_j_segment_start_string] - 1
        _j_segment_end = single_seq[_j_segment_end_string]

        # Locus identifier
        _locus_string = "locus" + suffix
        _productive_string = "productive" + suffix
        _complete_vdj_string = "complete_vdj" + suffix
        _locus = str(single_seq[_locus_string])
        _productive = single_seq[_productive_string]
        _complete_vdj = single_seq[_complete_vdj_string]

        # Collect species
        if "species" in single_seq.keys():
            species = single_seq["species"]
        else:
            warnings.warn(f"species is not in {single_seq.keys()}")
            species = "Unknown"

        # setup genbank object
        gb = GenBank(sequence, str(single_seq["sequence_id"]), description="AIRR annotation")

        # default qualifier dictionaries
        qd_v = {"gene": top_v, "species": species, "other_vgene": other_v}
        qd_d = {"gene": top_d, "species": species, "other_dgene": other_d}
        qd_j = {"gene": top_j, "species": species, "other_jgene": other_j}

        # FW1
        if all([isinstance(_fw1_start, int), isinstance(_fw1_end, int)]):
            feature = GenBankFeature(_fw1_start, _fw1_end, "FWR1", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # CDR1
        if all([isinstance(_cdr1_start, int), isinstance(_cdr1_end, int)]):
            feature = GenBankFeature(_cdr1_start, _cdr1_end, "CDR1", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # FW2
        if all([isinstance(_fw2_start, int), isinstance(_fw2_end, int)]):
            feature = GenBankFeature(_fw2_start, _fw2_end, "FWR2", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # CDR2
        if all([isinstance(_cdr2_start, int), isinstance(_cdr2_end, int)]):
            feature = GenBankFeature(_cdr2_start, _cdr2_end, "CDR2", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # FW3
        if all([isinstance(_fw3_start, int), isinstance(_fw3_end, int)]):
            feature = GenBankFeature(_fw3_start, _fw3_end, "FWR3", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # CDR3
        if all([isinstance(_cdr3_start, int), isinstance(_cdr3_end, int)]):
            feature = GenBankFeature(_cdr3_start, _cdr3_end, "CDR3")
            gb.add_feature(feature)

        # FW4
        if all([isinstance(_fw4_start, int), isinstance(_fw4_end, int)]):
            feature = GenBankFeature(_fw4_start, _fw4_end, "FWR4", qualifier_dict=qd_j)
            gb.add_feature(feature)

        # VGene
        if all([isinstance(_v_segment_start, int), isinstance(_v_segment_end, int)]):
            qd_v = {"gene": top_v, "species": species}
            feature = GenBankFeature(_v_segment_start, _v_segment_end, "VGene", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # DGene
        if all([isinstance(_d_segment_start, int), isinstance(_d_segment_end, int)]):
            qd_d = {"gene": top_d, "species": species}
            feature = GenBankFeature(_d_segment_start, _d_segment_end, "DGene", qualifier_dict=qd_d)
            gb.add_feature(feature)

        # JGene
        if all([isinstance(_j_segment_start, int), isinstance(_j_segment_end, int)]):
            qd_j = {"gene": top_j, "species": species}
            feature = GenBankFeature(_j_segment_start, _j_segment_end, "JGene", qualifier_dict=qd_j)
            gb.add_feature(feature)

        # Locus
        if all([isinstance(_v_segment_start, int), isinstance(_fw4_end, int)]):
            _p = "Productive"
            _vdj = "Complete"
            if not _productive:
                _p = "Non-Productive"
                # _locus = "_".join([_locus, _p])
            if not _complete_vdj:
                _vdj = "Non-Complete"
            qd_locus = {
                "locus": _locus,
                "species": species,
                "productive": _p,
                "VDJ": _vdj,
            }
            feature = GenBankFeature(_v_segment_start, _fw4_end, _locus, qualifier_dict=qd_locus)
            gb.add_feature(feature)
        return gb.record

    # IO
    def get_genbank(self):
        _genbanks = []
        # I'm no sure there is a better way to do this other than go through one by one
        for row in self._table.iterrows():
            _genbanks.append(AirrTable.parse_row_to_genbank(row))
        return _genbanks

    def to_genbank(self, filename: str, compression="gzip"):
        """write airrtable to genbank file

        Parameters
        ----------
        filename : str
            file name string path
        compression: ['gzip']
        """
        _records = self.get_genbank()
        SeqIO.write(_records, filename, "genbank")

    def get_json(self, indent=False) -> json:
        """Get json string of AirrTable

        Parameters
        ----------
        indent : bool, optional
            indent json string representation, by default False

        Returns
        -------
        json
            AirrTable as a json table
        """
        return self._table.to_json(indent=indent, orient="records")

    def to_json(self, filename: str, compression="gzip"):
        """write airrtable to json file

        Parameters
        ----------
        filename : str
            file name string path
        """
        self._table.to_json(filename, orient="records")

    def to_csv(self, path_or_buf=None, compression="infer", *args, **kwargs):
        """write airrtable to csv

        Overloaded pandas.to_csv(*args, **kwargs)
        """
        return self._table.to_csv(path_or_buf, **kwargs)

    def to_feather(self, path_or_buf=None, compression="infer", *args, **kwargs):
        """write airrtable to csv

        Overloaded pandas.to_csv(*args, **kwargs)
        """
        return self._table.to_feather(path_or_buf, **kwargs)

    @staticmethod
    def from_json(json_object: json) -> "AirrTable":
        """take in json object serialized and return AirrTable

        Returns
        -------
        AirrTable
            AirrTable object"""
        return AirrTable(pd.read_json(json_object))

    @staticmethod
    def read_csv(file: str) -> "AirrTable":
        """Read a csv file and return an airrtable object

        Returns
        -------
        AirrTable
            AirrTable object"""
        return AirrTable(pd.read_csv(file))

    @staticmethod
    def read_json(*args, **kwargs):
        """Read a json file and return an airrtable object

        Returns
        -------
        AirrTable
            AirrTable object"""
        return AirrTable(pd.read_json(*args, **kwargs))

    def _check_j_gene_liability(self, X):
        fw1 = X[0]
        cdr1 = X[1]
        fw2 = X[2]
        cdr2 = X[3]
        fw3 = X[4]
        cdr3 = X[5]
        fw4 = X[6]

        # If all fw and cdrs are reolsved
        if all([isinstance(i, str) for i in [fw1, cdr1, fw2, cdr2, fw3, cdr3, fw4]]):
            return False
        # if we have cdr2 and fwr3 and dont have cdr3, we probably have a problem
        if all([isinstance(i, str) for i in [cdr2, fw3]]) and (isinstance(cdr3, float) or cdr3 == ""):
            return True

        # fw1-cdr2 might be nan but we still get a fwr3-cdr3
        if any([isinstance(i, float) for i in [fw1, cdr1, fw2, cdr2]]) and (
            all([isinstance(cdr3, str), isinstance(fw3, str)])
        ):
            return False

        # we didn't get the fw4
        if all([isinstance(fw3, str), isinstance(cdr3, str)]) and (isinstance(fw4, float) or fw4 == ""):
            return True

        logger.debug(f" You shouldn't have gotten here {[fw1,cdr1, fw2, cdr2, fw3, cdr3, fw4]}")
        return True

    def get_sanitized_antibodies(self):
        """
        Get only entries with full length productive reads
        """
        return self.sanitized_antibodies

    def set_index(self, *args, **kwargs):
        return self._table.set_index(*args, **kwargs)

    def write_fasta(self, id_field: str, sequence_field: str, file_out: Path):
        """given an id field and sequence field, write out dataframe to fasata

        Parameters
        ----------
        id_field : str
            the id field from the dataframe to use in the fasta header
        sequence_field : str
            the seq field to use as the sqeuence
        file_out : Path
            the file output path to fasta
        """
        with open(file_out, "w") as f:
            for _, row in self._table.iterrows():
                f.write(f">{row[id_field]}\n{row[sequence_field]}\n")

    @property
    def empty(self) -> bool:
        return self._table.empty

    def __getitem__(self, col):
        return self._table[col]

    def __setitem__(self, key, val):
        self._table[key] = val

    def __len__(self):
        return len(self._table)

    def __repr__(self):
        return self._table.__repr__()

    def _repr_html_(self):
        """HTML reprsentation for jupyter notebooks"""
        return self._table._repr_html_()


class ScfvAirrTable:
    """Class for joined airr tables that correspond to heavy and light entries"""

    def __init__(
        self,
        heavy_airr_table: AirrTable,
        light_airr_table: AirrTable,
        join_on=("sequence_id", "sequence"),
    ):
        """constructor for heavy light chain airr table

        Args:
            heavy_airr_table (AirrTable): Heavy Chain Airr Table
            light_airr_table (AirrTable): Light Chain Airr Table

        Raises:
            TypeError: If heavy and light are not instances of airr table
        """
        if not all(
            [
                isinstance(heavy_airr_table, AirrTable),
                isinstance(light_airr_table, AirrTable),
            ]
        ):
            raise TypeError(
                f"heavy_airr_table {heavy_airr_table} and light_airr_table {light_airr_table} must be instance of {AirrTable}"
            )

        self.join_on = list(join_on)
        self._heavy = heavy_airr_table.set_index(self.join_on)
        self._light = light_airr_table.set_index(self.join_on)
        if sorted(self._heavy.index) != sorted(self._light.index):
            logger.warning("heavy airr and light airr don't share the same indexes, there maybe unpaired chains")
            # raise JoinAirrError(
            #     "heavy airr and light airr must have same sequence identification indexes so we can figure out pairing"
            # )
        if len(self._heavy) != len(self._light):
            logger.warning("heavy airr and light airr don't share the same indexes, there maybe unpaired chains")
            # raise JoinAirrError("heavy airr and light airr must be the same lenght")

        self._joined_table = self._heavy.join(
            self._light, rsuffix="_light", lsuffix="_heavy", how="outer"
        ).reset_index()
        self._joined_table_inner = self._heavy.join(
            self._light, rsuffix="_light", lsuffix="_heavy", how="inner"
        ).reset_index()

    @property
    def table(self) -> pd.DataFrame:
        return self._joined_table

    @property
    def table_inner(self) -> pd.DataFrame:
        return self._joined_table

    @property
    def heavy(self) -> pd.DataFrame:
        return self._heavy.reset_index()

    @property
    def light(self) -> pd.DataFrame:
        return self._light.reset_index()

    @property
    def iloc(self) -> pd.core.indexing._iLocIndexer:
        return self._joined_table.iloc

    @property
    def loc(self) -> pd.core.indexing._LocIndexer:
        return self._joined_table.loc

    def set_index(self, *args, **kwargs):
        return self._joined_table.set_index(*args, **kwargs)

    def get_sanitized_antibodies(self):
        """
        Get only productive full length scfv reads with heavy and light chains
        """
        _heavy_clean = AirrTable(self.heavy).sanitized_antibodies
        _light_clean = AirrTable(self.light).sanitized_antibodies
        _heavy_clean = _heavy_clean.set_index(self.join_on)
        _light_clean = _light_clean.set_index(self.join_on)
        return _heavy_clean.join(_light_clean, rsuffix="_light", lsuffix="_heavy", how="inner").reset_index()

    def get_genbank(self):
        _genbanks = []
        # I'm no sure there is a better way to do this other than go through one by one
        for row in self._joined_table.iterrows():
            _heavy_record = AirrTable.parse_row_to_genbank(row, suffix="_heavy")
            _light_record = AirrTable.parse_row_to_genbank(row, suffix="_light")
            _heavy_record.features += _light_record.features
            _genbanks.append(_heavy_record)
        return _genbanks

    def to_genbank(self, filename: str, compression="gzip"):
        """write airrtable to genbank file

        Parameters
        ----------
        filename : str
            file name string path
        compression: ['gzip']
        """
        _records = self.get_genbank()
        SeqIO.write(_records, filename, "genbank")

    def get_json(self, indent=False) -> json:
        """Get json string of HeavyLightAirrTable

        Parameters
        ----------
        indent : bool, optional
            indent json string representation, by default False

        Returns
        -------
        json
            HeavyLightAirrTable as a json table
        """
        return self._joined_table.to_json(indent=indent, orient="records")

    def to_json(self, filename: str, compression="gzip"):
        """write airrtable to json file

        Parameters
        ----------
        filename : str
            file name string path
        """
        self._joined_table.to_json(filename, orient="records", compression=compression)

    def to_csv(self, path_or_buf=None, compression="infer", *args, **kwargs):
        """write airrtable to csv

        Overloaded pandas.to_csv(*args, **kwargs)
        """
        self._joined_table.to_csv(path_or_buf, **kwargs)

    @staticmethod
    def from_json(json_object: json) -> "ScfvAirrTable":
        """take in json object serialized and ScfvAirrTable

        Returns
        -------
        ScFVAirrTable
            ScFVAirrTable object"""

        df = pd.read_json(json_object)
        joined_heavy = df[[i for i in df.columns if "_heavy" in i] + ["sequence_id"]].rename(
            columns=lambda x: x.replace("_heavy", "")
        )
        joined_light = df[[i for i in df.columns if "_light" in i] + ["sequence_id"]].rename(
            columns=lambda x: x.replace("_light", "")
        )

        return ScfvAirrTable(AirrTable(joined_heavy), AirrTable(joined_light))

    @staticmethod
    def deconstruct_scfv(dataframe: Union["ScfvAirrTable", pd.DataFrame]) -> Tuple[AirrTable, AirrTable]:
        """A static method to deconstruct an ScFV airr table into a heavy and light airr table

        Parameters
        ----------
        dataframe : Either an ScFVAirrTable or a pandas dataframe

        Returns
        -------
        Tuple[AirrTable, AirrTable]
            Heavy airr and light airr
        """
        non_airr = []
        if isinstance(dataframe, ScfvAirrTable):
            dataframe = dataframe.table
        for x in dataframe.columns:
            if "_heavy" in x or "_light" in x:
                continue
            non_airr.append(x)

        logger.debug(f"non-airr columns - {non_airr}")
        joined_heavy = dataframe[non_airr + [i for i in dataframe.columns if "_heavy" in i]].rename(
            columns=lambda x: x.replace("_heavy", "")
        )

        joined_light = dataframe[non_airr + [i for i in dataframe.columns if "_light" in i]].rename(
            columns=lambda x: x.replace("_light", "")
        )
        return AirrTable(joined_heavy), AirrTable(joined_light)

    @staticmethod
    def read_csv(file: str) -> "ScfvAirrTable":
        """Read a csv file and return an airrtable object

        Returns
        -------
        AirrTable
            AirrTable object"""
        _heavy, _light = ScfvAirrTable.deconstruct_scfv(pd.read_csv(file))
        return ScfvAirrTable(_heavy, _light)

    def __getitem__(self, col):
        return self._joined_table[col]

    def __len__(self):
        return len(self._joined_table)

    def __repr__(self):
        return self._joined_table.__repr__()

    def _repr_html_(self):
        """HTML reprsentation for jupyter notebooks"""
        return self._joined_table._repr_html_()

    def __setitem__(self, key, val):
        self._joined_table[key] = val


if __name__ == "__main__":
    pass
