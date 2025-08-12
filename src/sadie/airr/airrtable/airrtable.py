"""The AirrTable module"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from pprint import pformat
from typing import Any, Dict, Hashable, List, Optional, Tuple, Type, Union

import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
from Levenshtein import distance  # type: ignore
from numpy import nan

from sadie.airr.airrtable.constants import CONSTANTS_AIRR, IGBLAST_AIRR, OTHER_COLS
from sadie.airr.airrtable.genbank import GenBank, GenBankFeature
from sadie.airr.exceptions import MissingAirrColumns
from sadie.airr.models import AirrSeriesModel
from sadie.receptor.rearrangment import (
    AlignmentAnnotations,
    AlignmentPositions,
    InputSequence,
    JunctionLengths,
    PrimaryAnnotations,
    ReceptorChain,
    RegionPositions,
    RegionSequences,
)
from sadie.utility.io import SadieOutput
from sadie.utility.util import correct_alignment

logger = logging.getLogger("AIRRTable")


# DEPRECARED
# def get_reverse_compliment(seq: Seq) -> Seq:
#     return Seq(seq.reverse_complement())


# def _get_raw_seq(row: pd.Series) -> str: # type:
#     seq: Seq = Seq(row["sequence"])
#     rev_comp: bool = row["rev_comp"]
#     if rev_comp:
#         seq = get_reverse_compliment(seq)
#     return str(seq)


# def _get_seq_aa(seq: str) -> str:
#     return str(Seq(seq).translate())


def is_all_integer(row_1: int | float, row_2: int | float) -> bool:
    """Checks if all values in a series are integers"""
    if all([isinstance(row_1, int), isinstance(row_2, int)]):
        return True
    return False


class AirrSeries(pd.Series):  # type: ignore
    _metadata = ["meta"]  # add custom namespaces here

    def __init__(self, data: Any, copy: bool = False, *args, **kwargs):
        super(AirrSeries, self).__init__(data=data, copy=copy, *args, **kwargs)  # type: ignore
        # Only run verification if data is not from internal operations (SingleBlockManager)
        if not (data.__class__.__name__ == "SingleBlockManager"):
            if isinstance(data, pd.Series):
                self._verify()

    @property
    def _constructor(self) -> Type["AirrSeries"]:
        return AirrSeries

    @property
    def _constructor_expanddim(self) -> Type["AirrTable"]:
        return AirrTable

    def _verify(self) -> None:
        """Verifies that the AirrSeries is valid"""
        data: Dict[Any, Any] = AirrSeriesModel(**self).dict()  # type: ignore
        self.update(data)

    def to_input_sequence_object(self) -> InputSequence:
        """
        Converts the AirrSeries to a InputSequence
        """
        sequence_fields: pd.Index = pd.Index(InputSequence.get_airr_fields())
        model = InputSequence(**self[sequence_fields].to_dict())
        return model

    def to_primary_annotations_object(self) -> PrimaryAnnotations:
        """
        Converts the AirrSeries to a PrimaryAnnotations
        """
        sequence_fields: pd.Index = pd.Index(PrimaryAnnotations.get_airr_fields())
        model = PrimaryAnnotations(**self[sequence_fields].to_dict())
        return model

    def to_alignment_annotations_object(self) -> AlignmentAnnotations:
        """
        Converts the AirrSeries to a AlignmentAnnotations
        """
        sequence_fields: pd.Index = pd.Index(AlignmentAnnotations.get_airr_fields())
        model = AlignmentAnnotations(**self[sequence_fields].to_dict())
        return model

    def to_region_sequences_object(self) -> RegionSequences:
        """
        Converts the AirrSeries to a RegionSequences
        """
        sequence_fields: pd.Index = pd.Index(RegionSequences.get_airr_fields())
        model = RegionSequences(**self[sequence_fields].to_dict())
        return model

    def to_alignment_positions_object(self) -> AlignmentPositions:
        """
        Converts the AirrSeries to a AlignmentPositions
        """
        sequence_fields: pd.Index = pd.Index(AlignmentPositions.get_airr_fields())
        model = AlignmentPositions(**self[sequence_fields].to_dict())
        return model

    def to_region_positions_object(self) -> RegionPositions:
        """
        Converts the AirrSeries to a RegionPositions
        """
        sequence_fields: pd.Index = pd.Index(RegionPositions.get_airr_fields())
        model = RegionPositions(**self[sequence_fields].to_dict())
        return model

    def to_junction_lengths_object(self) -> JunctionLengths:
        """
        Converts the AirrSeries to a JunctionLengths
        """
        sequence_fields: pd.Index = pd.Index(JunctionLengths.get_airr_fields())
        model = JunctionLengths(**self[sequence_fields].to_dict())
        return model

    def to_receptor_chain_object(self) -> ReceptorChain:
        """
        Converts the AirrSeries to a ReceptorChain
        """
        model = ReceptorChain(
            input_sequence=self.to_input_sequence_object(),
            primary_annotations=self.to_primary_annotations_object(),
            alignment_annotations=self.to_alignment_annotations_object(),
            region_sequences=self.to_region_sequences_object(),
            alignment_positions=self.to_alignment_positions_object(),
            region_positions=self.to_region_positions_object(),
            junction_lengths=self.to_junction_lengths_object(),
        )
        return model


class AirrTable(pd.DataFrame):
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

    >>> at.columns
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
    >>> at = AirrTable.read_airr("tests/fixtures/kappa_lev_sample.tsv.gz")
    """

    # if you create a new property make sur eyou add it here or else it will be lost
    _metadata = ["_suffixes", "_islinked", "_key_column", "_verified"]

    def __init__(self, data: Any = None, key_column: str = "sequence_id", copy: bool = False):
        super(AirrTable, self).__init__(data=data, copy=copy)  # type: ignore
        # Only run initialization if data is not a BlockManager (internal pandas object)
        # BlockManager is passed during internal pandas operations like copy(), reindex(), etc.
        if not (data.__class__.__name__ == "BlockManager"):
            if self.__class__ is AirrTable:
                self._islinked: bool = False
                self._key_column: str = key_column
                self._suffixes: List[str] = []
                self._verify()
                key: str
                for key, value in IGBLAST_AIRR.items():
                    self[key] = self[key].astype(value)  # type: ignore[arg-type]

                # what about the non Airr columns?
                _have_other = self.columns.intersection(list(OTHER_COLS.keys()))
                for key in _have_other:
                    self[key] = self[key].astype(OTHER_COLS[key])  # type: ignore[arg-type]

                # what about the constant Airr columns?
                _have_constant = self.columns.intersection(list(CONSTANTS_AIRR.keys()))
                for key in _have_constant:
                    self[key] = self[key].astype(CONSTANTS_AIRR[key])  # type: ignore[arg-type]

    @property
    def _constructor_sliced(self) -> Type[AirrSeries]:
        return AirrSeries

    @property
    def verified(self) -> bool:
        if not hasattr(self, "_verified"):
            return False
        return self._verified

    @property
    def _constructor(self) -> "AirrTable":
        return AirrTable  # type: ignore[return-value]

    @property
    def key_column(self) -> str:
        return self._key_column

    @property
    def compliant_cols(self) -> List[str]:
        compliant_cols = list(set(IGBLAST_AIRR.keys()))
        if self.key_column not in compliant_cols:
            compliant_cols += [self.key_column]
        return compliant_cols

    @property
    def non_airr_columns(self) -> pd.Index:
        """Get an column index of non airr complient columns.

        These are columns that were probably joined on the airr dataframe, or can be just
        extra fields added to the airr dataframe


        Returns
        -------
        pd.Index
            The non-airr complient columns
        """
        return pd.Index(list(set(self.columns) - set(self.compliant_cols)))

    @property
    def airr_columns(self) -> pd.Index:
        """Get a column index of airr complient columns.

        Returns
        -------
        pd.Index
            airr complient column indexes
        """
        return pd.Index(list(set(self.columns).intersection(self.compliant_cols)))

    def get_sanitized_antibodies(self) -> "AirrTable":
        """Getter method for airr entires that have  all antibody segments

        Returns
        -------
        AirrTable
            AirrTable for just sanitized_antibodies
        """
        antibody_segments = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
        return AirrTable(self[self[antibody_segments].notna().all(axis=1)].copy())

    def to_fasta(
        self,
        filename: Path,
        sequence_field: str = "sequence_alignment",
    ) -> None:
        """given an id field and sequence field, write out dataframe to fasta

        Parameters
        ----------
        file_out : Path
            the file output path to fasta
        sequence_field : str (default: sequence_alignment)
            The seq field to use as the sqeuence for the fasta file
            Options:
                sequence
                sequence_alignment (default)
                germline_alignment
                sequence_alignment_aa
                germline_alignment_aa
                v_sequence_alignment
                v_sequence_alignment_aa
                v_germline_alignment
                v_germline_alignment_aa
                d_sequence_alignment
                d_sequence_alignment_aa
                d_germline_alignment
                d_germline_alignment_aa
                j_sequence_alignment
                j_sequence_alignment_aa
                j_germline_alignment
                j_germline_alignment_aa
                c_sequence_alignment
                c_sequence_alignment_aa
                c_germline_alignment
                c_germline_alignment_aa
                fwr1
                fwr1_aa
                cdr1
                cdr1_aa
                fwr2
                fwr2_aa
                cdr2
                cdr2_aa
                fwr3
                fwr3_aa
                fwr4
                fwr4_aa
                cdr3
                cdr3_aa
                junction
                junction_aa
                np1
                np2
                vdj_nt
                vdj_aa
        """
        seq_fields = [
            "sequence",
            "sequence_alignment",
            "germline_alignment",
            "sequence_alignment_aa",
            "germline_alignment_aa",
            "v_sequence_alignment",
            "v_sequence_alignment_aa",
            "v_germline_alignment",
            "v_germline_alignment_aa",
            "d_sequence_alignment",
            "d_sequence_alignment_aa",
            "d_germline_alignment",
            "d_germline_alignment_aa",
            "j_sequence_alignment",
            "j_sequence_alignment_aa",
            "j_germline_alignment",
            "j_germline_alignment_aa",
            "c_sequence_alignment",
            "c_sequence_alignment_aa",
            "c_germline_alignment",
            "c_germline_alignment_aa",
            "fwr1",
            "fwr1_aa",
            "cdr1",
            "cdr1_aa",
            "fwr2",
            "fwr2_aa",
            "cdr2",
            "cdr2_aa",
            "fwr3",
            "fwr3_aa",
            "fwr4",
            "fwr4_aa",
            "cdr3",
            "cdr3_aa",
            "junction",
            "junction_aa",
            "np1",
            "np2",
            "vdj_nt",
            "vdj_aa",
        ]
        records: List[SeqRecord.SeqRecord] = []
        description_fields = [x for x in self.columns if x not in seq_fields]
        if sequence_field not in seq_fields:
            raise ValueError(f"sequence_field must be one of {seq_fields}")
        with open(filename, "w") as f:
            for _, row in self.iterrows():
                records.append(
                    SeqRecord.SeqRecord(
                        Seq(row[sequence_field]),
                        id=row["sequence_id"],
                        name=row["sequence_id"],
                        description=", ".join([f"{x} {row[x]}" for x in description_fields]),
                        dbxrefs=None,
                        features=None,
                        annotations=None,
                        letter_annotations=None,
                    )
                )
            SeqIO.write(records, f, "fasta")

    def get_genbank(self) -> List[SeqRecord.SeqRecord]:
        """get genbank records from airrtable"""
        _genbanks = []
        # go through as an iterator
        for row in self.iterrows():
            _genbanks.append(AirrTable.parse_row_to_genbank(row))
        return _genbanks

    def to_genbank(self, filename: str) -> None:
        """write airrtable to genbank file

        Parameters
        ----------
        filename : str
            file name string path
        """
        _records = self.get_genbank()
        SeqIO.write(_records, filename, "genbank")

    def _verify(self) -> None:
        missing_columns = set(self.compliant_cols).difference(self.columns)
        if missing_columns:
            raise MissingAirrColumns(missing_columns)
        # drop any unnamed columns
        self.drop([i for i in self.columns if "Unnamed" in i], axis=1, inplace=True)

        # set boolean strings to boolelan types
        _to_boolean = ["productive", "vj_in_frame", "stop_codon", "rev_comp", "v_frameshift", "complete_vdj"]
        if self._islinked:
            _to_boolean_link = []
            for x in self._suffixes:
                _to_boolean_link += list(map(lambda y: y + x, _to_boolean))
            _to_boolean = _to_boolean_link
        self._set_boolean(_to_boolean)

        # a check that allows us to see if the junction presentation is wrong
        liability_keys = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
        if self._islinked:
            for suffix in self._suffixes:
                local_suffix = list(map(lambda x: x + suffix, liability_keys))
                self[f"liable{suffix}"] = self[local_suffix].apply(self._check_j_gene_liability, axis=1)
        else:
            self["liable"] = self[["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa"]].apply(
                self._check_j_gene_liability, axis=1
            )

        # get full nt vdj recombination
        vdj_nt_keys = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
        vdj_aa_keys = ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
        if self._islinked:
            for suffix in self._suffixes:
                _local_suffix_nt = list(map(lambda x: x + suffix, vdj_nt_keys))
                _local_suffix_aa = list(map(lambda x: x + suffix, vdj_aa_keys))
                self[f"vdj_nt{suffix}"] = self[_local_suffix_nt].apply(
                    lambda x: "".join([str(i) for i in x if not isinstance(i, float)]), axis=1
                )
                # get full aa vdj recombination
                self[f"vdj_aa{suffix}"] = self[_local_suffix_aa].apply(
                    lambda x: "".join([str(i) for i in x if not isinstance(i, float)]), axis=1
                )
        else:
            # get full nt vdj recombination
            self["vdj_nt"] = self[vdj_nt_keys].apply(
                lambda x: "".join([str(i) for i in x if not isinstance(i, float)]), axis=1
            )
            # get full aa vdj recombination
            self["vdj_aa"] = self[vdj_aa_keys].apply(
                lambda x: "".join([str(i) for i in x if not isinstance(i, float)]), axis=1
            )

        if self._islinked:
            for suffix in self._suffixes:
                # Insert the top call for each vdj call but right next to the N_call column using the insert method
                for call in ["v_call", "d_call", "j_call"]:
                    # pure light chain columns won't have a dcall
                    if call + suffix == "d_call_light":
                        continue

                    new_call = call + f"_top{suffix}"
                    # drop the volumn if it's already there, that helps with backwards compatibility
                    if new_call in self.columns:
                        self.drop(new_call, inplace=True, axis=1)

                    # Insert right next to the X_call airr_columns
                    self.insert(
                        self.columns.get_loc(call + suffix), new_call, self[call + suffix].str.split(",").str.get(0)
                    )

                # get mutation frequency rather than identity
                self.loc[:, f"v_mutation{suffix}"] = self[f"v_identity{suffix}"].apply(lambda x: (1 - x))
                # self.loc[:, f"v_identity{suffix}"] = self[f"v_identity{suffix}"].apply(lambda x: x / 100)

                # then get a percentage for AA by computing levenshtein
                self.loc[:, f"v_mutation_aa{suffix}"] = self[
                    [f"v_sequence_alignment_aa{suffix}", f"v_germline_alignment_aa{suffix}"]
                ].apply(lambda x: self._get_aa_distance(x), axis=1)

                # do the same for D and J gene segment portions
                self.loc[:, f"d_mutation{suffix}"] = self[f"d_identity{suffix}"].apply(
                    lambda x: (1 - x) if x else np.nan
                )
                # self.loc[:, f"d_identity{suffix}"] = self[f"d_identity{suffix}"].apply(lambda x: x / 100)

                self.loc[:, f"d_mutation_aa{suffix}"] = self[
                    [f"d_sequence_alignment_aa{suffix}", f"d_germline_alignment_aa{suffix}"]
                ].apply(lambda x: self._get_aa_distance(x), axis=1)
                self.loc[:, f"j_mutation{suffix}"] = self[f"j_identity{suffix}"].apply(lambda x: (1 - x))
                # self.loc[:, f"j_identity{suffix}"] = self[f"j_identity{suffix}"].apply(lambda x: x / 100)
                self.loc[:, f"j_mutation_aa{suffix}"] = self[
                    [f"j_sequence_alignment_aa{suffix}", f"j_germline_alignment_aa{suffix}"]
                ].apply(lambda x: self._get_aa_distance(x), axis=1)

        else:
            for call in ["v_call", "d_call", "j_call"]:
                # drop the volumn if it's already there, that helps with backwards compatibility
                if f"{call}_top" in self.columns:
                    self.drop(f"{call}_top", inplace=True, axis=1)
                # pure light chain columns won't have a dcall
                if call in self.columns:
                    # Insert right next to the X_call airr_columns
                    self.insert(
                        self.columns.get_loc(call), f"{call}_top", self[call].fillna("").str.split(",").str.get(0)
                    )

            # get mutation frequency rather than identity
            self.loc[:, "v_mutation"] = self["v_identity"].apply(lambda x: (1 - x))

            # then get a percentage for AA by computing levenshtein
            self.loc[:, "v_mutation_aa"] = self[["v_sequence_alignment_aa", "v_germline_alignment_aa"]].apply(
                lambda x: self._get_aa_distance(x), axis=1
            )

            # do the same for D and J gene segment portions
            self.loc[:, "d_mutation"] = self["d_identity"].apply(lambda x: (1 - x) if x else np.nan)
            self.loc[:, "d_mutation_aa"] = self[["d_sequence_alignment_aa", "d_germline_alignment_aa"]].apply(
                lambda x: self._get_aa_distance(x), axis=1
            )
            self.loc[:, "j_mutation"] = self["j_identity"].apply(lambda x: (1 - x))
            # self.loc[:, "j_identity"] = self["j_identity"].apply(lambda x: x / 100)
            self.loc[:, "j_mutation_aa"] = self[["j_sequence_alignment_aa", "j_germline_alignment_aa"]].apply(
                lambda x: self._get_aa_distance(x), axis=1
            )
        self._verified = True

    def _get_aa_distance(self, X: pd.Series) -> float:
        "get character levenshtein distance, will work with '-' on alignments"
        first: str = X[0]
        second: str = X[1]
        if not first or not second or isinstance(first, float) or isinstance(second, float):
            return nan
        d: int = distance(first, second)
        max_len: int = max(len(first), len(second))
        a: float = float(d / max_len)
        return a

    def _set_boolean(self, columns: List[str]) -> None:
        """Change 'F' and 'T' strings to boolean dataframe values"""
        for column in columns:
            self[column] = (
                self[column]
                .map(
                    {"T": True, "F": False, True: True, False: False, nan: False},
                    na_action="ignore",
                )
                .astype(bool)
            )

    def _unset_boolean(self, columns: List[str]) -> pd.DataFrame:
        """Change 'False' and 'True' boolean dataframe values to 'F' and 'T' strings"""
        # don't want to change object level df
        _copy_df = self.copy()
        for column in columns:
            _copy_df[column] = (
                _copy_df[column]
                .map(
                    {True: "T", False: "F"},
                    na_action="ignore",
                )
                .astype(str)
            )
        return _copy_df

    def _check_j_gene_liability(self, X: pd.Series) -> bool:
        """
        Check that the CDR3 was recombined with a FW4 J gene.
        A liable sequence can point to issues with the pentltes

        There are a few different liabilities:
        1. You have CDR1 and CDR2 but dont have CDR3, that can be a sign of a libaility sequence
            - Liability

        2. If we get FW3 + CDR3 but no FW4 we probably have a libility
            - Liability

        3. If we get only get FW3 and CDR3 + FW4, we probably have just sort read sequecning from
            - No liability

        Parameters
        ----------
        X : pd.Series
            ex. [fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4]

        Returns
        -------
        bool
           if the entry is liable or not based
        """
        fw1 = X[0]
        cdr1 = X[1]
        fw2 = X[2]
        cdr2 = X[3]
        fw3 = X[4]
        cdr3 = X[5]
        fw4 = X[6]
        isna = [isinstance(i, float) for i in [fw1, cdr1, fw2, cdr2, fw3, cdr3, fw4]]
        # if they are all nan:
        if all(isna):
            return False

        if sum(isna) == (len(isna) - 1):
            return False

        # If all fw and cdrs are reolsved
        if sum(isna) == 0:
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

        logger.debug(f"Liability shouldn't have gotten here {pformat(X.to_dict())}")
        return True

    def correct_indel(self) -> Union["AirrTable", "LinkedAirrTable"]:
        _fields = [
            ["sequence_alignment_aa", "germline_alignment_aa"],
            ["v_sequence_alignment_aa", "v_germline_alignment_aa"],
        ]
        if self.__class__ is LinkedAirrTable:
            _linked_fields: List[List[str]] = []
            for suffix in self._suffixes:
                for field_1, field_2 in _fields:
                    _linked_fields.append([f"{field_1}{suffix}", f"{field_2}{suffix}"])
            _fields = _linked_fields

        def _get_indel_index(field_1: str, field_2: str) -> pd.Index:
            index: pd.Index = self[
                (self[field_1].str.len() != self[field_2].str.len()) & (~self[field_1].isna()) & (~self[field_2].isna())
            ].index
            return index

        def _update_align(field_1: str, field_2: str) -> None:
            # get indels in sequence_germline_alignment_aa that were not accounted for in the total alignment
            indel_indexes = _get_indel_index(field_1, field_2)
            logger.info(f"Have {len(indel_indexes)} possible indels that are not in amino acid germline alignment")
            self.loc[:, f"{field_2}_corrected"] = False
            if not indel_indexes.empty:
                # Pass self through top level dataframe to remove validation checks when using apply.
                correction_alignments = pd.DataFrame(self.loc[indel_indexes, :]).apply(
                    lambda x: correct_alignment(x, field_1, field_2), axis=1
                )
                # Correction in place
                self.update(correction_alignments)  # type: ignore
            self.loc[indel_indexes, f"{field_2}_corrected"] = True

        for field_1, field_2 in _fields:
            _update_align(field_1, field_2)

        return self

    @staticmethod
    def parse_row_to_genbank(row: Tuple[Optional[Hashable], pd.Series], suffix: str = "") -> SeqRecord.SeqRecord:
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

        # v call could be nan or None
        if not isinstance(v_call, float) and v_call:
            top_v = v_call.split(",")[0]
            other_v = v_call.split(",")[1:]

        if not isinstance(d_call, float) and d_call:
            top_d = d_call.split(",")[0]
            other_d = d_call.split(",")[1:]

        if not isinstance(j_call, float) and j_call:
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
        if "refrence_name" in single_seq.keys():
            species = single_seq["reference_name"]
        else:
            warnings.warn(f"reference is not in {list(single_seq.keys())}")
            species = "Unknown"

        # setup genbank object
        gb = GenBank(sequence, str(single_seq["sequence_id"]), description="AIRR annotation")

        # default qualifier dictionaries
        qd_v = {"gene": top_v, "reference": species, "other_vgene": other_v}
        qd_d = {"gene": top_d, "reference": species, "other_dgene": other_d}
        qd_j = {"gene": top_j, "reference": species, "other_jgene": other_j}

        # FW1
        if is_all_integer(_fw1_start, _fw1_end):
            feature = GenBankFeature(_fw1_start, _fw1_end, "FWR1", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # CDR1
        if is_all_integer(_cdr1_start, _cdr1_end):
            feature = GenBankFeature(_cdr1_start, _cdr1_end, "CDR1", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # FW2
        if is_all_integer(_fw2_start, _fw2_end):
            feature = GenBankFeature(_fw2_start, _fw2_end, "FWR2", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # CDR2
        if is_all_integer(_cdr2_start, _cdr2_end):
            feature = GenBankFeature(_cdr2_start, _cdr2_end, "CDR2", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # FW3
        if is_all_integer(_fw3_start, _fw3_end):
            feature = GenBankFeature(_fw3_start, _fw3_end, "FWR3", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # CDR3
        if is_all_integer(_cdr3_start, _cdr3_end):
            feature = GenBankFeature(_cdr3_start, _cdr3_end, "CDR3")
            gb.add_feature(feature)

        # FW4
        if is_all_integer(_fw4_start, _fw4_end):
            feature = GenBankFeature(_fw4_start, _fw4_end, "FWR4", qualifier_dict=qd_j)
            gb.add_feature(feature)

        # VGene
        if is_all_integer(_v_segment_start, _v_segment_end):
            qd_v = {"gene": top_v, "species": species}
            feature = GenBankFeature(_v_segment_start, _v_segment_end, "VGene", qualifier_dict=qd_v)
            gb.add_feature(feature)

        # DGene
        if is_all_integer(_d_segment_start, _d_segment_end):
            qd_d = {"gene": top_d, "species": species}
            feature = GenBankFeature(_d_segment_start, _d_segment_end, "DGene", qualifier_dict=qd_d)
            gb.add_feature(feature)

        # JGene
        if is_all_integer(_j_segment_start, _j_segment_end):
            qd_j = {"gene": top_j, "species": species}
            feature = GenBankFeature(_j_segment_start, _j_segment_end, "JGene", qualifier_dict=qd_j)
            gb.add_feature(feature)

        # Locus
        if is_all_integer(_v_segment_start, _fw4_end):
            _p = "Productive"
            _vdj = "Complete"
            if not _productive:
                _p = "Non-Productive"
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

        gb.record.annotations["molecule_type"] = "DNA"
        return gb.record

    def to_airr(self, filename: str) -> None:
        """Output AIRR object to AIRR complient file

        Parameters
        ----------
        filename : str
            File name to output to, must end with .tsv or tsv.gz, tsv.bz which will be compressed
        """
        reverse_to_boolean = ["productive", "vj_in_frame", "stop_codon", "rev_comp", "v_frameshift", "complete_vdj"]
        suffixes = Path(filename).suffixes
        if len(suffixes) > 1 and suffixes[0] != ".tsv":
            raise ValueError("File name must end with .tsv or .tsv.gz or .tsv.bz")
        elif suffixes[0] != ".tsv":
            raise ValueError("File name must end with .tsv or .tsv.gz or .tsv.bz")
        self._unset_boolean(reverse_to_boolean).to_csv(filename, sep="\t")

    @staticmethod
    def read_airr(filename: Union[str, Path]) -> "AirrTable":
        """Read AIRR file into AIRR Table object

        Parameters
        ----------
        filename : Union[str, Path]
            File name to read in

        Returns
        -------
        AirrTable
            AIRR Table object
        """
        at = pd.read_csv(filename, sep="\t")
        return AirrTable(at)

    # Deprecated: no internal usage; falling back on built-in pandas __eq__ method
    # def __eq__(self, other: object) -> bool:  # type: ignore[override]
    #     """equals method for AirrTable"""
    #     # needs to be cast to dataframe and na needs to be fileld
    #     if not isinstance(other, AirrTable):
    #         raise NotImplementedError("Can only compare AirrTable to AirrTable")
    #     _dataframe = pd.DataFrame(self).fillna("")
    #     other_dataframe = pd.DataFrame(other).fillna("")
    #     _is_equal: bool = (_dataframe == other_dataframe).all().all()
    #     return _is_equal


class LinkedAirrTable(AirrTable):
    def __init__(
        self,
        data: Any = None,
        suffixes: List[str] = ["_heavy", "_light"],
        key_column: str = "sequence_id",
        copy: bool = False,
    ):
        super(LinkedAirrTable, self).__init__(data=data, copy=copy)
        # Only run initialization if data is not a BlockManager (internal pandas object)
        if not (data.__class__.__name__ == "BlockManager"):
            if isinstance(self, LinkedAirrTable):
                self._islinked = True
                self._key_column = key_column
                self._suffixes = suffixes
                if not self.verified:
                    self._verify()

    @property
    def key_column(self) -> str:
        return self._key_column

    @property
    def suffixes(self) -> List[str]:
        return self._suffixes

    @property
    def left_suffix(self) -> str:
        return self.suffixes[0]

    @property
    def right_suffix(self) -> str:
        return self.suffixes[1]

    @property
    def compliant_cols(self) -> List[str]:
        _complient_cols = []
        joined_keys = list(IGBLAST_AIRR.keys())
        for suffix in self._suffixes:
            _complient_cols += list(map(lambda x: x + suffix if x != f"{self.key_column}" else x, joined_keys))
        _complient_cols = list(set(_complient_cols))
        _complient_cols += [self.key_column]
        return _complient_cols

    def get_split_table(self) -> Tuple[AirrTable, AirrTable]:
        left_rows = [i for i in self.columns if self.left_suffix in i]
        right_rows = [i for i in self.columns if self.right_suffix in i]
        key_column = self.key_column
        common_columns = list(set(self.columns).difference(set(left_rows + right_rows)))
        left_table = self[common_columns + left_rows]
        left_airr_columns: List[str] = list(map(lambda x: x.replace(self._suffixes[0], ""), list(left_table.columns)))  # type: ignore[no-any-return]
        left_table.columns = left_airr_columns
        right_table = self[common_columns + right_rows]
        right_airr_columns = list(map(lambda x: x.replace(self._suffixes[1], ""), list(right_table.columns)))  # type: ignore[no-any-return]
        right_table.columns = right_airr_columns
        return AirrTable(left_table, key_column=key_column), AirrTable(right_table, key_column=key_column)
