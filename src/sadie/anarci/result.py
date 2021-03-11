# import gzip
# import bz2
# import json
# from pathlib import Path
# from typing import List, Union
import logging
import pandas as pd
from ast import literal_eval

# from ..antibody import chain
from .constants import ANARCI_RESULTS

# from .numbering import scheme_numbering

logger = logging.getLogger("ANARCI")


class AnarciResults(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        # use the __init__ method from DataFrame to ensure
        # that we're inheriting the correct behavior
        super(AnarciResults, self).__init__(*args, **kwargs)

    # this method is makes it so our methods return an instance
    # of ExtendedDataFrame, instead of a regular DataFrame
    @property
    def _constructor(self):
        return AnarciResults

    @staticmethod
    def read_csv(*args, **kwargs):
        return AnarciResults(
            pd.read_csv(
                *args,
                index_col=0,
                dtype=ANARCI_RESULTS,
                converters={"Numbering": literal_eval, "Insertion": literal_eval, "Numbered_Sequence": literal_eval},
                **kwargs
            )
        )

    # now define a custom method!
    # note that `self` is a DataFrame
    # def select_by_grade_and_teacher(self, grade, teacher):
    #     return self.query("grade == @grade").query("teacher == @teacher")

    # def __init__(self, summary_df: pd.DataFrame, scheme: str):
    #     self.main_df = summary_df
    #     self.scheme = scheme

    # @property
    # def main_df(self) -> pd.DataFrame:
    #     return self._main_df

    # @main_df.setter
    # def main_df(self, main_df: pd.DataFrame):
    #     if not isinstance(main_df, pd.DataFrame):
    #         raise TypeError(f"{main_df} must be of type {pd.DataFrame}")
    #     self._main_df = main_df

    # def __getattr__(self, item) -> property:
    #     """Get attribute if attribute is not found in self. For this class, its
    #     a pandas dataframe df attributes is looked up

    #     Parameters
    #     ----------
    #     item : str
    #         a pandas dataframe attribute

    #     Returns
    #     -------
    #     property
    #         property attribute of a dataframe
    #     """
    #     return self.main_df.__getattribute__(item)


# class AnarciResult:
#     """
#     Object class for a single AnarciResult

#     """

#     def __init__(self, summary: pd.DataFrame):
#         self.summary = summary

#     @property
#     def id(self) -> str:
#         """The main ID of the anarchi result derived from seq.id or manual ID

#         Returns
#         -------
#         str
#             the main id
#         """
#         return str(self.summary["Id"])

#     @property
#     def sequence(self) -> str:
#         return self.summary["sequence"]

#     @property
#     def start(self) -> int:
#         return int(self.summary["seqstart_index"])

#     @property
#     def end(self) -> int:
#         return int(self.summary["seqend_index"])

#     @property
#     def v_gene(self) -> str:
#         return self.summary["v_gene"]

#     @property
#     def j_gene(self) -> str:
#         return self.summary["j_gene"]

#     @property
#     def v_gene_identity(self) -> float:
#         return float(self.summary["v_identity"])

#     @property
#     def j_gene_identity(self) -> float:
#         return float(self.summary["j_identity"])

#     @property
#     def chain(self) -> str:
#         return self.summary["chain_type"]

#     @property
#     def scheme(self) -> str:
#         return self.summary["scheme"]

#     @property
#     def species(self) -> str:
#         return self.summary["identity_species"]

#     @property
#     def leader(self) -> str:
#         return self.sequence[0 : self.start]

#     @property
#     def tail(self) -> str:
#         return self.sequence[self.end + 1 :]

#     @property
#     def region_def(self) -> str:
#         return self.summary["region_definition"]

#     @property
#     def numbering_scheme(self) -> dict:
#         if self.chain in ["K", "L"]:
#             _chain = "light"
#         else:
#             _chain = "heavy"
#         return scheme_numbering[self.scheme][_chain][self.region_def]

#     @property
#     def summary(self) -> pd.Series:
#         return self._summary


#     @property
#     def alignment(self) -> pd.DataFrame:
#         return self._alignment

#     @alignment.setter
#     def alignment(self, alignment: pd.DataFrame):
#         if not isinstance(alignment, pd.DataFrame):
#             raise TypeError(f"{alignment} must be of type {pd.DataFrame}")
#         _ids = alignment["Id"].unique()
#         if len(_ids) > 1:
#             raise LookupError(f"summary dataframe has more than one Id {_ids}")
#         self._alignment = alignment

#     @property
#     def framework1_aa(self):
#         return self._lookup_segment("fwr1_aa_start", "fwr1_aa_end")

#     @property
#     def cdr1_aa(self):
#         return self._lookup_segment("cdr1_aa_start", "cdr1_aa_end")

#     @property
#     def framework2_aa(self):
#         return self._lookup_segment("fwr2_aa_start", "fwr2_aa_end")

#     @property
#     def cdr2_aa(self):
#         return self._lookup_segment("cdr2_aa_start", "cdr2_aa_end")

#     @property
#     def framework3_aa(self):
#         return self._lookup_segment("fwr3_aa_start", "fwr3_aa_end")

#     @property
#     def cdr3_aa(self):
#         return self._lookup_segment("cdr3_aa_start", "cdr3_aa_end")

#     @property
#     def framework4_aa(self):
#         return self._lookup_segment("fwr4_aa_start", "fwr4_aa_end")

#     @property
#     def vdj(self) -> str:
#         return "".join(list(self.alignment["Sequence"]))

#     @property
#     def segment_table(self):
#         return {
#             "Id": self.id,
#             "leader": self.leader,
#             "fwr1_aa": self.framework1_aa,
#             "cdr1_aa": self.cdr1_aa,
#             "fwr2_aa": self.framework2_aa,
#             "cdr2_aa": self.cdr2_aa,
#             "fwr3_aa": self.framework3_aa,
#             "cdr3_aa": self.cdr3_aa,
#             "fwr4_aa": self.framework4_aa,
#             "vdj": self.vdj,
#             "tail": self.tail,
#             "chain": self.chain,
#             "species": self.species,
#             "scheme": self.scheme,
#             "region_definition": self.region_def,
#             "v_gene": self.v_gene,
#             "v_gene_identity": self.v_gene_identity,
#             "j_gene": self.j_gene,
#             "j_gene_identity": self.j_gene_identity,
#         }

#     @property
#     def segment_table_no_gaps(self):
#         updated = {}
#         for key in self.segment_table:
#             value = self.segment_table[key]
#             if key in [
#                 "fwr1_aa",
#                 "cdr1_aa",
#                 "fwr2_aa",
#                 "cdr2_aa",
#                 "fwr3_aa",
#                 "cdr3_aa",
#                 "fwr4_aa",
#                 "vdj",
#             ]:
#                 updated[key] = value.replace("-", "")
#             else:
#                 updated[key] = value
#         return updated

#     def get_segment_table(self, remove_gaps=False):
#         if remove_gaps:
#             return self.segment_table_no_gaps
#         else:
#             return self.segment_table

#     def to_antibody(self):
#         _constructor = ""
#         if self.chain == "H":
#             _constructor = chain.HeavyChainAA
#         elif self.chain == "K":
#             _constructor = chain.KappaChainAA
#         elif self.chain == "L":
#             _constructor = chain.LambdaChainAA
#         else:
#             raise TypeError(f"{self.chain} is not H,K, or L.")
#         return _constructor(
#             name=self.id,
#             fwr1_aa=self.framework1_aa.replace("-", ""),
#             cdr1_aa=self.cdr1_aa.replace("-", ""),
#             fwr2_aa=self.framework2_aa.replace("-", ""),
#             cdr2_aa=self.cdr2_aa.replace("-", ""),
#             fwr3_aa=self.framework3_aa.replace("-", ""),
#             cdr3_aa=self.cdr3_aa.replace("-", ""),
#             fwr4_aa=self.framework4_aa.replace("-", ""),
#             v_gene=self.v_gene,
#             j_gene=self.j_gene,
#             species=self.species,
#             leader=self.leader,
#             tail=self.tail,
#             region_def=self.region_def,
#         )

#     def _lookup_segment(self, start_key, end_key):
#         return "".join(
#             list(
#                 self.alignment[
#                     (self.alignment["Numbering"] >= self.numbering_scheme[start_key])
#                     & (self.alignment["Numbering"] <= self.numbering_scheme[end_key])
#                 ]["Sequence"]
#             )
#         )


#     def get_json(self, indent=0):
#         return json.dumps(self.serialize, indent=indent)

#     def to_json(self, file: Path, compression="gzip"):
#         file = Path(file)
#         if compression == "gzip":
#             file = gzip.open(file, "wt")
#         json.dump(self.serialize, file)

#     @staticmethod
#     def from_json(from_json: json):
#         deserialize = json.loads(from_json)
#         return AnarciResult(pd.Series(deserialize["summary"]), pd.DataFrame(deserialize["alignment"]))

#     @staticmethod
#     def read_json(file: Path):
#         _filetype = filetype.guess(file)
#         if _filetype.extension == "gz":
#             file_buffer = gzip.open(file)
#         elif _filetype.extension == "bz2":
#             file_buffer = bz2.open(file)
#         else:
#             file_buffer = open(file)

#         deserialize = json.load(file_buffer)
#         return AnarciResult(pd.Series(deserialize["summary"]), pd.DataFrame(deserialize["alignment"]))

#     def __eq__(self, other):
#         _summary = (self.summary == other.summary).all()
#         _alignment = (self.alignment == other.alignment).all().all()
#         return _summary and _alignment

#     def __repr__(self):
#         return self.get_segment_table(True).__repr__()

#     def __str__(self):
#         return self.__repr__().__str__()


# class AnarciResults:
#     def __init__(self, **kwargs):
#         """Anarci Results

#         Can either construct with list of AnarciResult

#         of conrstruct with a summary,alignment and segment dataframe

#         Paramaters
#         ----------
#         results = [AnarciResult]

#         or

#         summary_table = [summary_table]
#             summarizes anarci aligments
#         alignemtn_table = alignment_table dataframe
#             aligment table of all anarci results
#         segment_table = segment_table dataframe
#             segment table from anarci, eg. CDR1

#         Raises
#         ------
#         ValueError
#             wrong paramamets
#         """
#         if "results" in kwargs.keys():
#             results = kwargs.get("results")
#             self._summary_table = pd.DataFrame([i.summary.to_dict() for i in results]).set_index("Id")
#             self._alignment_table = pd.concat([i.alignment for i in results]).set_index("Id")
#             self._segment_table = pd.DataFrame([i.segment_table for i in results]).set_index("Id")
#             # self.results = results
#         elif all([i in kwargs.keys() for i in ["summary_table", "alignment_table", "segment_table"]]):
#             self.summary_table = kwargs.get("summary_table")
#             self.alignment_table = kwargs.get("alignment_table")
#             self.segment_table = kwargs.get("segment_table")
#         else:
#             raise ValueError(
#                 f"Anarci results must contain either list of Anarci Results or must pass summary_table, alignment_table, segment_table. Passed {kwargs}"
#             )

#     @property
#     def summary_table(self) -> pd.DataFrame:
#         return self._summary_table

#     @summary_table.setter
#     def summary_table(self, summary_table_df: pd.DataFrame):
#         cols = [
#             "sequence",
#             "domain_no",
#             "hmm_species",
#             "chain_type",
#             "e-value",
#             "score",
#             "seqstart_index",
#             "seqend_index",
#             "identity_species",
#             "v_gene",
#             "v_identity",
#             "j_gene",
#             "j_identity",
#             "Chain",
#             "scheme",
#             "allowed_species",
#             "allowed_chains",
#         ]
#         if not isinstance(summary_table_df, pd.DataFrame):
#             TypeError(f"{summary_table_df} must be of type {type(pd.DataFrame)}")
#         if summary_table_df.index.name != "Id":
#             TypeError(f"Index of summary table is {summary_table_df.index.name} must be Id")
#         if [i in summary_table_df.columns for i in cols]:
#             TypeError(f"summary table must contain {cols}, missing {set(cols).difference(summary_table_df.columns)}")
#         self._summary_table = summary_table_df

#     @property
#     def segment_table(self) -> pd.DataFrame:
#         return self._segment_table

#     @segment_table.setter
#     def segment_table(self, segment_table_df: pd.DataFrame):
#         cols = [
#             "leader",
#             "fwr1_aa",
#             "cdr1_aa",
#             "fwr2_aa",
#             "cdr2_aa",
#             "fwr3_aa",
#             "cdr3_aa",
#             "fwr4_aa",
#             "vdj",
#             "tail",
#             "chain",
#             "species",
#             "scheme",
#             "v_gene",
#             "v_gene_identity",
#             "j_gene",
#             "j_gene_identity",
#         ]
#         if not isinstance(segment_table_df, pd.DataFrame):
#             TypeError(f"{segment_table_df} must be of type {type(pd.DataFrame)}")
#         if segment_table_df.index.name != "Id":
#             TypeError(f"Index of segment table is {segment_table_df.index.name} must be Id")
#         if [i in segment_table_df.columns for i in cols]:
#             TypeError(f"segment table must contain {cols}, missing {set(cols).difference(segment_table_df.columns)}")
#         self._segment_table = segment_table_df

#     @property
#     def alignment_table(self) -> pd.DataFrame:
#         return self._alignment_table

#     @alignment_table.setter
#     def alignment_table(self, alignment_table_df: pd.DataFrame):
#         cols = ["Chain", "Numbering", "Insertion", "Sequence", "scheme"]
#         if not isinstance(alignment_table_df, pd.DataFrame):
#             TypeError(f"{alignment_table_df} must be of type {type(pd.DataFrame)}")
#         if alignment_table_df.index.name != "Id":
#             TypeError(f"Index of alignment table is {alignment_table_df.index.name} must be Id")
#         if [i in alignment_table_df.columns for i in cols]:
#             TypeError(
#                 f"alignment table must contain {cols}, missing {set(cols).difference(alignment_table_df.columns)}"
#             )
#         self._alignment_table = alignment_table_df

#     @property
#     def segment_table_no_gaps(self) -> pd.DataFrame:
#         update_cols = [
#             "fwr1_aa",
#             "cdr1_aa",
#             "fwr2_aa",
#             "cdr2_aa",
#             "fwr3_aa",
#             "cdr3_aa",
#             "fwr4_aa",
#             "vdj",
#         ]
#         column_order = self.segment_table.columns
#         return pd.concat(
#             [
#                 self.segment_table[update_cols].apply(lambda x: x.str.replace("-", ""), axis=1),
#                 self.segment_table[column_order.difference(update_cols)],
#             ],
#             axis=1,
#         )[column_order]

#     def get_segment_table(self, remove_gaps=False) -> pd.DataFrame:
#         """Return a dataframe of the segment table with or without the IMGT gaps

#         Parameters
#         ----------
#         remove_gaps : bool, optional
#             remove the imgt gaps if present, by default False

#         Returns
#         -------
#         DataFrame
#             pandas dataframe of the segment table
#         """
#         if remove_gaps:
#             return self.segment_table_no_gaps
#         return self.segment_table

#     def get_result(self, id: Union[str, int]) -> AnarciResult:
#         """Get a single anarciresult

#         Parameters
#         ----------
#         id : Union[str, int]
#            The id of the result to return

#         Returns
#         -------
#         AnarciResult
#             Single Anarci result
#         """
#         _alignment = self.alignment_table.loc[id]
#         _summary = self.summary_table.loc[id]
#         return AnarciResult(_summary, _alignment.reset_index())

#     @property
#     def serialize(self) -> dict:
#         """Serialized Object for Anarci Results

#         Returns
#         -------
#         dict
#            dictionary of dictionaries serialized
#         """
#         return {
#             "segment_table": self.segment_table.reset_index().to_dict(),
#             "summary_table": self.summary_table.reset_index().to_dict(),
#             "alignment_table": self.alignment_table.reset_index().to_dict(),
#         }

#     def get_json(self, indent=0) -> json:
#         """Return json string representation of Anarci results

#         Parameters
#         ----------
#         indent : int, optional
#             if you want json indented, by default 0

#         Returns
#         -------
#         json
#            serializable string representation
#         """

#         return json.dumps(self.serialize, indent=indent)

#     @staticmethod
#     def from_json(json_str: str) -> "AnarciResults":
#         """Read json object and return AnarciResults

#         Returns
#         -------
#         AnarciResults
#             new AnarciResults method from json string
#         """
#         deserialize = json.loads(json_str)
#         return AnarciResults(
#             segment_table=pd.DataFrame.from_dict(deserialize["segment_table"]).set_index("Id"),
#             alignment_table=pd.DataFrame.from_dict(deserialize["alignment_table"]).set_index("Id"),
#             summary_table=pd.DataFrame.from_dict(deserialize["summary_table"]).set_index("Id"),
#         )

#     @staticmethod
#     def read_file(file: Path) -> "AnarciResults":
#         file = Path(file)
#         if not file.exists():
#             raise FileNotFoundError(f"can't find {file}")
#         _filetype = filetype.guess(str(file.absolute()))
#         if _filetype.extension == "gz":
#             file_buffer = gzip.open(file)
#         elif _filetype.extension == "bz2":
#             file_buffer = bz2.open(file)
#         else:
#             file_buffer = open(file)
#         deserialize = pickle.load(file_buffer)
#         return AnarciResults(
#             segment_table=deserialize["segment_table"],
#             alignment_table=deserialize["alignment_table"],
#             summary_table=deserialize["summary_table"],
#         )

#     def to_file(self, file: Path, compression="bzip"):
#         """Write out AnarciResults to a file
#         Parameters
#         ----------
#         file : Path
#             file path
#         """
#         if compression == "gzip":
#             file = Path(file + ".gz")
#             file_buffer = gzip.open(file, "wt")
#         elif compression == "bzip":
#             if not file.endswith(".bz2"):
#                 file = Path(file + ".bz2")
#             else:
#                 file = Path(file)
#             file_buffer = bz2.open(file, "wb")
#         else:
#             file_buffer = open(file, "wb")
#         with file_buffer:
#             pickle.dump(
#                 {
#                     "segment_table": self.segment_table,
#                     "alignment_table": self.alignment_table,
#                     "summary_table": self.summary_table,
#                 },
#                 file_buffer,
#             )
#         logger.info(f"Wrote file to {file.absolute()}")

#     def to_json(self, file: Path, compression="gzip"):
#         """Write out AnarciResults to a json file

#         Parameters
#         ----------
#         file : Path
#             file path
#         compression : str, optional
#             compress file object, by default "gzip"
#         """
#         file = Path(file)
#         if compression == "gzip":
#             gzip_file = gzip.open(file, "wt")
#             json.dump(self.serialize, gzip_file)
#         else:
#             json.dump(self.serialize, file)

#         logger.info(f"Wrote json file to {file.absolute()}")

#     @staticmethod
#     def read_json(file: Path):
#         _filetype = filetype.guess(file)
#         if _filetype.extension == "gz":
#             file_buffer = gzip.open(file)
#         elif _filetype.extension == "bz2":
#             file_buffer = bz2.open(file)
#         else:
#             file_buffer = open(file)
#         deserialize = json.load(file_buffer)
#         return AnarciResults(
#             segment_table=pd.DataFrame.from_dict(deserialize["segment_table"]).set_index("Id"),
#             alignment_table=pd.DataFrame.from_dict(deserialize["alignment_table"]).set_index("Id"),
#             summary_table=pd.DataFrame.from_dict(deserialize["summary_table"]).set_index("Id"),
#         )

#     @staticmethod
#     def concat(results: List["AnarciResults"]) -> "AnarciResults":
#         _summary = pd.concat([i.summary_table for i in results])
#         _alignment = pd.concat([i.alignment_table for i in results])
#         _segment = pd.concat([i.segment_table for i in results])
#         return AnarciResults(segment_table=_segment, alignment_table=_alignment, summary_table=_summary)

#     def __add__(self, other):
#         self._summary_table = pd.concat([self._summary_table, other._summary_table])
#         self._alignment_table = pd.concat([self._alignment_table, other._alignment_table])
#         self._segment_table = pd.concat([self._segment_table, other._segment_table])
#         return self

#     def __radd__(self, other):
#         if other == 0:
#             return self
#         else:
#             return self.__add__(other)

#     def __iter__(self):
#         for r in self.results:
#             yield r

#     def __repr__(self):
#         return self._segment_table.__repr__()

#     def _repr_html_(self):
#         """HTML reprsentation for jupyter notebooks"""
#         return self._segment_table._repr_html_()

#     def __str__(self):
#         return self.__repr__().__str__()

#     def __getitem__(self, id):
#         return self.get_result(id)

#     def __eq__(self, other):
#         return all(
#             [
#                 (self.alignment_table.sort_index() == other.alignment_table.sort_index()).all().all(),
#                 (self.summary_table.sort_index() == other.summary_table.sort_index()).all().all(),
#                 (self.segment_table.sort_index() == other.segment_table.sort_index()).all().all(),
#             ]
#         )
