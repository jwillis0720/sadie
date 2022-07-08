import logging
import pandas as pd
from ast import literal_eval

from .constants import NUMBERING_RESULTS
from sadie.numbering.scheme_numbering import scheme_numbering

logger = logging.getLogger("NUMBERING")


class NumberingResults(pd.DataFrame):
    def __init__(self, *args, scheme="", region_definition="", allowed_chains=[], allowed_species=[], **kwargs):
        # use the __init__ method from DataFrame to ensure
        # that we're inheriting the correct behavior
        super(NumberingResults, self).__init__(*args, **kwargs)
        # self["scheme"] = scheme
        # self["region_definition"] = region_definition
        # self["allowed_species"] = ",".join(allowed_species)
        # self["allowed_chains"] = ",".join(allowed_chains)
        # self._add_segment_regions()

    @property
    def _constructor(self):
        return NumberingResults

    def get_alignment_table(self) -> pd.DataFrame:
        """Get a numbered alignment table from the numbering and insertions

        Returns
        -------
        pd.DataFrame
            A dataframe with Id, chain_type, scheme and numbering. Values are the amino acid sequences
        """
        all_dataframes = []

        # I'm not sure if there is a more effiecient way to do this other than iterate through the df and pivot each row
        for index in range(len(self)):
            all_dataframes.append(self._pivot_alignment(self.iloc[index]))
        all_dataframes = pd.concat(all_dataframes)
        all_dataframes.columns = list(map(lambda x: str(x[0]) + x[1], all_dataframes.columns.values))
        all_dataframes = all_dataframes.reset_index()

        return self[["Id", "chain_type", "scheme"]].merge(all_dataframes, on="Id").copy()

    def _get_region(self, row, start: int, end: int, segment_name) -> pd.Series:
        with_segment = "".join(
            list(
                map(
                    lambda x: x[-1],
                    list(
                        filter(
                            lambda x: x[0] >= start and x[0] <= end,
                            list(
                                zip(
                                    row["Numbering"],
                                    row["Insertion"],
                                    row["Numbered_Sequence"],
                                )
                            ),
                        )
                    ),
                )
            )
        )
        without_segment = with_segment.replace("-", "")
        return pd.Series(
            {
                f"{segment_name}_gaps": with_segment,
                f"{segment_name}_no_gaps": without_segment,
            }
        )

    def _add_segment_regions(self) -> "NumberingResults":
        """Private method to delineate the framework and cdr boundaries from the numbering

        Returns
        -------
        NumberingResults
            Instance of NumberingResults
        """
        return_frames = []
        for group, sub_df in self.groupby(["scheme", "region_definition", "Chain"]):
            numbering = group[0]
            chain = {"H": "heavy", "KL": "light"}[group[-1]]
            boundaries = group[1]
            numbering_lookup = scheme_numbering[numbering][chain][boundaries]
            for region in [
                "fwr1_aa",
                "cdr1_aa",
                "fwr2_aa",
                "cdr2_aa",
                "fwr3_aa",
                "cdr3_aa",
                "fwr4_aa",
            ]:
                _start = numbering_lookup[f"{region}_start"]
                _end = numbering_lookup[f"{region}_end"]
                sub_df = sub_df.join(self.apply(lambda x: self._get_region(x, _start, _end, region), axis=1))
            return_frames.append(sub_df)
        segmented_df = pd.concat(return_frames).reset_index(drop=True)
        # everything preceding the antibody
        segmented_df["leader"] = segmented_df[["sequence", "seqstart_index"]].apply(lambda x: x[0][: x[1]], axis=1)

        # everything following the antibody. keyword tail will clash with pandas
        segmented_df["follow"] = segmented_df[["sequence", "seqend_index"]].apply(lambda x: x[0][x[1] + 1 :], axis=1)
        return segmented_df

    def _pivot_alignment(self, row: pd.Series) -> pd.DataFrame:
        """Private method to pivot a segmented row into an alignment series

        Parameters
        ----------
        row : pd.Series
            indidual Numbering result row

        Returns
        -------
            pivoted dataframe
        """
        pivoted_df = (
            pd.DataFrame(
                zip(row["Numbering"], row["Insertion"], row["Numbered_Sequence"]),
                columns=["numbering", "insertion", "sequence"],
            )
            .assign(Id=row["Id"])
            .pivot("Id", ["numbering", "insertion"], "sequence")
        )
        return pivoted_df

    # def get_sanatized_antibodies(self):
    #     # drop sequences that don't start at the first amino acid and dont end at the last amino acid.
    #     return self[(self["seqstart_index"] == 0) & (self["seqend_index"] == self["sequence"].str.len() - 1)]

    @staticmethod
    def read_csv(*args, **kwargs):
        return NumberingResults(
            pd.read_csv(
                *args,
                index_col=0,
                dtype=NUMBERING_RESULTS,
                converters={"Numbering": literal_eval, "Insertion": literal_eval, "Numbered_Sequence": literal_eval},
                **kwargs,
            )
        )

    # def drop_bad_numbering(self) -> "NumberingResults":
    #     return self[(self["seqstart_index"] == 0) & (self["seqend_index"] == self["sequence"].str.len() - 1)]
