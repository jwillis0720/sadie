from __future__ import annotations

import logging
import re
from typing import Any, Iterable, List, Optional, Union

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist
from sklearn.cluster import AgglomerativeClustering

from sadie.airr import AirrTable, LinkedAirrTable

logger = logging.getLogger("Cluster")

# ponytail: measured thread crossover; expose tuning only if other hardware disagrees.
_THREAD_MIN_ROWS = 150


class Cluster:
    """Main clustering class.

    This class is used to cluster a given set of data points.
    """

    def __init__(
        self,
        airrtable: Union[AirrTable, LinkedAirrTable],
        linkage: str = "complete",
        groupby: Optional[str] = None,
        lookup: List[str] = ["cdr1_aa", "cdr2_aa", "cdr3_aa"],
        pad_somatic: bool = False,
        include_only_v_gene: bool = False,
    ):
        """Initialize the clustering class.

        Arguments
        ---------
        airrtable (AirrTable, LinkedAirrTable): The airrtable to cluster.
        linkage (str): The linkage method to use. Default is complete. default is complete.
        groupby (str): The linkage method to use. Default is complete. default is complete.
        pad_somatic (bool): Whether to decrease the distance by 1 for every commons sommatic muttaion. Must run mutation analysis firsts

        Raises
        ------
        TypeError
            No airrtable was provided.
        ValueError
            groupby columns must be in the airrtable.
        ValueError
            lookup columns must be in the airrtable
        """
        if not isinstance(airrtable, (AirrTable, LinkedAirrTable)):
            raise TypeError("airrtable table must be a AirrTable or LinkedAirrTable")

        self.include_only_v_gene = include_only_v_gene
        if lookup == ["cdr1_aa", "cdr2_aa", "cdr3_aa"] and isinstance(airrtable, LinkedAirrTable):
            lookup = [i + "_heavy" for i in lookup] + [i + "_light" for i in lookup]

        if groupby is not None:
            diff = set(groupby).difference(set(airrtable.columns))
            if diff:
                raise ValueError(f"groupby column(s) {diff} not found in airrtable")
        if pad_somatic:
            if isinstance(airrtable, LinkedAirrTable):
                if "mutations_heavy" not in airrtable.columns or "mutations_light" not in airrtable.columns:
                    raise ValueError(
                        "pad_somatic requires mutations_heavy and mutations_light in columns. Run mutational analysis first with sadie.arirr.methods"
                    )
                else:
                    self.pad_somatic_values = ["mutations_heavy", "mutations_light"]

            else:
                if "mutations" not in airrtable.columns:
                    raise ValueError(
                        "pad_somatic requires mutations_heavy and mutations_light in columns. Run mutational analysis first with sadie.arirr.methods"
                    )
                else:
                    self.pad_somatic_values = ["mutations"]

        diff = set(lookup).difference(set(airrtable.columns))
        if diff:
            raise ValueError(f"lookup column(s) {diff} not found in airrtable")
        self.airrtable = airrtable
        self.linkage = linkage
        self.groupby = groupby
        self.lookup = lookup
        self.key_column = airrtable.key_column
        self.distance_df = None
        self.model = None
        self.pad_somatic = pad_somatic
        if isinstance(self.airrtable, LinkedAirrTable):
            self._type = "linked"
        else:
            self._type = "unlinked"

    def _get_v_gene_only(self, row: Iterable[str]) -> List[str]:
        row = list(row)
        if row:
            return [i for i in row if int(re.findall(r"\d+", i)[0]) < 94]
        return row

    def _get_distance_df(self, df: pd.DataFrame) -> Any:
        """Given a dataframe, get the N x N pairwise distances using Levenshtein distance of the lookup"""
        if self.pad_somatic:
            _lookup = self.lookup + self.pad_somatic_values
            if self.include_only_v_gene:
                logger.info(f"Including only V genes for {len(df)} rows")
                if self._type == "linked":
                    df["mutations_heavy"] = df["mutations_heavy"].apply(self._get_v_gene_only)
                    df["mutations_light"] = df["mutations_light"].apply(self._get_v_gene_only)
                else:
                    df["mutations"] = df["mutations"].apply(self._get_v_gene_only)
        else:
            _lookup = self.lookup
        df_lookup = df[_lookup].to_dict(orient="index")

        rows = list(df_lookup.values())
        workers = -1 if len(df) >= _THREAD_MIN_ROWS else 1
        distances = np.zeros((len(df), len(df)), dtype=np.float64)
        for metric in self.lookup:
            values = [str(row[metric]) for row in rows]
            distances += cdist(
                values,
                values,
                scorer=Levenshtein.distance,
                workers=workers,
                dtype=np.int32,
            )

        if self.pad_somatic:
            for metric in self.pad_somatic_values:
                rows_by_mutation: dict[str, list[int]] = {}
                for row_index, row in enumerate(rows):
                    mutations: Iterable[str] = row[metric]
                    for mutation in set(mutations):
                        rows_by_mutation.setdefault(mutation, []).append(row_index)
                for row_indexes in rows_by_mutation.values():
                    distances[np.ix_(row_indexes, row_indexes)] -= 1
            np.maximum(distances, 0, out=distances)
        return distances

    def cluster(self, distance_threshold: int = 3) -> Union[AirrTable, LinkedAirrTable]:
        """Cluster the data.

        This method clusters the data using the specified linkage and affinity
        methods.

        Arguments
        ---------
            distance_threshold (int): The maximum distance between two points to be. Default is 3.
        """
        if self.groupby is None:
            self.distance_df = self._get_distance_df(self.airrtable)
            model = AgglomerativeClustering(
                linkage=self.linkage, metric="precomputed", distance_threshold=distance_threshold, n_clusters=None
            )
            model.fit(self.distance_df)
            self.model = model
            # Create the data frame
            self.airrtable["cluster"] = model.labels_
        else:
            cluster_catcher = []
            for g, g_df in self.airrtable.groupby(self.groupby):
                logger.info(f"Clustering group {g}")
                # sub_df = g_df
                sub_df = g_df.copy()
                self.distance_df = self._get_distance_df(sub_df)
                # Calculate the linkage matrix
                model = AgglomerativeClustering(
                    linkage=self.linkage,
                    metric="precomputed",
                    distance_threshold=distance_threshold,
                    n_clusters=None,
                )
                if len(sub_df) == 1:
                    _labels = [0]
                else:
                    model.fit(self.distance_df)
                    _labels = model.labels_
                # Create the data frame
                if isinstance(g, str):
                    labels = list(map(lambda x: f"{g}_{str(x)}", _labels))
                elif isinstance(g, (list, tuple)):
                    _sub_labels = "_".join([str(i) for i in g])
                    labels = list(map(lambda x: f"{_sub_labels}_{str(x)}", _labels))
                else:
                    raise ValueError("groupby must be a string or a list/tuple of strings")
                sub_df["cluster"] = labels
                cluster_catcher.append(sub_df)
            self.airrtable = pd.concat(cluster_catcher)
        if self._type == "unlinked":
            return AirrTable(self.airrtable, key_column=self.key_column)
        return LinkedAirrTable(self.airrtable, key_column=self.key_column)
