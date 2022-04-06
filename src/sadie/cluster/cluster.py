import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances
from Levenshtein._levenshtein import distance as lev_distance
from sadie.airr import AirrTable, LinkedAirrTable
from typing import Any, List, Optional, Union
import numpy as np
import logging
import numpy.typing as npt


logger = logging.getLogger("Cluster")


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

    def _get_distance_df(self, df: pd.DataFrame) -> Any:
        """Given a dataframe, get the N x N pairwise distances using Levenshtein distance of the lookup"""
        if self.pad_somatic:
            _lookup = self.lookup + self.pad_somatic_values
        else:
            _lookup = self.lookup
        df_lookup = df[_lookup].to_dict(orient="index")

        def calc_lev(x: npt.ArrayLike, y: npt.ArrayLike) -> float:
            dist = 0
            for metric in self.lookup:
                dist += lev_distance(str(df_lookup[x[0]][metric]), str(df_lookup[y[0]][metric]))  # type: ignore[index]
            if self.pad_somatic and x[0] != y[0]:  # type: ignore[index]
                if len(self.pad_somatic_values) == 2:
                    _mutations_1_heavy = df_lookup[x[0]][self.pad_somatic_values[0]]  # type: ignore[index]
                    _mutations_2_heavy = df_lookup[y[0]][self.pad_somatic_values[0]]  # type: ignore[index]
                    _mutations_1_light = df_lookup[x[0]][self.pad_somatic_values[1]]  # type: ignore[index]
                    _mutations_2_light = df_lookup[y[0]][self.pad_somatic_values[1]]  # type: ignore[index]
                    subtract_heavy = len(np.intersect1d(_mutations_1_heavy, _mutations_2_heavy))
                    subtract_light = len(np.intersect1d(_mutations_1_light, _mutations_2_light))
                    subtract_all = subtract_heavy + subtract_light
                else:
                    _mutations_1 = df_lookup[x[0]][self.pad_somatic_values[0]]  # type: ignore[index]
                    _mutations_2 = df_lookup[y[0]][self.pad_somatic_values[0]]  # type: ignore[index]
                    subtract_all = len(np.intersect1d(_mutations_1, _mutations_2))
                dist -= subtract_all
            return max(dist, 0)

        X: npt.ArrayLike = np.array(df.index).reshape(-1, 1)
        return pairwise_distances(X, metric=calc_lev, n_jobs=-1)

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
                linkage=self.linkage, affinity="precomputed", distance_threshold=distance_threshold, n_clusters=None
            )
            model.fit(self.distance_df)
            self.model = model
            # Create the data frame
            self.airrtable["cluster"] = model.labels_
        else:
            cluster_catcher = []
            for g, g_df in self.airrtable.groupby(self.groupby):
                # sub_df = g_df
                sub_df = g_df.copy()
                self.distance_df = self._get_distance_df(sub_df)
                # Calculate the linkage matrix
                model = AgglomerativeClustering(
                    linkage=self.linkage,
                    affinity="precomputed",
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
