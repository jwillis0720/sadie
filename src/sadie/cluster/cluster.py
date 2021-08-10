import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances
from Levenshtein._levenshtein import distance as lev_distance
from sadie.airr import AirrTable, LinkedAirrTable
from typing import Union
import numpy as np


class Cluster:
    """Main clustering class.

    This class is used to cluster a given set of data points.
    """

    def __init__(
        self,
        airrtable: Union[AirrTable, LinkedAirrTable],
        linkage="complete",
        groupby=None,
        lookup=["cdr1_aa", "cdr2_aa", "cdr3_aa"],
    ):
        """Initialize the clustering class.

        Args:
            data (pandas.DataFrame): The data to be clustered.
            n_clusters (int): The number of clusters to be formed.
            linkage (str): The linkage algorithm to be used.
            affinity (str): The affinity algorithm to be used.
            metric (str): The distance metric to be used.
            metric_params (dict): The parameters for the distance metric.
            n_jobs (int): The number of jobs to be used for the clustering.
        """
        if not isinstance(airrtable, (AirrTable, LinkedAirrTable)):
            raise TypeError("airrtable table must be a AirrTable or LinkedAirrTable")
        if groupby is not None:
            diff = set(groupby).difference(set(airrtable.columns))
            if diff:
                raise ValueError(f"groupby column(s) {diff} not found in airrtable")
        diff = set(lookup).difference(set(airrtable.columns))
        if diff:
            raise ValueError(f"lookup column(s) {diff} not found in airrtable")
        self.airrtable = airrtable
        self.linkage = linkage
        self.groupby = groupby
        self.lookup = lookup
        self.key_column = airrtable.key_column
        if isinstance(self.airrtable, LinkedAirrTable):
            self._type = "linked"
        else:
            self._type = "unlinked"

    def _get_distance_df(self, df):
        """Given a dataframe, get the N x N pairwise distances using Levenshtein distance of the lookup"""
        df_lookup = df[self.lookup].to_dict(orient="index")

        def calc_lev(x, y):
            dist = 0
            for metric in self.lookup:
                dist += lev_distance(str(df_lookup[x[0]][metric]), str(df_lookup[y[0]][metric]))
            return dist

        X = np.array(df.index).reshape(-1, 1)
        return pairwise_distances(X, metric=calc_lev, n_jobs=-1)

    def cluster(self, distance_threshold=3):
        """Cluster the data.

        This method clusters the data using the specified linkage and affinity
        methods.
        """
        if self.groupby is None:
            distance_df = self._get_distance_df(self.airrtable)
            model = AgglomerativeClustering(
                linkage=self.linkage, affinity="precomputed", distance_threshold=distance_threshold, n_clusters=None
            )
            model.fit(distance_df)

            # Create the data frame
            self.airrtable["cluster"] = model.labels_
        else:
            cluster_catcher = []
            for g, g_df in self.airrtable.groupby(self.groupby):
                distance_df = self._get_distance_df(g_df)
                # Calculate the linkage matrix
                model = AgglomerativeClustering(
                    linkage=self.linkage,
                    affinity="precomputed",
                    distance_threshold=distance_threshold,
                    n_clusters=None,
                )
                if len(g_df) == 1:
                    _labels = [0]
                else:
                    model.fit(distance_df)
                    _labels = model.labels_
                # Create the data frame
                if isinstance(g, str):
                    labels = list(map(lambda x: f"{g}_{str(x)}", _labels))
                elif isinstance(g, (list, tuple)):
                    _sub_labels = "_".join([str(i) for i in g])
                    labels = list(map(lambda x: f"{_sub_labels}_{str(x)}", _labels))
                else:
                    raise ValueError("groupby must be a string or a list/tuple of strings")
                g_df["cluster"] = labels
                cluster_catcher.append(g_df)
            self.airrtable = pd.concat(cluster_catcher)
        if self._type == "unlinked":
            return AirrTable(self.airrtable, key_column=self.key_column)
        return LinkedAirrTable(self.airrtable, key_column=self.key_column)
