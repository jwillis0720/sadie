from unittest.mock import patch

import numpy as np
import pandas as pd
from Levenshtein import distance as lev_distance
from rapidfuzz.process import cdist as rapidfuzz_cdist

from sadie.airr import AirrTable, LinkedAirrTable
from sadie.cluster import Cluster
from tests.conftest import SadieFixture


def test_threaded_distance_matches_existing_value_normalization():
    cluster = Cluster(AirrTable(pd.DataFrame({"cdr1_aa": ["AAA"]})), lookup=["cdr1_aa"])
    values = (["AAA", "AAB", None, np.nan, pd.NA, 123, "Å"] * 22)[:151]
    df = pd.DataFrame({"cdr1_aa": values}, index=np.arange(10, 312, 2), dtype=object)
    lookup = df[["cdr1_aa"]].to_dict(orient="index")
    rows = list(lookup.values())
    expected = np.array(
        [[lev_distance(str(left["cdr1_aa"]), str(right["cdr1_aa"])) for right in rows] for left in rows],
        dtype=np.float64,
    )

    result = cluster._get_distance_df(df)

    np.testing.assert_array_equal(result, expected)
    assert result.dtype == np.float64
    np.testing.assert_array_equal(result, result.T)
    np.testing.assert_array_equal(np.diag(result), np.zeros(len(df)))


def test_threaded_somatic_distance_matches_legacy_padding():
    values = (["AAA", "AAB", None, np.nan, pd.NA, 123, "Å"] * 22)[:151]
    mutations = (
        [
            np.array(["A1T", "A1T", "C94G"], dtype=object),
            ["A1T", "G93C"],
            [],
            np.array(["G93C", "T120A"], dtype=object),
            ["C94G"],
            ["Å93B"],
            ["X2Y", "T120A"],
        ]
        * 22
    )[:151]
    source = pd.DataFrame({"cdr1_aa": values, "mutations": mutations}, dtype=object)

    for include_only_v_gene in (False, True):
        df = source.copy(deep=True)
        cluster = Cluster(
            AirrTable(df.iloc[:1].copy()),
            lookup=["cdr1_aa"],
            pad_somatic=True,
            include_only_v_gene=include_only_v_gene,
        )
        with patch("sadie.cluster.cluster.cdist", wraps=rapidfuzz_cdist) as cdist_spy:
            result = cluster._get_distance_df(df)
        rows = list(df[["cdr1_aa", "mutations"]].to_dict(orient="index").values())
        expected = np.array(
            [
                [
                    max(
                        lev_distance(str(left["cdr1_aa"]), str(right["cdr1_aa"]))
                        - len(np.intersect1d(left["mutations"], right["mutations"])),
                        0,
                    )
                    for right in rows
                ]
                for left in rows
            ],
            dtype=np.float64,
        )

        np.testing.assert_array_equal(result, expected)
        assert result.dtype == np.float64
        assert {call.kwargs["workers"] for call in cdist_spy.call_args_list} == {-1}


def test_linked_somatic_padding_counts_each_chain(fixture_setup: SadieFixture):
    source = pd.read_feather(fixture_setup.get_catnap_joined_with_mutational_analysis()).iloc[:2]
    linked = LinkedAirrTable(source, key_column="cellid")
    linked["cdr3_aa_heavy"] = ["AAAA", "BBBB"]
    linked.at[0, "mutations_heavy"] = np.array(["A1T", "A1T", "C94G"], dtype=object)
    linked.at[1, "mutations_heavy"] = np.array(["A1T", "C94G"], dtype=object)
    linked.at[0, "mutations_light"] = np.array(["A1T", "D93E"], dtype=object)
    linked.at[1, "mutations_light"] = np.array(["A1T", "D93E"], dtype=object)
    cluster = Cluster(
        linked,
        lookup=["cdr3_aa_heavy"],
        pad_somatic=True,
        include_only_v_gene=True,
    )

    result = cluster._get_distance_df(linked)

    np.testing.assert_array_equal(result, np.array([[0.0, 1.0], [1.0, 0.0]]))
    assert linked.at[0, "mutations_heavy"] == ["A1T", "A1T"]


def test_linked_cluster_preserves_exact_output_with_reconstruction(fixture_setup: SadieFixture):
    ids = ["VRC26.05", "PCT64-18B", "VRC26.06", "PCT64-18C", "VRC26.08", "PCT64-18D"]
    source = pd.read_feather(fixture_setup.get_catnap_joined_with_mutational_analysis())
    linked = LinkedAirrTable(source.set_index("cellid").loc[ids].reset_index(), key_column="cellid")
    expected = pd.DataFrame(linked).loc[[1, 3, 5, 0, 2, 4]].copy()
    expected["cluster"] = [
        "IGHV3-15*01_IGKV3-20*01_1",
        "IGHV3-15*01_IGKV3-20*01_0",
        "IGHV3-15*01_IGKV3-20*01_0",
        "IGHV3-30*03_IGLV1-51*02_2",
        "IGHV3-30*03_IGLV1-51*02_1",
        "IGHV3-30*03_IGLV1-51*02_0",
    ]
    cluster = Cluster(linked, groupby=["v_call_top_heavy", "v_call_top_light"])

    result = cluster.cluster(10)

    pd.testing.assert_frame_equal(pd.DataFrame(result), expected)
    assert result is not cluster.airrtable
    assert isinstance(result, LinkedAirrTable)
    assert result.key_column == "cellid"
    assert result.suffixes == ["_heavy", "_light"]
    assert result.verified
    assert "cluster" not in linked.columns


def test_cluster(heavy_catnap_airrtable, light_catnap_airrtable):
    for table in [heavy_catnap_airrtable, light_catnap_airrtable]:
        cluster = Cluster(table)
        clustered_df = cluster.cluster(10)
        assert "cluster" in clustered_df.columns
        assert isinstance(clustered_df, AirrTable)

    linked = LinkedAirrTable(
        heavy_catnap_airrtable.merge(light_catnap_airrtable, on="cellid", suffixes=["_heavy", "_light"]),
        key_column="cellid",
    )
    cluster = Cluster(
        linked,
        groupby=["v_call_top_heavy", "v_call_top_light"],
        lookup=["cdr1_aa_heavy", "cdr2_aa_heavy", "cdr3_aa_heavy", "cdr1_aa_light", "cdr2_aa_light", "cdr3_aa_light"],
    )
    cluster_df_linked = cluster.cluster(10)
    assert isinstance(cluster_df_linked, LinkedAirrTable)
    assert "cluster" in cluster_df_linked.columns


def test_cluster_with_somatic_pad(fixture_setup: SadieFixture):
    light_airrtable = AirrTable(pd.read_feather(fixture_setup.get_catnap_light_with_mutational_analysis()))
    cluster_api = Cluster(light_airrtable, pad_somatic=True)
    cluster_df_with_pad = cluster_api.cluster(5)
    heavy_airrtable = AirrTable(pd.read_feather(fixture_setup.get_catnap_heavy_with_mutational_analysis()))
    cluster_api = Cluster(heavy_airrtable, pad_somatic=True)
    cluster_df_with_pad = cluster_api.cluster(5)

    joined_airrtable = LinkedAirrTable(
        pd.read_feather(fixture_setup.get_catnap_joined_with_mutational_analysis()), key_column="cellid"
    )
    cluster_api = Cluster(joined_airrtable, pad_somatic=True)
    cluster_df_with_pad = cluster_api.cluster(5)
    assert len(cluster_df_with_pad[cluster_df_with_pad["cellid"].str.startswith("CH0")]["cluster"].unique()) == 1

    # test somatic pad with v_gene_only
    cluster_api = Cluster(joined_airrtable, pad_somatic=True, include_only_v_gene=True)
    cluster_df_with_pad = cluster_api.cluster(5)
    assert len(cluster_df_with_pad[cluster_df_with_pad["cellid"].str.startswith("CH0")]["cluster"].unique()) == 1
