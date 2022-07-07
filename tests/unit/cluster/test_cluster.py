from sadie.cluster import Cluster
from sadie.airr import AirrTable, LinkedAirrTable
import pandas as pd

from tests.conftest import SadieFixture


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


def test_cluster_with_somatic_pad(fixture_setup:SadieFixture):
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
