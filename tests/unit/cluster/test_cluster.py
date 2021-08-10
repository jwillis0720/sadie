from sadie.cluster import Cluster
from sadie.airr import AirrTable, LinkedAirrTable


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
