from sadie.cluster.cluster import Cluster
from sadie.airr import AirrTable, LinkedAirrTable
import pandas as pd


def test_cluster():
    heavy = AirrTable(pd.read_csv("tests/integration/airr/fixtures/catnap_heavy_airrtable.csv.gz", index_col=0))
    light = AirrTable(pd.read_csv("tests/integration/airr/fixtures/catnap_light_airrtable.csv.gz", index_col=0))
    cluster = Cluster(heavy)
    clustered_df = cluster.cluster(10)
    assert "cluster" in clustered_df.columns
    cluster = Cluster(heavy, groupby=["v_call_top"])
    clustered_df = cluster.cluster(10)
    cluster = Cluster(light, groupby=["v_call_top"])
    clustered_df = cluster.cluster(10)
    assert "cluster" in clustered_df.columns

    linked = LinkedAirrTable(
        pd.read_csv("tests/unit/airr/fixtures/airr_tables/linked_airr_table_dummy.csv.gz", index_col=0)
    )
    cluster = Cluster(
        linked,
        groupby=["v_call_top_heavy", "v_call_top_light"],
        lookup=["cdr1_aa_heavy", "cdr2_aa_heavy", "cdr3_aa_heavy", "cdr1_aa_light", "cdr2_aa_light", "cdr3_aa_light"],
    )
    cluster_df_linked = cluster.cluster(10)
    assert isinstance(cluster_df_linked, LinkedAirrTable)
    assert "cluster" in cluster_df_linked.columns
