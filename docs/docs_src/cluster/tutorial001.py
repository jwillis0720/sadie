import pandas as pd

from sadie.airr.airrtable import LinkedAirrTable
from sadie.cluster import Cluster

# read in the bNAb data
df = pd.read_feather("hiv_bnabs.feather")

# this will create a linked heavy and light airr
# table where the unique name is called the cellid
linked_airr_table = LinkedAirrTable(df, key_column="cellid")


# setup api
cluster_api = Cluster(
    linked_airr_table,
    linkage="single",
    groupby=None,
    lookup=["cdr1_aa_heavy", "cdr2_aa_heavy", "cdr3_aa_heavy", "cdr1_aa_light", "cdr2_aa_light", "cdr3_aa_light"],
    pad_somatic=True,
)

# run the clustering with distance cutoff of 5
clustered_table = cluster_api.cluster(5)
# a column named cluster has been added to the airr table

# print out clusters but first sort by biggest cluster
gb = pd.DataFrame(clustered_table.set_index("cellid")).groupby(["cluster"])
for cluster in gb.size().sort_values(ascending=False).index:
    members = clustered_table.query(f"cluster=={cluster}")["cellid"].to_list()
    print("cluster name:", cluster, "contains", members, end="\n\n")
