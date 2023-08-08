from sadie.airr import Airr
from sadie.airr.methods import run_mutational_analysis
from sadie.cluster import Cluster

fasta_file = "catnap.fasta"
airr_api = Airr("human", adaptable=True)
airr_table = airr_api.run_fasta(fasta_file)
airr_table_mutational = run_mutational_analysis(airr_table, "kabat", run_multiproc=False)

cluster_api = Cluster(airr_table_mutational, lookup=["cdr1_aa", "cdr2_aa", "cdr3_aa"], pad_somatic=True)

airr_table_with_cluster = cluster_api.cluster(5)
airr_table_with_cluster.sort_values("cluster")
print(airr_table_with_cluster.sort_values("cluster")[["sequence_id", "cluster"]].head(5))
