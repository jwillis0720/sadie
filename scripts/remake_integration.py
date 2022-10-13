# Use this script if you ever have changes with major integration changes
from sadie.airr import Airr

input_heavy = "tests/data/fixtures/fasta_inputs/catnap_nt_heavy.fasta"
input_light = "tests/data/fixtures/fasta_inputs/catnap_nt_light.fasta"
output_heavy = "tests/data/fixtures/airr_tables/catnap_heavy_airrtable.feather"
output_light = "tests/data/fixtures/airr_tables/catnap_light_airrtable.feather"
airr_api = Airr("human", adaptable=True)
catnap_heavy = airr_api.run_fasta(input_heavy)
catnap_light = airr_api.run_fasta(input_light)
catnap_heavy.to_feather(output_heavy)
catnap_light.to_feather(output_light)
