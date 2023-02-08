import pandas as pd

from sadie.airr import AirrTable, LinkedAirrTable
from sadie.airr import methods as airr_methods

input_feather = "tests/data/fixtures/airr_tables/bum_igl_assignment_macaque.feather"
output_feather = "tests/data/fixtures/airr_tables/igl_out.feather"
airr_table = AirrTable(pd.read_feather(input_feather))
out_airr_single = airr_methods.run_termini_buffers(airr_table)
igl_df = airr_methods.run_igl_assignment(out_airr_single)
igl_df.to_feather(output_feather)

input_feather = "tests/data/fixtures/airr_tables/bum_link_input.feather"
output_feather = "tests/data/fixtures/airr_tables/bum_link_solution.feather"
lat = LinkedAirrTable(pd.read_feather(fixture_setup.get_bum_link_igl_assignment()))
out_lat = airr_methods.run_termini_buffers(lat)
igl_df = airr_methods.run_igl_assignment(out_lat)
igl_df.to_feather(output_feather)
