import pandas as pd

from sadie.airr import AirrTable

# write airr table to a csv
airr_table_1 = AirrTable(pd.read_csv("PG9 AIRR.csv"))

# write to a json file
airr_table_2 = AirrTable(pd.read_json("PG9 AIRR.json", orient="records"))

# write to an excel file
airr_table_3 = AirrTable(pd.read_excel("PG9 AIRR.xlsx"))

# write to a parquet file that is read by spark
airr_table_4 = AirrTable(pd.read_parquet("PG9 AIRR.parquet"))

# write to a feather file that has rapid IO
airr_table_5 = AirrTable(pd.read_feather("PG9 AIRR.feather"))
