import pandas as pd

from sadie.airr import AirrTable

# use AirrTable method to convert AirrTable.tsv to an AirrTable object
airr_table = AirrTable.read_airr("PG9 AIRR.tsv.gz")
print(type(airr_table), isinstance(airr_table, AirrTable))

# or use pandas read_csv method
airr_table_from_pandas = AirrTable(pd.read_csv("PG9 AIRR.tsv.gz", sep="\t"))
print(type(airr_table_from_pandas), isinstance(airr_table_from_pandas, AirrTable))
print(airr_table_from_pandas == airr_table)
