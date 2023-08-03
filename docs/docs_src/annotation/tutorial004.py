import pandas as pd

from sadie.airr import AirrTable

# use AirrTable method to convert AirrTable.tsv to an AirrTable object
pg9_path = "PG9 AIRR.tsv.gz"


airr_table = AirrTable.read_airr(pg9_path)
print(type(airr_table), isinstance(airr_table, AirrTable))

# or use pandas read_csv method
airr_table_from_pandas = AirrTable(pd.read_csv(pg9_path, sep="\t"))
print(type(airr_table_from_pandas), isinstance(airr_table_from_pandas, AirrTable))
