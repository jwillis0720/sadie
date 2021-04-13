# %%
import pandas as pd
from sadie.airr.airrtable.constants import IGBLAST_AIRR
from sadie.airr.exceptions import MissingAirrColumns
from sadie.airr.airrtable import AirrTable

file_to_test = pd.read_csv("../tests/unit/airr/fixtures/airr_tables/dog_igh.csv.gz")
table = AirrTable(file_to_test)

# %%
