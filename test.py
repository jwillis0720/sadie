from sadie.airr import Airr
import pandas as pd

file = "tests/data/fixtures/airr_tables/dog_igh.tsv.gz"
dog_df = pd.read_csv(file, sep="\t")
airr_api = Airr("dog")
unjoined_df = airr_api.run_dataframe(dog_df, "sequence_id", "sequence")
print(unjoined_df)
