from sadie.airr import Airr

# define a single sequence
pg9_seq = "CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGTCCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGTGGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCCAGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATTTTTGTGTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTATTATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC"

# setup API object
airr_api = Airr("human")

# run sequence and return airr table with sequence_id and sequence
airr_table = airr_api.run_single("PG9", pg9_seq)

# write airr table to a csv
airr_table.to_csv("PG9 AIRR.csv")

# write to a json file
airr_table.to_json("PG9 AIRR.json", orient="records")

# write to a browser friendly html file
airr_table.to_html("PG9 AIRR.html")

# write to an excel file
airr_table.to_excel("PG9 AIRR.xlsx")

# write to a parquet file that is read by spark
airr_table.to_parquet("PG9 AIRR.parquet")

# write to a feather file that has rapid IO
airr_table.to_feather("PG9 AIRR.feather")
