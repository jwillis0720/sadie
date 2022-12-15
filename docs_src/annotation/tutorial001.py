# use Airr module
# import pandas for dataframe handling
import pandas as pd

from sadie.airr import Airr

# define a single sequence
pg9_seq = "CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGTCCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGTGGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCCAGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATTTTTGTGTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTATTATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC"

# setup API  object
airr_api = Airr("human")

# run sequence and return airr table with sequence_id and sequence
airr_table = airr_api.run_single("PG9", pg9_seq)

# output object types
print(type(airr_table))
print(isinstance(airr_table, pd.DataFrame))
