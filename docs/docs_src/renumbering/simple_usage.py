# Use Renumbering module
# import pandas for dataframe handling
import pandas as pd

from sadie.renumbering import Renumbering

# define a single sequence
vrc26_seq = "QKQLVESGGGVVQPGRSLTLSCAASQFPFSHYGMHWVRQAPGKGLEWVASITNDGTKKYHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFCVRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDIWGHGTMVIVS"

# setup API  object
renumbering_api = Renumbering(scheme="chothia", region_assign="imgt", run_multiproc=True)

# run sequence and return airr table with sequence_id and sequence
numbering_table = renumbering_api.run_single("VRC26.27", vrc26_seq)

# output object types
print(numbering_table)
