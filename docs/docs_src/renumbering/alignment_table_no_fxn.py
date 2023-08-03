# Use Renumbering module
import pandas as pd

from sadie.renumbering import Renumbering

# define a single sequence
vrc26_seq = "QKQLVESGGGVVQPGRSLTLSCAASQFPFSHYGMHWVRQAPGKGLEWVASITNDGTKKYHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFCVRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDIWGHGTMVIVS"

# setup API  object
renumbering_api = Renumbering(scheme="chothia", region_assign="imgt", run_multiproc=True)
# run sequence and return airr table with sequence_id and sequence
numbering_table = renumbering_api.run_single("VRC26.27", vrc26_seq)

# get the handy dandy alignment table
numbering_table.get_alignment_table()
