# Use Renumbering module
from sadie.renumbering import Renumbering
from multiprocessing import freeze_support

# define a single sequence
vrc26_seq = "QKQLVESGGGVVQPGRSLTLSCAASQFPFSHYGMHWVRQAPGKGLEWVASITNDGTKKYHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFCVRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDIWGHGTMVIVS"

# setup API  object
renumbering_api = Renumbering(scheme="chothia", region_assign="imgt", run_multiproc=False)

# run sequence and return airr table with sequence_id and sequence
numbering_table = renumbering_api.run_single("VRC26.27", vrc26_seq)

# get the handy dandy alignment table
numbering_table.get_alignment_table()
