# Use Renumbering module
import pandas as pd

from sadie.renumbering import Renumbering

# define a single sequence
vrc26_seq = "QKQLVESGGGVVQPGRSLTLSCAASQFPFSHYGMHWVRQAPGKGLEWVASITNDGTKKYHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFCVRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDIWGHGTMVIVS"


# We wrap these in a function so we can use multiprocessing
def run() -> pd.DataFrame:
    # setup API  object
    renumbering_api = Renumbering(scheme="chothia", region_assign="imgt", run_multiproc=True)

    # run sequence and return airr table with sequence_id and sequence
    numbering_table = renumbering_api.run_single("VRC26.27", vrc26_seq)

    # get the handy dandy alignment table
    return numbering_table.get_alignment_table()


if __name__ == "__main__":
    print(run())
