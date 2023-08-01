# Use Renumbering module
from sadie.renumbering import Renumbering
import pandas as pd

# define a fasta file
catnap_fasta = "catnap_aa_heavy_sample.fasta"


# We wrap these in a function so we can use multiprocessing
def run() -> pd.DataFrame:
    # setup API  object
    renumbering_api = Renumbering(scheme="chothia", region_assign="imgt", run_multiproc=True)

    # run the renumbering on a file
    numbering_table = renumbering_api.run_file(catnap_fasta)

    return numbering_table


if __name__ == "__main__":
    print(run())
