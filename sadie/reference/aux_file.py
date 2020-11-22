"""Auxillary level function interfaces"""

import logging
import re
import os

import pandas as pd
from Bio import SeqIO

# module
from .settings import MOTIF_LOOKUP

logger = logging.getLogger(__file__)


def make_auxillary_file(engine, aux_path):
    """
    Given the imgt file structure object and aux path, make the aux data
    """

    js_only = pd.read_sql("j_segment_imgt", con=engine, index_col="index")
    for species, species_df in js_only.groupby("common"):
        aux_file_species = os.path.join(aux_path, "{}_gl.aux".format(species))
        logger.info("Writing J Auxiliary data to dataframe %s", aux_file_species)
        species_df["marker"] = species_df["gene"].str.split("-").str.get(0).str[0:4].str[::-1].str[:2]
        species_df.loc[:, "reading_frame"] = species_df["reading_frame"].astype(int) - 1
        species_df[["gene", "reading_frame", "marker", "end_cdr3_nt_index"]].to_csv(
            aux_file_species, sep="\t", header=None, index=False
        )
