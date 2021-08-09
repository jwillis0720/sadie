"""Auxillary level function interfaces"""

import logging
import os
from pathlib import Path

from sadie.reference import Reference

logger = logging.getLogger(__file__)


def _determine_left_over(X):
    nt = X[0]
    reading_frame = X[1]
    return len(nt[reading_frame:]) % 3


def make_auxillary_file(reference: Reference, outpath: Path):
    """
    Given the imgt file structure object and aux path, make the aux data
    """
    database = reference.get_dataframe()
    for group, group_df in database.groupby("source"):
        receptor_aux_dir = os.path.join(outpath, f"{group}/aux_db")

        if not os.path.exists(receptor_aux_dir):
            logger.info(f"Creating {receptor_aux_dir}")
            os.makedirs(receptor_aux_dir)
        for species, species_df in group_df.groupby("species"):
            sub_species_keys = species_df["sub_species"].unique()
            if len(sub_species_keys) > 1:
                chimera = True
            else:
                chimera = False

            aux_file_species = os.path.join(receptor_aux_dir, f"{species}_gl.aux")
            common_df = species_df[species_df["gene_segment"] == "J"].copy()
            bad_remainders = common_df[(common_df["imgt.remainder"].isna())]
            if not bad_remainders.empty:
                logger.warning(f"Had to drop {bad_remainders.shape[0]} rows due to bad remainder for {group}-{species}")
                common_df.drop(bad_remainders.index, inplace=True)
            common_df = common_df[(common_df["imgt.cdr3_end"] != "")]
            common_df.loc[:, "reading_frame"] = common_df["imgt.reading_frame"].astype(int)
            common_df.loc[:, "left_over"] = common_df["imgt.remainder"].astype(int)
            common_df.loc[:, "end"] = common_df["imgt.cdr3_end"].astype(int) - 1
            common_df["marker"] = common_df["gene"].str.split("-").str.get(0).str[0:4].str[::-1].str[:2]
            if chimera:
                common_df["gene"] = common_df["common"] + "|" + common_df["gene"]
            common_df[["gene", "reading_frame", "marker", "end", "left_over"]].to_csv(
                aux_file_species, sep="\t", header=None, index=False
            )
            logger.info(f"Wrote aux to {aux_file_species}")
