"""Auxillary level function interfaces"""

import logging
import os
import pandas as pd
from ..antibody.exception import BadGene
from ..reference import get_loaded_database
from yaml import load as yml_load

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


# module

logger = logging.getLogger(__file__)


def get_databases_types(database_json):
    return list(set(map(lambda x: x["source"], database_json)))


def get_species_from_database(database_json):
    return list(set(map(lambda x: x["common"], database_json)))


def get_filtered_data(database_json, common, source, segment):
    return list(
        filter(
            lambda x: x["source"] == source and x["common"] == common and x["gene_segment"] == segment,
            database_json,
        )
    )


def _determine_left_over(X):
    nt = X[0]
    reading_frame = X[1]
    return len(nt[reading_frame:]) % 3


def make_auxillary_file(reference, outpath):
    """
    Given the imgt file structure object and aux path, make the aux data
    #"""
    logger.debug("Generating from IMGT Internal Database File")
    reference = yml_load(open(reference), Loader=Loader)
    database = get_loaded_database()
    for db_type in reference:
        receptor_aux_dir = os.path.join(outpath, f"{db_type}/aux_db")

        if not os.path.exists(receptor_aux_dir):
            logger.info(f"Creating {receptor_aux_dir}")
            os.makedirs(receptor_aux_dir)
        # for functional in set(reference_database.yaml[database].keys()):
        for common in reference[db_type]:
            dataset = []
            sub_species_keys = reference[db_type][common].keys()
            if len(sub_species_keys) > 1:
                chimera = True
            else:
                chimera = False
            for sub_species in sub_species_keys:
                requested_j = list(filter(lambda x: x[3] == "J", reference[db_type][common][sub_species]))
                have_j = list(filter(lambda x: x["common"] == sub_species, database[db_type]))
                new = list(filter(lambda x: x["gene"] in requested_j, have_j))
                if len(new) != len(requested_j):
                    raise BadGene(sub_species, requested_j, have_j)
                dataset += new

            aux_file_species = os.path.join(receptor_aux_dir, f"{common}_gl.aux")
            common_df = pd.json_normalize(dataset).fillna("")
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
