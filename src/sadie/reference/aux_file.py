"""Auxillary level function interfaces"""
from .yaml import YamlRef

import logging
import os
import gzip
import json
import pandas as pd

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


def make_auxillary_file(database, outdir):
    """
    Given the imgt file structure object and aux path, make the aux data
    #"""
    ig_database = json.load(gzip.open(database, "rt"))
    reference_database = YamlRef()
    for database in set(reference_database.yaml.keys()):
        receptor_aux_dir = os.path.join(outdir, f"{database}/aux_db")
        if not os.path.exists(receptor_aux_dir):
            logger.info("Creating %s", receptor_aux_dir)
            os.makedirs(receptor_aux_dir)
        # for functional in set(reference_database.yaml[database].keys()):
        commons = set(reference_database.yaml[database]["all"].keys())
        for common in commons:
            dataset = []
            sub_species_keys = reference_database.get_sub_species(database, "all", common)
            if len(sub_species_keys) > 1:
                chimera = True
            else:
                chimera = False
            for sub_species in sub_species_keys:
                ref_lookup = reference_database.get_gene_segment(database, "all", common, sub_species, "J")
                pre_filetered = get_filtered_data(ig_database, sub_species, database, "J")
                dataset += list(filter(lambda x: x["gene"] in ref_lookup, pre_filetered))
            aux_file_species = os.path.join(receptor_aux_dir, f"{common}_gl.aux")
            common_df = pd.json_normalize(dataset)
            common_df = common_df[common_df["imgt.end_cdr3_nt"] != ""]
            common_df.loc[:, "reading_frame"] = common_df["imgt.reading_frame"].astype(int) - 1
            common_df.loc[:, "left_over"] = common_df[["imgt.j_gene_nt", "reading_frame"]].apply(
                _determine_left_over, axis=1
            )
            common_df["marker"] = common_df["gene"].str.split("-").str.get(0).str[0:4].str[::-1].str[:2]
            common_df["end"] = common_df["imgt.end_cdr3_nt"].astype(int) - 1
            if chimera:
                common_df["gene"] = common_df["common"] + "|" + common_df["gene"]
            common_df[["gene", "reading_frame", "marker", "end", "left_over"]].to_csv(
                aux_file_species, sep="\t", header=None, index=False
            )
            logger.info(f"Wrote aux to {aux_file_species}")
    # ig_database = json.load(gzip.open(database, "rt"))
    # for common, source in itertools.product(
    #     get_species_from_database(ig_database),
    #     get_databases_types(ig_database),
    # ):
    #     filtered_data = get_filtered_data(ig_database, common, source, "J")
    #     if not filtered_data:
    #         logger.debug(f"No info for {common}-{source} j chain")
    #         continue
    #     logger.info("Writing J Auxiliary data to dataframe %s", aux_file_species)
