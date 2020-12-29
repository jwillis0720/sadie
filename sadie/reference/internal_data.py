# Std lib
import os
import logging
import gzip
import json

# Third party
import pandas as pd

# This module
from .blast import write_blast_db

logger = logging.getLogger(__name__)


def make_blast_db_for_internal(df, dboutput):
    """Make a blast database from dataframe"""
    out_fasta = dboutput + ".fasta"
    logger.debug("Writing fasta to {}".format(out_fasta))
    with open(out_fasta, "w") as f:
        for id_, seq in zip(df["gene"], df["sequence"]):
            f.write(">{}\n{}\n".format(id_, seq))
    out_db = out_fasta.split(".fasta")[0]
    write_blast_db(out_fasta, out_db)


def get_databases_types(database_json):
    return list(set(map(lambda x: x["source"], database_json)))


def get_species_from_database(database_json):
    return list(set(map(lambda x: x["common"], database_json)))


def get_filtered_data(database_json, source, common, segment):
    return list(
        filter(
            lambda x: x["source"] == source and x["common"] == common and x["gene_segment"] == segment,
            database_json,
        )
    )


def generate_internal_annotaion_file_from_db(database, outpath, only_functional):
    logger.debug("Generating from IMGT Internal Database File")
    ig_database = json.load(gzip.open(database, "rt"))

    available_species = get_species_from_database(ig_database)
    logger.debug("have the following species %s", available_species)

    # The internal data file structure goes internal_path/{species}/
    # Interate through species and make
    for db_type in get_databases_types(ig_database):
        for common in available_species:
            species_internal_db_path = os.path.join(outpath, db_type, "Ig", "internal_data", common)
            logger.debug("Found species %s, using imgt database file", common)
            if not os.path.exists(species_internal_db_path):
                logger.info("Creating {}".format(species_internal_db_path))
                os.makedirs(species_internal_db_path)

            filtered_json = get_filtered_data(ig_database, db_type, common, "V")
            if not filtered_json:
                logger.info(f"No entries for {db_type}-{common}-V")
                continue

            # normalize will flatten nested json
            filt_df = pd.json_normalize(filtered_json)

            index_df = filt_df[
                [
                    "gene",
                    "imgt.fwr1_nt_index_start",
                    "imgt.fwr1_nt_index_end",
                    "imgt.cdr1_nt_index_start",
                    "imgt.cdr1_nt_index_end",
                    "imgt.fwr2_nt_index_start",
                    "imgt.fwr2_nt_index_end",
                    "imgt.cdr2_nt_index_start",
                    "imgt.cdr2_nt_index_end",
                    "imgt.fwr3_nt_index_start",
                    "imgt.fwr3_nt_index_end",
                ]
            ]
            genes_df = filt_df[["gene", "imgt.v_gene_nt"]].rename({"imgt.v_gene_nt": "sequence"}, axis=1)
            scheme = "imgt"
            internal_annotations_file_path = os.path.join(species_internal_db_path, f"{common}.ndm.{scheme}")
            segment = [i.split("-")[0][0:4][::-1][:2] for i in index_df["gene"]]
            index_df["segment"] = segment
            index_df["weird_buffer"] = 0
            logger.info("Writing to annothation file {}".format(internal_annotations_file_path))
            index_df.to_csv(internal_annotations_file_path, sep="\t", header=False, index=False)
            logger.info("Wrote to annothation file {}".format(internal_annotations_file_path))
            # blast reads these suffixes depending on receptor
            suffix = "V"
            # suffix = "TV_V"
            DB_OUTPATH = os.path.join(species_internal_db_path, f"{common}_{suffix}")
            # Pass the dataframe and write out the blast database
            make_blast_db_for_internal(genes_df, DB_OUTPATH)
