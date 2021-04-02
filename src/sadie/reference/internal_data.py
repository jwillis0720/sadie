# Std lib
import os
import logging
import gzip
import json

# Third party
import pandas as pd

# This module
from .blast import write_blast_db
from .yaml import YamlRef
from ..antibody.exception import BadGene

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


def get_filtered_data(database_json, source, segment, subset):
    if subset == "all":
        function = ["ORF", "F", "P", "I"]
    elif subset == "functional":
        function = ["F"]
    else:
        raise Exception(f"{subset} not a valide option, ['all', 'funtional']")
    return list(
        filter(
            lambda x: x["source"] == source and x["gene_segment"] == segment and x["functional"] in function,
            database_json,
        )
    )


def generate_internal_annotaion_file_from_db(database, outpath):
    logger.debug("Generating from IMGT Internal Database File")
    ig_database = json.load(gzip.open(database, "rt"))
    reference_database = YamlRef()

    # The internal data file structure goes {db_type}/{all|filtered}/Ig/internal_path/{species}/
    # Interate through species and make
    for db_type in reference_database.get_reference_types():
        for subset in reference_database.get_functional_keys(db_type):
            for common in reference_database.get_species_keys(db_type, subset):
                species_internal_db_path = os.path.join(outpath, db_type, subset, "Ig", "internal_data", common)
                logger.debug(f"Found species {common}, using {db_type} database file")
                if not os.path.exists(species_internal_db_path):
                    logger.info(f"Creating {species_internal_db_path}")
                    os.makedirs(species_internal_db_path)

                filtered_data = get_filtered_data(ig_database, db_type, "V", subset)
                sub_species_keys = reference_database.get_sub_species(db_type, subset, common)
                requested_entries = []
                for sub_species in sub_species_keys:
                    gene_segments = reference_database.get_gene_segment(db_type, subset, common, sub_species, "V")
                    request_list = list(
                        filter(lambda x: x["common"] == sub_species and x["gene"] in gene_segments, filtered_data)
                    )
                    if len(request_list) != len(gene_segments):
                        raise BadGene(
                            sub_species, gene_segments, list(map(lambda x: (x["common"], x["gene"]), filtered_data))
                        )
                    requested_entries += request_list

                # if len(sub_species_keys) > 1:
                #     for sub_species in sub_species_keys:

                # else:
                #     filtered_json = get_filtered_data(ig_database, db_type, common, "V", subset)
                #     requested_vs = reference_database.get_gene_segment(db_type, subset, common, common, "V")

                # normalize will flatten nested json
                filt_df = pd.json_normalize(requested_entries)
                if len(filt_df["common"].unique()) > 1:
                    filt_df["gene"] = filt_df["common"] + "|" + filt_df["gene"]

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
                if len(filt_df["common"].unique()) > 1:
                    segment = [i.split("|")[-1].split("-")[0][0:4][::-1][:2] for i in index_df["gene"]]
                else:
                    segment = [i.split("-")[0][0:4][::-1][:2] for i in index_df["gene"]]
                index_df["segment"] = segment
                index_df["weird_buffer"] = 0
                logger.info("Writing to annotation file {}".format(internal_annotations_file_path))
                index_df.to_csv(internal_annotations_file_path, sep="\t", header=False, index=False)
                logger.info("Wrote to annotation file {}".format(internal_annotations_file_path))
                # blast reads these suffixes depending on receptor
                suffix = "V"
                # suffix = "TV_V"
                DB_OUTPATH = os.path.join(species_internal_db_path, f"{common}_{suffix}")
                # Pass the dataframe and write out the blast database
                make_blast_db_for_internal(genes_df, DB_OUTPATH)
