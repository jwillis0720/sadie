# Std lib
import os
import logging

# Third party
import pandas as pd

# This module
from .blast import write_blast_db
from ..antibody.exception import BadGene
from ..reference import get_loaded_database
from yaml import load as yml_load

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader
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


def get_filtered_data(data, segment):
    return list(
        filter(
            lambda x: x["gene_segment"] == segment,
            data,
        )
    )


def generate_internal_annotaion_file_from_db(reference, outpath):
    logger.debug("Generating from IMGT Internal Database File")
    reference = yml_load(open(reference), Loader=Loader)
    database = get_loaded_database()

    # The internal data file structure goes {db_type}/Ig/internal_path/{species}/
    # Interate through species and make
    for db_type in reference.keys():
        # db_type eg. cutom, imgt

        # get a filtered database for V genes
        filtered_data = get_filtered_data(database[db_type], "V")

        for common in reference[db_type]:
            # species level database
            common_data = reference[db_type][common]

            species_internal_db_path = os.path.join(outpath, db_type, "Ig", "internal_data", common)
            logger.debug(f"Found species {common}, using {db_type} database file")
            if not os.path.exists(species_internal_db_path):
                logger.info(f"Creating {species_internal_db_path}")
                os.makedirs(species_internal_db_path)

            # maybe we requested a chimeric speicies that has more than one sub species
            sub_species_keys = common_data.keys()

            # here are the requested entries we are asking to put in our blast database
            requested_entries = []
            for sub_species in sub_species_keys:
                # only get the common species from our database
                sub_filtered = list(filter(lambda x: x["common"] == sub_species, filtered_data))
                gene_segments = list(filter(lambda x: x[3] == "V", common_data[sub_species]))
                request_list = list(filter(lambda x: x["gene"] in gene_segments, sub_filtered))
                if len(request_list) != len(gene_segments):
                    accepted_genes = list(map(lambda x: x["gene"], sub_filtered))
                    raise BadGene(sub_species, gene_segments, accepted_genes)
                requested_entries += request_list

                # normalize will flatten nested json
                filt_df = pd.json_normalize(requested_entries)

                if filt_df.empty:
                    logger.warning(f"{common}:{sub_species} has no V genes in {db_type} database")
                    continue
                # if we have hybrid species we shall name them with <species>|gene
                if len(filt_df["common"].unique()) > 1:
                    filt_df["gene"] = filt_df["common"] + "|" + filt_df["gene"]

                index_df = filt_df[
                    [
                        "gene",
                        "imgt.fwr1_start",
                        "imgt.fwr1_end",
                        "imgt.cdr1_start",
                        "imgt.cdr1_end",
                        "imgt.fwr2_start",
                        "imgt.fwr2_end",
                        "imgt.cdr2_start",
                        "imgt.cdr2_end",
                        "imgt.fwr3_start",
                        "imgt.fwr3_end",
                    ]
                ].copy()
                index_df = (index_df.set_index("gene") + 1).astype("Int64").reset_index()
                index_df = index_df.drop(index_df[index_df.isna().any(axis=1)].index)
                genes_df = filt_df.copy()
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
