# Std lib
import logging
import os
from pathlib import Path

# Third party
from sadie.reference import Reference

# This module
from sadie.reference.blast import write_blast_db

logger = logging.getLogger("reference")


def make_blast_db_for_internal(df, dboutput):
    """Make a blast database from dataframe"""
    out_fasta = dboutput + ".fasta"
    logger.debug("Writing fasta to {}".format(out_fasta))
    with open(out_fasta, "w") as f:
        for id_, seq in zip(df["gene"], df["sequence"]):
            f.write(">{}\n{}\n".format(id_, seq))
    out_db = out_fasta.split(".fasta")[0]
    write_blast_db(out_fasta, out_db)


def make_internal_annotaion_file(reference: Reference, outpath: Path):
    logger.debug(f"Generating internal annotation file at {outpath}")
    # The internal data file structure goes {db_type}/Ig/internal_path/{species}/

    database = reference.get_dataframe()
    for group, group_df in database.groupby("source"):
        # db_type eg. cutom, imgt

        # get a filtered database for V genes
        filtered_data = group_df.loc[group_df["gene_segment"] == "V"]

        # the species is the actual entity we are using for the annotation, e.g se09 or human
        for species, species_df in filtered_data.groupby("species"):
            # species level database
            species_internal_db_path = os.path.join(outpath, group, "Ig", "internal_data", species)
            logger.debug(f"Found species {species}, using {group} database file")
            if not os.path.exists(species_internal_db_path):
                logger.info(f"Creating {species_internal_db_path}")
                os.makedirs(species_internal_db_path)

            # maybe we requested a chimeric speicies that has more than one sub species, e.g a mouse model
            sub_species_keys = species_df["sub_species"].unique()

            # for sub_species in sub_species_keys:
            #     # only get the common species from our database
            #     sub_filtered = species_df.loc[species_df["sub_species"] == sub_species].copy()

            # if we have hybrid species we shall name them with <species>|gene
            gene_df = species_df.copy()
            if len(sub_species_keys) > 1:
                gene_df["gene"] = gene_df["common"] + "|" + gene_df["gene"]

            index_df = gene_df[
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

            # makes everything an integer. sets gene to index so its not affected
            # add +1 to so we get 1-based indexing
            index_df = (index_df.set_index("gene") + 1).astype("Int64").reset_index()

            # drop anything where there is an na in the annotation idnex
            index_df = index_df.drop(index_df[index_df.isna().any(axis=1)].index)
            scheme = "imgt"
            internal_annotations_file_path = os.path.join(species_internal_db_path, f"{species}.ndm.{scheme}")
            if len(sub_species_keys) > 1:
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
            DB_OUTPATH = os.path.join(species_internal_db_path, f"{species}_{suffix}")
            # Pass the dataframe and write out the blast database
            make_blast_db_for_internal(gene_df, DB_OUTPATH)
