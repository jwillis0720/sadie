# Std lib
import os
import logging

# Third party
import pandas as pd

# This module
from .settings import IMGT_LOOKUP
from .settings import IMGT_DEF_nt
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


def generate_internal_annotaion_file_from_db(engine, INTERNAL_DATA_PATH, only_functional):
    species = IMGT_LOOKUP.keys()
    logger.debug("Generating from IMGT Internal Database File")

    imgt_db_df = pd.read_sql("v_segment_imgt", con=engine, index_col="index")
    available_species = list(imgt_db_df["common"].unique())
    logger.debug("have the following species %s", available_species)

    # The internal data file structure goes internal_path/{species}/
    # Interate through species and make
    for species, species_df in imgt_db_df.groupby("common"):
        species_internal_db_path = os.path.join(INTERNAL_DATA_PATH, species)

        logger.debug("Found species %s, using imgt database file", species)
        if not os.path.exists(species_internal_db_path):
            logger.info("Creating {}".format(species_internal_db_path))
            os.makedirs(species_internal_db_path)

        # this will be used in making our blast database, since we want everything that exists in blast to have an internal reference
        species_for_blast = species_df[species_df["region_definiton"] == "imgt"][["gene", "v_gene_nt"]]
        species_for_blast.loc[:, "receptor"] = species_for_blast["gene"].str[0:2]
        species_for_blast.rename({"v_gene_nt": "sequence"}, axis=1, inplace=1)
        species_for_blast_gb = species_for_blast.groupby(["gene"])
        more_than_one_gene = species_for_blast_gb.size()[species_for_blast_gb.size() > 1].index
        if not more_than_one_gene.empty:
            logger.warning(
                f"Warning: {species}-{list(more_than_one_gene)} contains multiple entries. The most likely cause is more than one latin name sharing this common name"
            )

        scheme = "imgt"
        internal_annotations_file_path = os.path.join(
            species_internal_db_path, f"{species}.ndm.{scheme}".format(species)
        )
        internal_df = species_df.groupby("gene").head(1)[
            [
                "gene",
                "fwr1_nt_index_start",
                "fwr1_nt_index_end",
                "cdr1_nt_index_start",
                "cdr1_nt_index_end",
                "fwr2_nt_index_start",
                "fwr2_nt_index_end",
                "cdr2_nt_index_start",
                "cdr2_nt_index_end",
                "fwr3_nt_index_start",
                "fwr3_nt_index_end",
            ]
        ]
        logger.info("Writing to annothation file {}".format(internal_annotations_file_path))
        internal_df.loc[:, "segment"] = internal_df["gene"].str.split("-").str.get(0).str[0:4].str[::-1].str[:2]
        internal_df.loc[:, "weird_buffer"] = 0
        internal_df.to_csv(internal_annotations_file_path, sep="\t", header=False, index=False)
        logger.info("Wrote to annothation file {}".format(internal_annotations_file_path))
        logger.info("Making internal files for other schemes")
        for scheme in ["kabat", "abm", "contact", "chothia", "scdr"]:
            scheme_df = pd.read_sql(f"v_segment_{scheme}", con=engine, index_col="index")
            scheme_df_species = scheme_df[scheme_df["common"] == species]
            if scheme_df_species.empty:
                logger.warning(f"{scheme} for {species} V segment annotations is empty...very sad.")
                continue

            #     # # anotations file path
            internal_annotations_file_path = os.path.join(
                species_internal_db_path, f"{species}.ndm.{scheme}".format(species)
            )
            scheme_df_species.loc[:, "segment"] = (
                scheme_df_species["gene"].str.split("-").str.get(0).str[0:4].str[::-1].str[:2]
            )
            scheme_df_species.loc[:, "weird_buffer"] = 0
            logger.debug("Writing to annothation file {}".format(internal_annotations_file_path))
            if len(scheme_df_species["gene"].isin(internal_df["gene"])) != len(scheme_df_species):
                raise Exception(
                    f"There is some error in the {scheme} for {species}. All of the {scheme} is not found in the imgt_df."
                )

            # print(internal_df, scheme_df_species)
            scheme_df_species[internal_df.columns].to_csv(
                internal_annotations_file_path, sep="\t", header=False, index=False
            )
            logger.info("Wrote to annothation file {}".format(internal_annotations_file_path))
        # Unfortunately we have to regropup by receptor when we make the fasta file
        for receptor, receptor_df in species_for_blast.groupby("gene").head(1).groupby("receptor"):
            # blast reads these suffixes depending on receptor
            if receptor == "IG":
                suffix = "V"
            else:
                suffix = "TR_V"
            DB_OUTPATH = os.path.join(species_internal_db_path, species + "_{}".format(suffix))
            # Pass the dataframe and write out the blast database
            make_blast_db_for_internal(receptor_df, DB_OUTPATH)
