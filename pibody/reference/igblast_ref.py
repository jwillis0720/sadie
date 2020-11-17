import logging
import os

# third party
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# module
from .settings import BLAST_CONVENTION
from .blast import write_blast_db

logger = logging.getLogger(__file__)


def write_out_fasta(sequences, directory, species, gene_seg):
    logger = logging.getLogger(__file__)
    output_fasta = os.path.join(directory, "{}_{}.fasta".format(species, gene_seg))
    logger.debug("output fasta {}".format(output_fasta))

    # Im sure this will come back to haunt me, but if, we've seen the name, save that
    seen_short_names = []
    with open(output_fasta, "w") as f:
        for sequence in sequences:
            name = sequence.name
            if "|" not in name:
                # we may have preparsed this if its coming from database
                short_name = name
            else:
                short_name = name.split("|")[1]
            if short_name in seen_short_names:
                logger.warning(
                    "DANGER!!! DANGER!! We have already seen this gene name once, probably a duplicate species {}-{}-{}\n{}!!".format(
                        species, gene_seg, short_name, sequence.description
                    )
                )
                continue
            seq = str(sequence.seq).replace(".", "")
            f.write(">{}\n{}\n".format(short_name, seq))
            seen_short_names.append(short_name)

    return output_fasta


def make_igblast_ref_database(engine, blast_dir, only_functional):
    """[summary]

    Parameters
    ----------
    engine : [type]
        [description]
    blast_dir : [type]
        [description]
    only_functional : [type]
        [description]
    """
    # The blast DB groups by V,D and J
    imgt_db_unique = pd.read_sql("imgt_unique", con=engine, index_col="index")

    # clean
    v_segments_imgt = pd.read_sql("v_segment_imgt", con=engine, index_col="index")
    j_segments_imgt = pd.read_sql("j_segment_imgt", con=engine, index_col="index")
    imgt_db_unique["receptor"] = imgt_db_unique["gene"].str[0:2]
    imgt_db_unique["segment"] = imgt_db_unique["gene"].str.split("-").str.get(0).str[0:4].str[-1]

    for receptor in ["IG", "TR"]:
        recptor_translate_name = BLAST_CONVENTION[receptor]
        receptor_blast_dir = os.path.join(blast_dir, recptor_translate_name)
        if not os.path.exists(receptor_blast_dir):
            logger.info("Creating %s", receptor_blast_dir)
            os.makedirs(receptor_blast_dir)

        for species in imgt_db_unique["common"].unique():
            species_blast_dir = os.path.join(receptor_blast_dir, species)
            if not os.path.exists(species_blast_dir):
                logger.info("Creating %s", species_blast_dir)
                os.makedirs(species_blast_dir)

            ## V segment
            species_segment_df = v_segments_imgt[v_segments_imgt["common"] == species]
            species_segment_df = species_segment_df.groupby(["gene"]).head(1)
            joined_seqs = []
            for v, seq in zip(species_segment_df["gene"], species_segment_df["v_gene_nt"]):
                joined_seqs.append(SeqRecord(Seq(seq), name=v, id=v))
            fasta_file = write_out_fasta(joined_seqs, species_blast_dir, species, "V")
            write_blast_db(fasta_file, fasta_file.split(".fasta")[0])

            # D segment - Pull from imgt db unique
            species_segment_df = imgt_db_unique[
                (imgt_db_unique["common"] == species) & (imgt_db_unique["segment"] == "D")
            ]
            if species_segment_df.empty:
                logger.warning(f"{receptor}-{species}-D is empty...sad")
            else:
                species_segment_df = species_segment_df.groupby(["gene"]).head(1)
                joined_seqs = []
                for d, seq in zip(species_segment_df["gene"], species_segment_df["nt_sequence_no_gaps"]):
                    joined_seqs.append(SeqRecord(Seq(seq), name=d, id=d))
                fasta_file = write_out_fasta(joined_seqs, species_blast_dir, species, "D")
                write_blast_db(fasta_file, fasta_file.split(".fasta")[0])

            # J segment - Pull from imgt db unique
            species_segment_df = j_segments_imgt[(j_segments_imgt["common"] == species)]
            if species_segment_df.empty:
                logger.warning(f"{receptor}{species}-J is empty...sad")
            else:
                species_segment_df = species_segment_df.groupby(["gene"]).head(1)
                joined_seqs = []
                for j, seq in zip(species_segment_df["gene"], species_segment_df["nt_sequence_no_gaps"]):
                    joined_seqs.append(SeqRecord(Seq(seq), name=j, id=j))
                fasta_file = write_out_fasta(joined_seqs, species_blast_dir, species, "J")
                write_blast_db(fasta_file, fasta_file.split(".fasta")[0])