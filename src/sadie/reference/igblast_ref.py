import logging
import os
from pathlib import Path
from typing import List

# third party
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sadie.reference import Reference

# module
from .blast import write_blast_db

logger = logging.getLogger(__file__)


def write_out_fasta(sequences: List[SeqRecord], outpath: Path) -> Path:
    logger = logging.getLogger(__file__)
    output_fasta = outpath + ".fasta"
    logger.debug("output fasta {}".format(output_fasta))

    # Im sure this will come back to haunt me, but if, we've seen the name, save that
    seen_short_names = []
    with open(output_fasta, "w") as f:
        for sequence in sequences:
            name = sequence.name
            seq = str(sequence.seq).replace(".", "")
            f.write(">{}\n{}\n".format(name, seq))
            seen_short_names.append(name)
    return output_fasta


def get_databases_types(database_json):
    return list(set(map(lambda x: x["source"], database_json)))


def get_species_from_database(database_json):
    return list(set(map(lambda x: x["common"], database_json)))


def make_igblast_ref_database(reference: Reference, outpath: Path):
    # The blast DB groups by V,D and J
    logger.debug("Generating from IMGT Internal Database File")
    database = reference.get_dataframe()
    for group, group_df in database.groupby("source"):
        for species, species_df in group_df.groupby("species"):
            receptor_blast_dir = os.path.join(outpath, f"{group}/Ig/blastdb/")
            sub_species_keys = species_df["sub_species"].unique()
            if not os.path.exists(receptor_blast_dir):
                os.makedirs(receptor_blast_dir)
            for segment, segment_df in species_df.groupby("gene_segment"):
                if len(sub_species_keys) > 1:
                    chimera = True
                else:
                    chimera = False

                # for sub_species in sub_species_keys:
                #     genes = reference[db_type][common][sub_species]
                #     requested_genes_segments = list(filter(lambda x: x[3] == segment, genes))
                #     common_data = list(filter(lambda x: x["common"] == sub_species, database[db_type]))
                #     actual_genes += list(filter(lambda x: x["gene"] in requested_genes_segments, common_data))

                # if not actual_genes:
                #     logger.warning(f"Nothing found for {species}-{segment}-{database}")
                #     continue
                out_segment = os.path.join(receptor_blast_dir, f"{species}_{segment}")
                if chimera:
                    seqs = segment_df.apply(
                        lambda x: SeqRecord(Seq(str(x["sequence"])), name=x["sub_species"] + "|" + x["gene"]), axis=1
                    ).to_list()
                else:
                    seqs = segment_df.apply(
                        lambda x: SeqRecord(Seq(str(x["sequence"])), name=x["gene"]), axis=1
                    ).to_list()
                # returns a full fasta path
                fasta_file = write_out_fasta(seqs, out_segment)
                write_blast_db(fasta_file, fasta_file.split(".fasta")[0])
                logger.info(f"Wrote blast for {fasta_file}")
