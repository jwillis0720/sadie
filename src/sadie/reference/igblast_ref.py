import logging
import os

# third party
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# module
from .blast import write_blast_db
from ..reference import get_loaded_database
from yaml import load as yml_load

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

logger = logging.getLogger(__file__)


def write_out_fasta(sequences, outpath):
    logger = logging.getLogger(__file__)
    output_fasta = outpath + ".fasta"
    logger.debug("output fasta {}".format(output_fasta))

    # Im sure this will come back to haunt me, but if, we've seen the name, save that
    seen_short_names = []
    with open(output_fasta, "w") as f:
        for sequence in sequences:
            name = sequence.name
            if name in seen_short_names:
                logger.warning(
                    "DANGER!!! DANGER!! We have already seen this gene name once, probably a duplicate species {}!!".format(
                        name
                    )
                )
                continue
            seq = str(sequence.seq).replace(".", "")
            f.write(">{}\n{}\n".format(name, seq))
            seen_short_names.append(name)

    return output_fasta


def get_databases_types(database_json):
    return list(set(map(lambda x: x["source"], database_json)))


def get_species_from_database(database_json):
    return list(set(map(lambda x: x["common"], database_json)))


def make_igblast_ref_database(reference, outpath):
    """[summary]

    Parameters
    ----------
    """
    # The blast DB groups by V,D and J
    logger.debug("Generating from IMGT Internal Database File")
    reference = yml_load(open(reference), Loader=Loader)
    database = get_loaded_database()

    for db_type in reference:
        for common in reference[db_type]:
            receptor_blast_dir = os.path.join(outpath, f"{db_type}/Ig/blastdb/")
            if not os.path.exists(receptor_blast_dir):
                os.makedirs(receptor_blast_dir)
            for segment in ["V", "D", "J"]:
                sub_species_keys = reference[db_type][common]
                if len(sub_species_keys) > 1:
                    chimera = True
                else:
                    chimera = False
                actual_genes = []
                for sub_species in sub_species_keys:
                    genes = reference[db_type][common][sub_species]
                    requested_genes_segments = list(filter(lambda x: x[3] == segment, genes))
                    common_data = list(filter(lambda x: x["common"] == sub_species, database[db_type]))
                    actual_genes += list(filter(lambda x: x["gene"] in requested_genes_segments, common_data))

                if not actual_genes:
                    logger.warning(f"Nothing found for {common}-{segment}-{database}")
                    continue
                out_segment = os.path.join(receptor_blast_dir, f"{common}_{segment}")
                if chimera:
                    seqs = list(
                        map(
                            lambda x: (x["common"] + "|" + x["gene"], x["sequence"]),
                            actual_genes,
                        )
                    )
                else:
                    seqs = list(map(lambda x: (x["gene"], x["sequence"]), actual_genes))
                joined_seqs = [SeqRecord(Seq(seq), name=n, id=n) for n, seq in seqs]
                fasta_file = write_out_fasta(joined_seqs, out_segment)
                write_blast_db(fasta_file, fasta_file.split(".fasta")[0])
                logger.info(f"Wrote blast for {fasta_file}")
