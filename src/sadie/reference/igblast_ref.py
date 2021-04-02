import logging
import os
import json
import gzip

# third party
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# module
from .blast import write_blast_db
from .yaml import YamlRef

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


def get_filtered_data(database_json, source, common, receptor, segment, subset):
    if subset == "all":
        function = ["ORF", "F", "P", "I"]
    elif subset == "functional":
        function = ["F"]
    else:
        raise Exception(f"{subset} not a valide option, ['all', 'funtional']")
    return list(
        filter(
            lambda x: x["source"] == source
            and x["common"] == common
            and x["gene_segment"] == segment
            and x["receptor"] == receptor
            and x["functional"] in function,
            database_json,
        )
    )


def make_igblast_ref_database(database, outdir):
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
    ig_database = json.load(gzip.open(database, "rt"))
    reference_database = YamlRef()

    for database in reference_database.yaml:
        for functional in reference_database.yaml[database]:
            for common in reference_database.yaml[database][functional]:
                receptor_blast_dir = os.path.join(outdir, f"{database}/{functional}/Ig/blastdb/")
                for segment in ["V", "D", "J"]:
                    sub_species_keys = reference_database.get_sub_species(database, functional, common)
                    if len(sub_species_keys) > 1:
                        chimera = True
                        actual_genes = []
                        for sub_species in reference_database.yaml[database][functional][common]:
                            genes = reference_database.get_gene_segment(
                                database, functional, common, sub_species, segment
                            )
                            matched = get_filtered_data(ig_database, database, sub_species, "Ig", segment, functional)

                            # filter again
                            actual_genes += list(filter(lambda x: x["gene"] in genes, matched))
                    else:
                        chimera = False
                        genes = reference_database.get_gene_segment(database, functional, common, common, segment)
                        matched = get_filtered_data(ig_database, database, common, "Ig", segment, functional)

                        # filter again
                        actual_genes = list(filter(lambda x: x["gene"] in genes, matched))

                    if not actual_genes:
                        logger.warning(f"Nothing found for {common}-{segment}-{database}")
                        continue
                    out_segment = os.path.join(receptor_blast_dir, f"{common}_{segment}")
                    if chimera:
                        seqs = list(
                            map(
                                lambda x: (x["common"] + "|" + x["gene"], x["imgt"][f"{segment.lower()}_gene_nt"]),
                                actual_genes,
                            )
                        )
                    else:
                        seqs = list(map(lambda x: (x["gene"], x["imgt"][f"{segment.lower()}_gene_nt"]), actual_genes))
                    joined_seqs = [SeqRecord(Seq(seq), name=n, id=n) for n, seq in seqs]
                    fasta_file = write_out_fasta(joined_seqs, out_segment)
                    write_blast_db(fasta_file, fasta_file.split(".fasta")[0])
                    logger.info(f"Wrote blast for {fasta_file}")

    # for receptor, subset, common, source in itertools.product(
    #     ["Ig", "TCR"],
    #     ["all", "functional"],
    #     get_species_from_database(ig_database),
    #     get_databases_types(ig_database),
    # ):
    #     receptor_blast_dir = os.path.join(outdir, f"{source}/{subset}/{receptor}/blastdb/")
    #     if not os.path.exists(receptor_blast_dir):
    #         logger.info("Creating %s", receptor_blast_dir)
    #         os.makedirs(receptor_blast_dir)
    #     for segment in list("VDJ"):
    #         filtered_df = get_filtered_data(ig_database, source, common, receptor, segment, subset)
    #         if not filtered_df:
    #             logger.info(f"No entries for {source}-{subset}-{receptor}-{common}-{segment}")
    #             continue
    #         out_segment = os.path.join(receptor_blast_dir, f"{common}_{segment}")
    #         genes = list(map(lambda x: x["gene"], filtered_df))
    #         seqs = list(map(lambda x: x["imgt"][f"{segment.lower()}_gene_nt"], filtered_df))
    #         joined_seqs = [SeqRecord(Seq(seq), name=v, id=v) for v, seq in zip(genes, seqs)]
