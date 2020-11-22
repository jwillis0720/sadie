import datetime
import glob
import logging
import os
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import Seq, SeqIO, SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature

# Module level
from .settings import IMGT_LOOKUP, REVERSE_IMGT_LOOKUP

# Module level logger
logger = logging.getLogger(__name__)


def generate_gb_record(sequence, name, description):

    # Create a sequence
    # sequence_string = str(sequence)
    sequence_object = Seq.Seq(str(sequence))

    # Create a record
    record = SeqRecord.SeqRecord(
        sequence_object,
        id=name,
        name=name,
        description=description,
    )
    return record


def generate_v(imgt_db, outpath):
    # VGene
    imgt_db_df = pd.read_csv(imgt_db, index_col=0)

    # records is a dictionary of lists
    for species, species_df in imgt_db_df.groupby("common"):
        species_df = species_df.set_index("gene")
        records = defaultdict(list)
        for gene in species_df.index:
            # lookup gene
            v_gene_series = species_df.loc[gene, ["FW1", "CDR1", "FW2", "CDR2", "FW3", "CDR3"]]

            # Is it functional
            functional = species_df.loc[gene, "functional"]

            # Latin name
            latin = species_df.loc[gene, "latin"]

            # if we just return a series
            if isinstance(v_gene_series, pd.Series):
                v_gene_seq = v_gene_series.str.cat()

            # if we reutrn a dataframe
            else:
                v_gene_series = v_gene_series.iloc[-1]
                v_gene_seq = v_gene_series.str.cat()
                functional = functional.iloc[-1]
                latin = latin.iloc[-1]

            functional = functional.replace("(", "")
            functional = functional.replace(")", "")
            # sequence, name description
            v_gene_record = generate_gb_record(v_gene_seq, gene, latin)

            # iterate through features
            for feature in v_gene_series.index:
                feature_string = v_gene_series[feature]
                if feature_string is np.nan:
                    continue

                if feature == "CDR3":
                    start = len(v_gene_seq) - len(feature_string)
                    end = len(v_gene_seq)
                else:
                    start = v_gene_seq.index(feature_string)
                    end = start + len(feature_string)

                qualifiers_dict = {
                    "gene": gene,
                    "latin": latin,
                    "organism": species,
                    "functional": functional,
                    "scheme": "IMGT",
                }
                ##annotate each part of the v gene
                v_gene_record.features.append(
                    SeqFeature(FeatureLocation(start=start, end=end), type=feature, qualifiers=qualifiers_dict)
                )
                ##But also generate the segment record
                segment_record_name = "{}|{}".format(gene, species)
                segment_record = generate_gb_record(feature_string, segment_record_name, latin)

                ##this is to catch a locus name that is too long
                if len(segment_record_name) > 16:
                    if len(segment_record_name) + 1 + len(str(len(segment_record))) > 28:
                        segment_record_name = "{}|{}".format(gene[0 : 15 - len(species)], species)
                        segment_record = generate_gb_record(feature_string, segment_record_name, latin)
                seq_feature = SeqFeature(
                    FeatureLocation(start=0, end=len(feature_string)), type=feature, qualifiers=qualifiers_dict
                )
                seq_feature.qualifiers["codon_start"] = [1]
                seq_feature.qualifiers["translation"] = Seq.Seq(feature_string).translate()
                segment_record.features.append(seq_feature)
                segment_record.annotations = {
                    "organism": latin,
                    "source": species,
                    "date": (datetime.date.today().strftime("%d-%b-%Y")),
                }
                ##Add to the segment dictionary
                records[feature].append(segment_record)
            # annotate the whole v gene
            v_gene_record.annotations = {
                "organism": latin,
                "source": species,
                "date": (datetime.date.today().strftime("%d-%b-%Y")),
            }

            v_gene_sq_feature = SeqFeature(
                FeatureLocation(start=0, end=len(v_gene_seq)), type="VGene", qualifiers=qualifiers_dict
            )

            v_gene_sq_feature.qualifiers["codon_start"] = [1]
            v_gene_sq_feature.qualifiers["translation"] = Seq.Seq(v_gene_seq).translate()
            ##annotate the vgene itself
            v_gene_record.features.append(v_gene_sq_feature)
            records["v_gene"].append(v_gene_record)

        _v_gene_path = os.path.join(outpath, "VGenes")
        if not os.path.exists(_v_gene_path):
            logger.info("Creating %s", _v_gene_path)
            os.makedirs(_v_gene_path)
        ###done V gene
        SeqIO.write(records["v_gene"], _v_gene_path + f"/{species}_V.gb", "genbank")
        SeqIO.write(records["FW1"], _v_gene_path + f"/{species}_FW1.gb", "genbank")
        SeqIO.write(records["FW2"], _v_gene_path + f"/{species}_FW2.gb", "genbank")
        SeqIO.write(records["FW3"], _v_gene_path + f"/{species}_FW3.gb", "genbank")
        SeqIO.write(records["CDR1"], _v_gene_path + f"/{species}_CDR1.gb", "genbank")
        SeqIO.write(records["CDR2"], _v_gene_path + f"/{species}_CDR2.gb", "genbank")


def generate_j(aux_path, imgt_fasta, outpath):

    for files in glob.glob(aux_path + "/*"):
        species = files.split("/")[-1].split("_")[0]
        j_records = []
        fw4_records = []
        imgt_lookup = IMGT_LOOKUP[species]
        j_files = glob.glob("{}/{}/*/*J.fasta".format(imgt_fasta, imgt_lookup))
        if not j_files:
            logger.critical("Cant find any J genes in %s fror %s", aux_path, species)
            continue
        lookup_table = pd.read_csv(files, delim_whitespace=True, names=["j", "reading_frame", "end_cdr3"])
        for fasta_file in j_files:
            for sequence in SeqIO.parse(fasta_file, "fasta"):
                short_name = sequence.id.split("|")[1]
                functional = sequence.description.split("|")[3]
                functional = functional.replace("(", "")
                sequence_strip = str(sequence.seq).replace(".", "")
                j_record = generate_gb_record(sequence_strip, short_name, sequence.description)
                try:
                    gene_lookup = lookup_table.loc[short_name]
                    if isinstance(gene_lookup, pd.DataFrame):
                        gene_lookup = gene_lookup.iloc[0]

                except KeyError:
                    logger.warning("didnt find {} for species {}".format(short_name, species))
                    continue

                cdr3_end = int(gene_lookup["end_cdr3"])
                j_record.features.append(
                    SeqFeature(
                        FeatureLocation(start=0, end=cdr3_end),
                        type="CDR3",
                        qualifiers={
                            "gene": short_name,
                            "organism": species,
                            "functional": functional,
                            "latin": imgt_lookup,
                        },
                    )
                )

                # FW4
                fw4_string = sequence_strip[cdr3_end + 1 :]
                fw4_record = generate_gb_record(fw4_string, f"FW4|{species}|{short_name}", sequence.description)
                fw4_record.features.append(
                    SeqFeature(
                        FeatureLocation(start=0, end=len(fw4_string)),
                        type="FW4",
                        qualifiers={
                            "gene": short_name,
                            "organism": species,
                            "functional": functional,
                            "latin": imgt_lookup,
                        },
                    )
                )

                j_record.features.append(
                    SeqFeature(
                        FeatureLocation(start=cdr3_end, end=len(sequence_strip)),
                        type="FW4",
                        qualifiers={
                            "gene": short_name,
                            "organism": species,
                            "functional": functional,
                            "latin": imgt_lookup,
                        },
                    )
                )

                j_record.features.append(
                    SeqFeature(
                        FeatureLocation(start=0, end=len(sequence_strip)),
                        type="JGene",
                        qualifiers={
                            "gene": short_name,
                            "organism": species,
                            "functional": functional,
                            "latin": imgt_lookup,
                        },
                    )
                )

                j_record.annotations = {
                    "organism": imgt_lookup,
                    "source": species,
                    "date": (datetime.date.today().strftime("%d-%b-%Y")),
                }
                fw4_record.annotations = {
                    "organism": imgt_lookup,
                    "source": species,
                    "date": (datetime.date.today().strftime("%d-%b-%Y")),
                }
                j_records.append(j_record)
                fw4_records.append(fw4_record)

        _j_gene_path = os.path.join(outpath, "JGenes")
        if not os.path.exists(_j_gene_path):
            logger.info("Creating %s", _j_gene_path)
            os.makedirs(_j_gene_path)
        SeqIO.write(j_records, _j_gene_path + f"/{species}_J.gb", "genbank")
        SeqIO.write(fw4_records, _j_gene_path + f"/{species}_FW4.gb", "genbank")


def generate_d(imgt_fasta, aux_path, outpath):
    logger.info("Attempting D Genes")
    for species_directory in glob.glob(imgt_fasta + "/*"):
        records = []
        # file directory is the latin name
        latin = os.path.basename(species_directory)

        common_name = REVERSE_IMGT_LOOKUP[latin]
        d_files = glob.glob("{}/{}/*/*D.fasta".format(imgt_fasta, latin))
        if not d_files:
            logger.critical("Cant find any d genes in %s fror %s", aux_path, latin)
            continue
        for fasta_file in d_files:
            for sequence in SeqIO.parse(fasta_file, "fasta"):
                short_name = sequence.id.split("|")[1]
                functional = sequence.description.split("|")[3]
                functional = functional.replace("(", "")
                functional = functional.replace(")", "")
                record_name = functional + "|" + short_name
                sequence_strip = str(sequence.seq).replace(".", "")
                record = generate_gb_record(sequence_strip, record_name, sequence.description)
                record.features.append(
                    SeqFeature(
                        FeatureLocation(start=0, end=len(sequence.seq)),
                        type="DGene",
                        qualifiers={
                            "gene": short_name,
                            "organism": common_name,
                            "latin": latin,
                            "functional": functional,
                        },
                    )
                )

                record.annotations = {
                    "organism": latin,
                    "source": common_name,
                    "date": (datetime.date.today().strftime("%d-%b-%Y")),
                }
                records.append(record)

        _d_gene_path = os.path.join(outpath, "DGenes")
        if not os.path.exists(_d_gene_path):
            logger.info("Creating %s", _d_gene_path)
            os.makedirs(_d_gene_path)
        SeqIO.write(records, _d_gene_path + f"/{common_name}_D.gb", "genbank")


def generate_genbank(imgt_fasta, imgt_db, aux_path, outpath):
    """Generate genbank files from your database

    Arguments:
        imgt_fasta {path} -- path to iMGT vquest fasta files
        imgt_db {file} -- where is the IMGT DataBase File
        outpath {path} -- where to output
    """
    generate_v(imgt_db, outpath)
    generate_j(aux_path, imgt_fasta, outpath)
    generate_d(imgt_fasta, aux_path, outpath)
    # do d Genes
