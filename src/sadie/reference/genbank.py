import datetime
import itertools
import gzip
import json
import logging
import os

from Bio import Seq, SeqIO, SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature

# Module level logger
logger = logging.getLogger(__name__)


def get_databases_types(database_json):
    return list(set(map(lambda x: x["source"], database_json)))


def get_species_from_database(database_json):
    return list(set(map(lambda x: x["common"], database_json)))


def get_filtered_data(database_json, source, common, receptor, segment):
    return list(
        filter(
            lambda x: x["source"] == source
            and x["common"] == common
            and x["gene_segment"] == segment
            and x["receptor"] == receptor,
            database_json,
        )
    )


def generate_gb_record(sequence, name, description):

    # Create a sequence
    # sequence_string = str(sequence)
    sequence_object = Seq.Seq(str(sequence))

    # Create a record
    record = SeqRecord.SeqRecord(sequence_object, id=name, name=name, description=description)
    return record


def generate_v(filtered_df):
    # all SeqRecords
    records = []
    for entry in filtered_df:
        # A signle json entry
        segment_dictionary = entry["imgt"]
        gene = entry["gene"]
        source = entry["source"]
        common = entry["common"]
        functional = entry["functional"]
        latin = entry["latin"]
        sequence = segment_dictionary["v_gene_nt"]
        aa_sequecne = segment_dictionary["v_gene_aa"]
        name = "|".join([common, source, gene])
        v_gene_record = generate_gb_record(sequence, name, latin)

        for segment in ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3"]:
            if not segment_dictionary[f"{segment}_nt"]:
                logger.debug(f"{common}-{gene}-{segment} nt is empty")
                continue
            if segment == "cdr3":
                segment_start = segment_dictionary["fwr3_nt_index_end"] + 1
                segment_end = len(sequence)

            else:
                segment_start = segment_dictionary[f"{segment}_nt_index_start"]
                segment_end = segment_dictionary[f"{segment}_nt_index_end"]
            translation = segment_dictionary[f"{segment}_aa"]
            if segment[0:2] == "fw":
                seg_type = "FR"
            else:
                seg_type = "CDR"
            feature = SeqFeature(
                FeatureLocation(start=segment_start - 1, end=segment_end),
                type=seg_type,
                qualifiers={
                    "gene": gene,
                    "latin": latin,
                    "organism": common,
                    "functional": functional,
                    "scheme": "imgt",
                    "source": source,
                    "translation": translation,
                    "standard_name": segment.upper().replace("W", ""),
                },
            )
            v_gene_record.features.append(feature)

        # annotate the whole v gene
        v_gene_record.annotations = {
            "organism": common,
            "latin": latin,
            "common": common,
            "source": source,
            "molecule_type": "nucleotide",
            "date": (datetime.date.today().strftime("%d-%b-%Y")),
        }

        v_gene_seq_feature = SeqFeature(
            FeatureLocation(start=0, end=len(sequence)),
            type="v_segment",
            qualifiers={
                "gene": gene,
                "latin": latin,
                "organism": common,
                "functional": functional,
                "scheme": "imgt",
                "source": source,
                "translation": aa_sequecne,
                "standard_name": f"{source}|{common}|{gene}",
            },
        )
        imgt_numbering = segment_dictionary["v_gene_imgt_numbering"]
        for index in range(0, len(imgt_numbering)):
            number = index * 3
            imgt_feature = SeqFeature(
                FeatureLocation(start=number, end=number + 3),
                type="IMGT",
                qualifiers={"standard_name": imgt_numbering[index]},
            )
            v_gene_record.features.append(imgt_feature)
        # annotate the vgene itself
        v_gene_record.features.append(v_gene_seq_feature)
        records.append(v_gene_record)
    return records


def generate_j(filtered_df):
    # all SeqRecords
    records = []
    for entry in filtered_df:
        # A signle json entry
        gene = entry["gene"]
        source = entry["source"]
        common = entry["common"]
        functional = entry["functional"]
        latin = entry["latin"]
        name = "|".join([common, source, gene])

        # Segment dictionary
        segment_dictionary = entry["imgt"]
        j_gene_nt = segment_dictionary["j_gene_nt"]
        # take_a_look = j_gene_nt
        j_gene_aa = segment_dictionary["j_gene_aa"]
        # cdr3_nt = segment_dictionary["cdr3_nt"]
        cdr3_aa = segment_dictionary["cdr3_aa"]
        fwr4_nt = segment_dictionary["fwr4_nt"]
        fwr4_aa = segment_dictionary["fwr4_aa"]
        reading_frame = segment_dictionary["reading_frame"]
        end_cdr3_nt = segment_dictionary["end_cdr3_nt"]

        exonucleate_this = len(j_gene_nt[int(reading_frame) - 1 :]) % 3
        if exonucleate_this:
            j_gene_nt = j_gene_nt[:-exonucleate_this]
            fwr4_nt = fwr4_nt[:-exonucleate_this]

        # J gene record

        j_gene_record = generate_gb_record(j_gene_nt, name, latin)
        if end_cdr3_nt:
            # CDR3
            cdr3_nt_feature = SeqFeature(
                FeatureLocation(start=0, end=end_cdr3_nt),
                type="CDR",
                qualifiers={
                    "gene": gene,
                    "latin": latin,
                    "organism": common,
                    "functional": functional,
                    "scheme": "imgt",
                    "source": source,
                    "translation": cdr3_aa,
                    "standard_name": "CDR3",
                },
            )
            j_gene_record.features.append(cdr3_nt_feature)

            # FW4
            fwr4_nt_feature = SeqFeature(
                FeatureLocation(start=end_cdr3_nt, end=len(j_gene_nt)),
                type="FR",
                qualifiers={
                    "gene": gene,
                    "latin": latin,
                    "organism": common,
                    "functional": functional,
                    "scheme": "imgt",
                    "source": source,
                    "translation": fwr4_aa,
                    "standard_name": "FR4",
                },
            )
            j_gene_record.features.append(fwr4_nt_feature)

        # annotate the whole j gene
        j_gene_record.annotations = {
            "organism": common,
            "latin": latin,
            "common": common,
            "source": source,
            "molecule_type": "nucleotide",
            "date": (datetime.date.today().strftime("%d-%b-%Y")),
        }
        j_gene_seq_feature = SeqFeature(
            FeatureLocation(start=0, end=len(j_gene_nt)),
            type="j_segment",
            qualifiers={
                "gene": gene,
                "latin": latin,
                "organism": common,
                "functional": functional,
                "scheme": "imgt",
                "source": source,
                "translation": j_gene_aa,
                "standard_name": f"{source}|{common}|{gene}",
            },
        )
        j_gene_record.features.append(j_gene_seq_feature)
        fwr4_index_start = j_gene_nt.index(fwr4_nt)
        imgt_start = 118
        for index in range(fwr4_index_start, fwr4_index_start + len(fwr4_nt), 3):
            imgt_feature = SeqFeature(
                FeatureLocation(start=index, end=index + 3),
                type="IMGT",
                qualifiers={"standard_name": imgt_start},
            )
            imgt_start += 1
            j_gene_record.features.append(imgt_feature)
        records.append(j_gene_record)
    return records


def generate_d(filtered_df):
    # all SeqRecords
    records = []
    for entry in filtered_df:
        # A signle json entry
        gene = entry["gene"]
        source = entry["source"]
        common = entry["common"]
        functional = entry["functional"]
        latin = entry["latin"]
        name = "|".join([common, source, gene])
        d_seq = entry["imgt"]["d_gene_nt"]
        d_gene_record = generate_gb_record(d_seq, name, latin)
        d_gene_seq_feature = SeqFeature(
            FeatureLocation(start=0, end=len(d_seq)),
            type="d_segment",
            qualifiers={
                "gene": gene,
                "latin": latin,
                "organism": common,
                "functional": functional,
                "scheme": "imgt",
                "source": source,
                "standard_name": f"{source}|{common}|{gene}",
            },
        )
        d_gene_record.annotations = {
            "organism": common,
            "latin": latin,
            "common": common,
            "source": source,
            "molecule_type": "nucleotide",
            "date": (datetime.date.today().strftime("%d-%b-%Y")),
        }
        d_gene_record.features.append(d_gene_seq_feature)
        records.append(d_gene_record)
    return records


def generate_genbank(database, outdir):
    """Generate genbank files from your database

    Arguments:
        database {file} -- where is the database file
        outpath {path} -- where to output
    """
    ig_database = json.load(gzip.open(database, "rt"))
    logger = logging.getLogger("Genebank")
    for receptor, common, source in itertools.product(
        ["Ig", "TCR"],
        get_species_from_database(ig_database),
        get_databases_types(ig_database),
    ):
        receptor_gb_dir = os.path.join(outdir, f"{source}/{receptor}")
        if not os.path.exists(receptor_gb_dir):
            logger.info(f"Making directory {receptor_gb_dir}")
            os.makedirs(receptor_gb_dir)
        for segment in list("VDJ"):
            out_file = f"{receptor_gb_dir}/{common}_{segment}.gb"
            logger.info(f"Writing files {out_file}")
            filtered_df = get_filtered_data(ig_database, source, common, receptor, segment)
            if not filtered_df:
                continue

            if segment == "V":
                records = generate_v(filtered_df)
                SeqIO.write(records, open(out_file, "w"), "gb")
            if segment == "J":
                records = generate_j(filtered_df)
                SeqIO.write(records, open(out_file, "w"), "gb")
            if segment == "D":
                records = generate_d(filtered_df)
                SeqIO.write(records, open(out_file, "w"), "gb")
            logger.info(f"Wrote files {out_file}")
    return outdir
