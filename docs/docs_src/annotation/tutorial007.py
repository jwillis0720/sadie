import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from sadie.airr import Airr


class GenBankProcessor:
    def __init__(self, email, gene_id, airr_species="human"):
        Entrez.email = email
        self.gene_id = gene_id
        self.airr_api = Airr(airr_species)

    def fetch_genebank_record(self):
        try:
            handle = Entrez.efetch(db="nucleotide", id=self.gene_id, rettype="gb")
            genebank_record = list(SeqIO.parse(handle, "genbank"))[0]
            return genebank_record
        except Exception as e:
            print(f"Error fetching GenBank record: {e}")
            return None

    def get_airr_table(self, gene_symbol, sequence):
        return self.airr_api.run_single(gene_symbol, sequence)


class GenBankFeatureAdder:
    def __init__(self, genebank_record, map_dict):
        self.genebank_record = genebank_record
        self.map_dict = map_dict

    def add_feature(self, feature, airr_table):
        feature_type = feature
        if feature in ["VGene", "JGene", "DGene"]:
            start = int(airr_table[f"{feature_type.lower()[0]}_sequence_start"][0]) - 1
            end = int(airr_table[f"{feature_type.lower()[0]}_sequence_end"][0])
            qualifier_dict = {
                "gene": [airr_table[self.map_dict[feature_type]][0]],
                "species": [airr_table["reference_name"][0]],
            }
        else:
            start = int(airr_table[f"{feature_type.lower()}_start"][0]) - 1
            end = int(airr_table[f"{feature_type.lower()}_end"][0])
            qualifier_dict = {
                "gene": [airr_table[self.map_dict[feature_type]][0]],
                "reference": [airr_table["reference_name"][0]],
            }

        location = FeatureLocation(start, end)
        _feature = SeqFeature(location, type=feature_type, qualifiers=qualifier_dict)
        self.genebank_record.features.append(_feature)


def get_record_from_seq(seq, seq_id, seq_name):
    seq = Seq(seq)
    record = SeqRecord(id=seq_id, seq=seq, name=seq_name)
    record.annotations["molecule_type"] = "DNA"
    return record


# genbank_record = get_record_from_seq(seq, seq_id, seq_name)


def main(email=None, gene_id=None, seq=None, seq_id=None, seq_name=None):
    map_dict = {
        "FWR1": "v_call",
        "FWR2": "v_call",
        "FWR3": "v_call",
        "FWR4": "j_call",
        "CDR1": "v_call",
        "CDR2": "v_call",
        "CDR3": "d_call",
        "VGene": "v_call",
        "JGene": "j_call",
        "DGene": "d_call",
    }

    if seq is not None:
        processor = GenBankProcessor(email, gene_id)
        genbank_record = get_record_from_seq(seq, seq_id, seq_name)
        # print(genbank_record.seq)
    else:
        processor = GenBankProcessor(email, gene_id)
        genbank_record = processor.fetch_genebank_record()
        # print('Also got here')

    feature_types = ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "CDR3", "FWR4", "VGene", "DGene", "JGene"]

    for feature in feature_types:
        airr_table = processor.get_airr_table(feature, genbank_record.seq)
        feature_adder = GenBankFeatureAdder(genbank_record, map_dict)
        feature_adder.add_feature(feature, airr_table)
    print(genbank_record.format("genbank"))
    return genbank_record


if __name__ == "__main__":
    email = "example@mail.com"
    gene_id = "GU272045.1"
    pg9_seq = "CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGTCCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGTGGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCCAGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATTTTTGTGTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTATTATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC"
    id = "PG9-Antibody"
    seq_name = "PG9"

    # Fetch the genbank file and annotate
    genebank_record = main(seq=pg9_seq, seq_id=id, seq_name=seq_name)
    with open(f"{id}_airr.gb", "w") as handle:
        SeqIO.write(genebank_record, handle, "genbank")
