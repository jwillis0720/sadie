import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from sadie.renumbering import Renumbering


class ProteinSequenceProcessor:
    def __init__(
        self,
        sequence,
        seq_name,
        molecule_type="protein",
        organism="Human",
        scheme="imgt",
        region_assign="imgt",
        run_multiproc=False,
    ):
        self.sequence = sequence
        self.scheme = scheme
        self.organism = organism
        self.name = seq_name
        self.type = molecule_type
        self.region_assign = region_assign
        self.run_multiproc = run_multiproc
        self.map_number = {
            "FWR1": "v_gene",
            "CDR1": "v_gene",
            "FWR2": "v_gene",
            "CDR2": "v_gene",
            "FWR3": "v_gene",
            "CDR3": "d_gene",
            "FWR4": "j_gene",
            "VGene": "v_gene",
            "DGene": "d_gene",
            "JGene": "j_gene",
        }

    def process_sequence(self):
        renumbering_api = Renumbering(
            scheme=self.scheme, region_assign=self.region_assign, run_multiproc=self.run_multiproc
        )

        # Run sequence and return renumbering table with sequence_id and sequence
        numbering_table = renumbering_api.run_single(self.name, self.sequence)

        seq = Seq(self.sequence)
        record = SeqRecord(id=self.name, seq=seq, name=self.name)
        record.annotations["molecule_type"] = self.type
        record.annotations["organism"] = numbering_table.hmm_species[0]

        for feature in self.map_number.keys():
            _feature = self._get_feature_numbering(numbering_table, feature)
            if _feature:
                record.features.append(_feature)

        return record

    def _get_start_stop(self, feature, numbering_table):
        start = self.sequence.find(numbering_table[f"{feature.lower()}_aa_no_gaps"][0])
        end = start + len(numbering_table[f"{feature.lower()}_aa_no_gaps"][0])
        return start, end

    def _get_feature_numbering(self, numbering_table, feature):
        feature_type = feature
        if feature in ["VGene", "JGene", "DGene"]:
            return None  # The numbering scheme is missing these details (will try and get FW and CDR later)
        else:
            start, end = self._get_start_stop(feature, numbering_table)
            qualifier_dict = {"reference": numbering_table.hmm_species[0]}
            try:
                qualifier_dict["gene"] = [numbering_table[self.map_number[feature_type]][0]]
            except KeyError:
                qualifier_dict["gene"] = "Missing"
            location = FeatureLocation(start, end)
            _feature = SeqFeature(location, type=feature_type, qualifiers=qualifier_dict)

        return _feature


if __name__ == "__main__":
    seg_name = "PG9"
    pg9_aa = "QRLVESGGGVVQPGSSLRLSCAASGFDFSRQGMHWVRQAPGQGLEWVAFIKYDGSEKYHADSVWGRLSISRDNSKDTLYLQMNSLRVEDTATYFCVREAGGPDYRNGYNYYDFYDGYYNYHYMDVWGKGTTVTVSS"
    processor = ProteinSequenceProcessor(pg9_aa, seg_name)
    genbank_record = processor.process_sequence()
    print(genbank_record.format("genbank"))
