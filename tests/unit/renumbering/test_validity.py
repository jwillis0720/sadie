import json

import pandas as pd
from Bio import SeqIO

from sadie.renumbering import Renumbering


def test_alignments(fixture_setup):
    renumbering_api = Renumbering(
        use_numbering_hmms=True,
        allowed_species=["human"],
    )

    with open(fixture_setup.alignment_data / "catnap_aa_heavy_sample_mutually-exclusive.json") as f:
        anarci_alignments = json.load(f)

    for seqrecord in SeqIO.parse(fixture_setup.fasta_inputs / "catnap_aa_heavy_sample.fasta", format="fasta"):
        seq_id = seqrecord.id
        seq = str(seqrecord.seq)

        anarci_align = anarci_alignments[seq_id]

        renumber_align = [
            list(v)
            for v in list(
                renumbering_api.run_single(seq_id="0", seq=seq)
                .get_alignment_table()
                .iloc[:, 3:]
                .to_dict("records")[0]
                .items()
            )
        ]

        # print(renumbering_api.run_single(seq_id="0", seq=seq).to_dict('records'))

        anarci_p = "".join([position for position, _ in anarci_align])
        anarci_s = "".join([seq for _, seq in anarci_align])

        renumber_p = "".join([position for position, _ in renumber_align])
        renumber_s = "".join([seq for _, seq in renumber_align])

        # print(seq_id)
        # print(seq)
        # print(anarci_s)
        # print(renumber_s)

        assert anarci_s == renumber_s
        assert anarci_p == renumber_p


# TODO: the alignment changes with mutiple need to have a seperate ANARCI align to a complete fasta
def test_alignments_file(fixture_setup):
    anarci_table = pd.read_csv(fixture_setup.alignment_data / "catnap_aa_heave_sample_H.csv")
    anarci_table["Id"] = anarci_table["Id"].str.replace("<unknown description>", "").str.strip()
    anarci_table = anarci_table[anarci_table.domain_no == 0].iloc[:, [0, *list(range(13, len(anarci_table.columns)))]]

    renumbering_api = Renumbering(
        use_numbering_hmms=True,
        allowed_species=["human"],
        run_multiproc=False,
    )

    renumbering_results = renumbering_api.run_file(fixture_setup.fasta_inputs / "catnap_aa_heavy_sample.fasta")
    renumbering_table = renumbering_results.get_alignment_table()
    renumbering_table = renumbering_table.iloc[:, [0, *list(range(3, len(renumbering_table.columns)))]]

    anarci_table = anarci_table.sort_values("Id").reset_index(drop=True)
    renumbering_table = renumbering_table.sort_values("Id").reset_index(drop=True)

    pd.testing.assert_frame_equal(
        anarci_table, renumbering_table, check_dtype=False, check_index_type=False, check_frame_type=False
    )
