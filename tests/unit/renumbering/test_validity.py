import json

import pandas as pd

from sadie.renumbering import Renumbering


def test_alignments(fixture_setup):

    renumbering_api = Renumbering(
        use_numbering_hmms=True,
        allowed_species=["alpaca", "human", "mouse", "pig", "rabbit", "rat", "rhesus"],
    )

    with open(fixture_setup.alignment_data / "anarci-alignments.json") as f:
        alignments = json.load(f)

    for _seq_id, align_seq in alignments.items():
        anarci_align = align_seq["align"]
        seq = align_seq["seq"]

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

        anarci_p = "".join([position for position, _ in anarci_align])
        anarci_s = "".join([seq for _, seq in anarci_align])

        renumber_p = "".join([position for position, _ in renumber_align])
        renumber_s = "".join([seq for _, seq in renumber_align])

        # print(anarci_s)
        # print(renumber_s)

        assert anarci_s == renumber_s
        assert anarci_p == renumber_p


# TODO: the alignment changes with mutiple need to have a seperate ANARCI align to a complete fasta
def test_alignments_file(fixture_setup):

    anarci_table = pd.read_csv(fixture_setup.alignment_data / "anarci-ali-complete.csv")
    anarci_table = anarci_table[anarci_table.domain_no == 0].iloc[:, [0, *list(range(13, len(anarci_table.columns)))]]

    renumbering_api = Renumbering(
        use_numbering_hmms=True,
        allowed_species=["alpaca", "human", "mouse", "pig", "rabbit", "rat", "rhesus"],
        run_multiproc=False,
    )

    renumbering_results = renumbering_api.run_file(fixture_setup.alignment_data / "anarci-ali.fasta")
    renumbering_table = renumbering_results.get_alignment_table()
    renumbering_table = renumbering_table.iloc[:, [0, *list(range(3, len(renumbering_table.columns)))]]

    anarci_table = anarci_table.sort_values("Id").reset_index(drop=True)
    renumbering_table = renumbering_table.sort_values("Id").reset_index(drop=True)

    pd.testing.assert_frame_equal(
        anarci_table, renumbering_table, check_dtype=False, check_index_type=False, check_frame_type=False
    )
