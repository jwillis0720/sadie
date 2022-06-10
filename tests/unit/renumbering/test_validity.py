import json

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
# def test_alignments_file(fixture_setup):

#     renumbering_api = Renumbering(
#         use_numbering_hmms=True,
#         allowed_species=["alpaca", "human", "mouse", "pig", "rabbit", "rat", "rhesus"],
#         run_multiproc=False
#     )

#     with open(fixture_setup.alignment_data / "anarci-alignments.json") as f:
#         alignments = json.load(f)

#     renumbering_results = renumbering_api.run_file(fixture_setup.alignment_data / "anarci-ali.fasta")
#     # from IPython import embed; embed()
#     # print(renumbering_results.get_alignment_table().head(1).to_dict('records'))
#     renumbering_alignments_obj = renumbering_results.get_alignment_table().fillna('').drop(['chain_type', 'scheme'], axis=1).to_dict("records")

#     renumbering_alignments = {}
#     for renumbering_align in renumbering_alignments_obj:
#         seq_id = renumbering_align.pop('Id')
#         renumbering_align = renumbering_align.items()
#         renumbering_alignments[seq_id] = renumbering_align

#     for seq_id, align_seq in alignments.items():
#         anarci_align = align_seq['align']
#         # seq = align_seq['seq']

#         renumbering_align = renumbering_alignments.get(seq_id)
#         # for debugging
#         if not renumbering_align:
#             continue

#         print(seq_id)

#         anarci_p = "".join([position for position, _ in anarci_align])
#         anarci_s = "".join([seq for _, seq in anarci_align])

#         renumber_p = "".join([position for position, _ in renumbering_align])
#         renumber_s = "".join([seq for _, seq in renumbering_align])

#         print(anarci_s)
#         print(renumber_s)

#         # dataframe will append '-' to the end of the sequence so we have to trim them off
#         assert anarci_s == renumber_s[:len(anarci_s)]
#         assert anarci_p == renumber_p[:len(anarci_p)]
