# def test_anarci_multiprocess():
#     main_file = fixture_file("scfv_main_file.fasta")
#     anarci_api = Anarci(allowed_species=["dog", "cat"])
#     results = anarci_api.run_file(main_file, multi=True)
#     assert isinstance(results, AnarciResults)
#     other_results = AnarciResults.read_json(fixture_file("mp_result.json.gz"))
#     assert_frame_equal(other_results.segment_table.sort_index(), results.segment_table.sort_index())
#     assert_frame_equal(other_results.summary_table.sort_index(), results.summary_table.sort_index())
#     assert_frame_equal(other_results.alignment_table.sort_index(), results.alignment_table.sort_index())
