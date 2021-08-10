from sadie.utility.util import get_consensus_of_paired_end_abi


def test_abi_consensus_check(fixture_setup):
    pe_files = fixture_setup.get_abi_pe_files()
    abi_file_1_path = pe_files[0]
    abi_file_2_path = pe_files[1]
    assert all([abi_file_1_path.exists(), abi_file_2_path.exists()])
    consensus_seq = get_consensus_of_paired_end_abi(abi_file_1_path, abi_file_2_path)
    assert consensus_seq == fixture_setup.get_pe_consensus_seq()
