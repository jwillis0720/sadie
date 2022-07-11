from pathlib import Path
from sadie.utility.util import get_consensus_of_paired_end_abi, split_fasta
from tests.conftest import SadieFixture


def test_abi_consensus_check(fixture_setup: SadieFixture):
    pe_files = fixture_setup.get_abi_pe_files()
    abi_file_1_path = pe_files[0]
    abi_file_2_path = pe_files[1]
    assert all([abi_file_1_path.exists(), abi_file_2_path.exists()])
    consensus_seq = get_consensus_of_paired_end_abi(abi_file_1_path, abi_file_2_path)
    assert consensus_seq == fixture_setup.get_pe_consensus_seq()


def test_fasta_splitter(fixture_setup: SadieFixture):
    test_file: Path = fixture_setup.get_oas_fasta()
    tempdir = Path(fixture_setup.tmp_path)
    split_fasta(test_file, 20, tempdir / Path("split_fasta"), "fasta")
