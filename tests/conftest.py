import bz2
import gzip
import shutil
from pathlib import Path
from typing import List
import pytest
from Bio import SeqIO

from sadie.airr import Airr
from sadie.airr.airrtable import AirrTable


class SadieFixture:
    def __init__(self, tmp_path_factory):
        self.base_datadir = Path("tests/data/fixtures/")
        self.fasta_inputs = self.base_datadir / "fasta_inputs/"
        self.fastq_inputs = self.base_datadir / "fastq_inputs/"
        self.airr_table_inputs = self.base_datadir / "airr_tables/"
        self.tmp_path = tmp_path_factory.mktemp("sadie_fixture")

    def get_card(self) -> Path:
        return self.base_datadir / "card.png"

    def get_dummy_scfv_table(self) -> Path:
        return self.airr_table_inputs / "linked_airr_table_dummy.csv.gz"

    def get_dog_airrdataframe_file(self) -> Path:
        return self.airr_table_inputs / "dog_igh.csv.gz"

    def get_json_as_dataframe(self) -> Path:
        return self.airr_table_inputs / "airrtable.json.gz"

    def get_uncompressed_json_as_df(self) -> Path:
        return self._get_uncompressed_file(self.get_json_as_dataframe(), "gz")

    def get_busted_airrtable(self) -> Path:
        return self.airr_table_inputs / "busted.csv.gz"

    def get_scfv_fasta(self) -> Path:
        return self.fasta_inputs / "scfv.fasta"

    def get_long_scfv_fastq(self) -> Path:
        return SeqIO.parse(self._get_uncompressed_file(self.fastq_inputs / "long_scfv.fastq.gz", "gz"), "fastq")

    def get_scfv_sequences(self) -> SeqIO.SeqRecord:
        return SeqIO.parse(self.get_scfv_fasta(), "fasta")

    def get_pg9_heavy_fasta(self) -> Path:
        return self.fasta_inputs / "PG9_H.fasta"

    def get_pg9_heavy_sequence(self) -> SeqIO.SeqRecord:
        return SeqIO.read(self.get_pg9_heavy_fasta(), "fasta")

    def get_pg9_light_fasta(self) -> Path:
        return self.fasta_inputs / "PG9_L.fasta"

    def get_pg9_light_sequence(self) -> SeqIO.SeqRecord:
        return SeqIO.read(self.get_pg9_light_fasta(), "fasta")

    def get_pg9_heavy_multiple_fasta(self) -> Path:
        return self.fasta_inputs / "PG9_H_multiple.fasta"

    def get_pg9_light_multiple_fasta(self) -> Path:
        return self.fasta_inputs / "PG9_L_multiple.fasta"

    def get_pg9_heavy_fasta_compressed(self, compression) -> Path:
        return self._get_file_compressed(self.get_pg9_heavy_fasta(), compression)

    def get_pg9_light_fasta_compressed(self, compression) -> Path:
        return self._get_file_compressed(self.get_pg9_light_fasta(), compression)

    def get_fasta_files(self) -> List[Path]:
        return [
            self.get_scfv_fasta(),
            self.get_pg9_heavy_multiple_fasta(),
            self.get_pg9_light_fasta(),
            self.get_pg9_light_multiple_fasta(),
        ]

    def get_adaptable_pentalty_test_seq(self) -> str:
        test_sequence = """
        GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACC
        ATCACTTGCCAGGCGAGTCAGGACATTAGCAACTATTTAAATTGGTATCAGCAGAAACCA
        GGGAAAGCCCCTAAGCTCCTGATCTACGATGCATCCAATTTGGAAACAGGGGTCCCATCA
        AGGTTCAGTGGAAGTGGATCTGGGACAGATTTTACTTTCACCATCAGCAGCCTGCAGCCT
        GAAGATATTGCAACATATTACTGTCAACAGTATGATAATTTCGGCGGAGGGACCAAGGTG
        GACATCAAAC""".replace(
            "\n", ""
        )
        return test_sequence

    def get_adaptable_pentalty_test_seq_scfv(self) -> str:
        scfv = """GGCGGCCGAGCTCGACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACCGTGTCACCATCACTT
            GCCGGGGCAAGTCAAATCATTAGCAACTATCTAAACTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTA
            TGCTACATCCAATTTGCACAGTGGGGTCCCATCACGCTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCA
            GCAGTCTGCAACCTGAAGATTTTGCAACTTACTACTGTCAACAGAGTTACAGTTTCCCGCCCACTTTCGGCGGAGGTACC
            AAGGTGGAGATCAAACGTACGGTTGCCGCTCCTTCTGTATTCATATTTCCGCCCTCCGATGAACAGCTTAAATCGGGCAC
            TGCTTCGGTAGTCTGCCTTCTGAATAATTTCTATCCCCGCGAGGCCAAGGTGCAATGGAAAGTCGACAATGCACTGCAAA
            GTGGAAACTCGCAAGAAAGCGTCACCGAACAGGACAGTAAAGATTCCACCTATAGCCTGTCATCGACACTTACCCTGAGT
            AAGGCTGATTACGAAAAGCACAAGGTTTACGCTTGCGAAGTAACTCACCAGGGCCTCTCAAGCCCTGTTACAAAGTCATT
            TAACAGAGGGGAATGCTAATTCTAGATAATTATCAAGGAGACAGTCATAATGAAATACCTATTGCCTACGGCAGCCGCTG
            GATTGTTATTACTCGCTGCCCAACCAGCCATGGCCGAGGTGCAGCTGCTCGAGGTGCAGCTGTTGGAGTCTGGGGGAGGC
            TTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGACTGACGATCTCATCCTACGCAATGTCATGGGT
            CCGCCAGGCTCCAGGGAANGGGGCTGGAGT"""
        return scfv.replace("\n", "").replace(" ", "")

    def get_catnap_heavy_nt(self) -> Path:
        return self.fasta_inputs / "catnap_nt_heavy.fasta"

    def get_catnap_light_nt(self) -> Path:
        return self.fasta_inputs / "catnap_nt_light.fasta"

    def get_monkey_edge_seq(self) -> str:
        return "TCCAGTCCCTGCAGGCCGGGAGGCAGGTGACCTCTGCCTCAGACCCCCACTCCAGACACCAGACAGAGGGGCAGGCCCCCCAGAACCAAAGTGGAGGGACGACCCGTCAAGGACAAACCAGACCAAGGGACACTGAGCCCAGCACGGGAAGGTCCCCAGATAGACCAGGAGGTTTCTGGAGGTGTCTGTGCCACAGTGGGGTATAGCAGCAGATCCGACTACGGTAGCAACTTTTGGGACTACTGGGGCCAGGGAGTCCTGGTCACCGTCTCCTCAGCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTCCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGC"

    def _get_file_compressed(self, uncompressed_file, compression) -> Path:
        tmp_file_name = uncompressed_file.name + f".{compression}"
        return_path = self.tmp_path / tmp_file_name
        with open(uncompressed_file, "rb") as f_in:
            if compression == "bz2":
                fn = bz2.open(return_path, "wb")
            elif compression == "gz":
                fn = gzip.open(return_path, "wb")
            with fn as f_out:
                shutil.copyfileobj(f_in, f_out)
        return return_path

    def _get_uncompressed_file(self, compressed_file, compression) -> Path:
        tmp_file_name = compressed_file.stem
        return_path = self.tmp_path / tmp_file_name
        if compression == "bz2":
            fn = bz2.open(compressed_file, "rb")
        elif compression == "gz":
            fn = gzip.open(compressed_file, "rb")
        with fn as f_in:
            with open(return_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return return_path


@pytest.fixture(scope="session", autouse=True)
def fixture_setup(tmp_path_factory):
    return SadieFixture(tmp_path_factory)


@pytest.fixture(scope="session", autouse=False)
def heavy_catnap_airrtable(fixture_setup) -> AirrTable:
    """ A permanant fixture of the catnap heavy airr table run through adaptable airrtable"""
    airr_api = Airr("human", adaptable=True)
    airrtable_heavy = airr_api.run_fasta(fixture_setup.get_catnap_heavy_nt())
    airrtable_heavy["cellid"] = airrtable_heavy["sequence_id"].str.split("_").str.get(0)
    return airrtable_heavy


@pytest.fixture(scope="session", autouse=False)
def light_catnap_airrtable(fixture_setup) -> AirrTable:
    """ A permanant fixture of the catnap light airr table run through adaptable airrtable"""
    airr_api = Airr("human", adaptable=True)
    airrtable_light = airr_api.run_fasta(fixture_setup.get_catnap_light_nt())
    airrtable_light["cellid"] = airrtable_light["sequence_id"].str.split("_").str.get(0)
    return airrtable_light
