"""Pytest conftest with all the fixture classes"""
import bz2
import gzip
import json
import shutil
from pathlib import Path
from typing import List

import pytest
from Bio import SeqIO
from sadie.airr import Airr
from sadie.airr.airrtable import AirrTable


def _get_file_compressed(tmp_path: Path, uncompressed_file: Path, compression: str) -> Path:
    """helper funcitons for compressing files"""
    tmp_file_name = uncompressed_file.name + f".{compression}"
    return_path = tmp_path / tmp_file_name
    with open(uncompressed_file, "rb") as f_in:
        if compression == "bz2":
            fn = bz2.open(return_path, "wb")
        elif compression == "gz":
            fn = gzip.open(return_path, "wb")
        with fn as f_out:
            shutil.copyfileobj(f_in, f_out)
    return return_path


def _get_uncompressed_file(tmp_path: Path, compressed_file: Path, compression: str) -> Path:
    """helper funcitons for uncompressing files"""
    tmp_file_name = compressed_file.stem
    return_path = tmp_path / tmp_file_name
    if compression == "bz2":
        fn = bz2.open(compressed_file, "rb")
    elif compression == "gz":
        fn = gzip.open(compressed_file, "rb")
    with fn as f_in:
        with open(return_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return return_path


class AirrSequences:
    """
    Organization of sequence fixtures that are used by Airr and Anarci primarily.
    Will be sub-classed by SadieFixtures
    """

    def __init__(self, tmp_path: Path, base_datadir: Path):
        self.base_datadir = base_datadir
        self.tmp_path = tmp_path
        self.fasta_inputs = self.base_datadir / "fasta_inputs/"
        self.fastq_inputs = self.base_datadir / "fastq_inputs/"
        self.airr_table_inputs = self.base_datadir / "airr_tables/"
        self.single_seqs_json = json.load(open(self.base_datadir / "single-sequences.json"))

    def get_pg9_heavy_fasta(self) -> Path:
        """get the file path for PG9 heavy chain fasta"""
        return self.fasta_inputs / "PG9_H.fasta"

    def get_pg9_heavy_sequence(self) -> SeqIO.SeqRecord:
        """get the sequence recored for the PG9 heavy chain"""
        return SeqIO.read(self.get_pg9_heavy_fasta(), "fasta")

    def get_pg9_light_fasta(self) -> Path:
        """get the file path for PG9 light chain fasta"""
        return self.fasta_inputs / "PG9_L.fasta"

    def get_pg9_light_sequence(self) -> SeqIO.SeqRecord:
        """get the sequence recored for the PG9 light chain"""
        return SeqIO.read(self.get_pg9_light_fasta(), "fasta")

    def get_pg9_heavy_multiple_fasta(self) -> Path:
        """get multiple fasta file path for multiple PG9 heavy chains"""
        return self.fasta_inputs / "PG9_H_multiple.fasta"

    def get_pg9_heavy_multiple_fasta_seqs(self) -> List[SeqIO.SeqRecord]:
        """get multiple SeqRecords for multiple PG9 heavy chains"""
        return list(SeqIO.parse(self.get_pg9_heavy_multiple_fasta(), "fasta"))

    def get_pg9_light_multiple_fasta(self) -> Path:
        """get multiple fasta file path for multiple PG9 light chains"""
        return self.fasta_inputs / "PG9_L_multiple.fasta"

    def get_pg9_light_multiple_fasta_seqs(self) -> List[SeqIO.SeqRecord]:
        """get multiple SeqRecords for multiple PG9 light chains"""
        return list(SeqIO.parse(self.get_pg9_light_multiple_fasta(), "fasta"))

    def get_pg9_heavy_fasta_compressed(self, compression: str) -> Path:
        """Get a pg9 heavy fasta file compressed with a given compression

        Parameters
        ----------
        compression : str
            the compression to use., eg. "gz" or "bz2"

        Returns
        -------
        Path
            file path for the compressed file
        """
        return _get_file_compressed(self.tmp_path, self.get_pg9_heavy_fasta(), compression)

    def get_pg9_light_fasta_compressed(self, compression) -> Path:
        """Get a pg9 light fasta file compressed with a given compression

        Parameters
        ----------
        compression : str
            the compression to use., eg. "gz" or "bz2"

        Returns
        -------
        Path
            file path for the compressed file
        """
        return _get_file_compressed(self.tmp_path, self.get_pg9_light_fasta(), compression)

    def get_scfv_fasta(self) -> Path:
        """get a fasta file path for the scfv file that contains duplicated scfv files"""
        return self.fasta_inputs / "scfv.fasta"

    def get_scfv_sequences(self) -> SeqIO.SeqRecord:
        """get the SeqRecords for the scfv fasta file that contains duplicated scfv files"""
        return SeqIO.parse(self.get_scfv_fasta(), "fasta")

    def get_long_scfv_fastq(self) -> List[SeqIO.SeqRecord]:
        """get a fastq file path for the scfv file records that contains many unique scfv files"""
        return list(
            SeqIO.parse(_get_uncompressed_file(self.tmp_path, self.fastq_inputs / "long_scfv.fastq.gz", "gz"), "fastq")
        )

    def get_catnap_heavy_nt(self) -> Path:
        """get a file path for the catnap heavy nt file"""
        return self.fasta_inputs / "catnap_nt_heavy.fasta"

    def get_catnap_light_nt(self) -> Path:
        """get a file path for the catnap light nt file"""
        return self.fasta_inputs / "catnap_nt_light.fasta"

    def get_catnap_heavy_aa(self) -> Path:
        """get a file path for the catnap heavy aa file"""
        return self.fasta_inputs / "catnap_aa_heavy_sample.fasta"

    def get_catnap_light_aa(self) -> Path:
        """get a file path for the catnap light aa file"""
        return self.fasta_inputs / "catnap_aa_light_sample.fasta"

    def get_oas_fasta(self) -> Path:
        """Get a random 1000 subsample fasta file path of the OAS data set."""
        return self.fasta_inputs / "OAS_subsample_1000.fasta"

    def get_dog_aa_seqs(self) -> Path:
        "A fasta file path containing random canine AA sequences"
        return self.fasta_inputs / "random_dog_contigs_aa.fasta"

    def get_fasta_files(self) -> List[Path]:
        """Get a list of different fasta files for CLI testing. scfv, pg9_h, pg9_h_multi, pg9_l, pg9_l_multi"""
        return [
            self.get_scfv_fasta(),
            self.get_pg9_heavy_multiple_fasta(),
            self.get_pg9_light_fasta(),
            self.get_pg9_light_multiple_fasta(),
        ]

    def get_adaptable_pentalty_test_seq(self) -> str:
        """get a signle sequence for the testing sequence for adaptable pentalties"""
        return self.single_seqs_json["adaptible_pentalty_test_seq"]

    def get_adaptable_pentalty_test_seq_scfv(self) -> str:
        """get a single scfv sequence for the testing sequence for adaptable pentalties"""
        return self.single_seqs_json["adaptible_pentalty_test_seq_scfv"]

    def get_monkey_edge_seq(self) -> str:
        """get a single sequence for the testing sequence for testing the weird macaque edge case"""
        return self.single_seqs_json["monkey_edge_case"]


class AirrTables:
    """A class for organization of airrtable related fixtures"""

    def __init__(self, tmp_path: Path, base_datadir: Path):
        self.base_datadir = base_datadir
        self.tmp_path = tmp_path
        self.airr_table_inputs = self.base_datadir / "airr_tables/"

    def get_dog_airrtable(self) -> Path:
        """get a file path for the dog airr table csv.gzjj"""
        return self.airr_table_inputs / "dog_igh.csv.gz"

    def get_linked_airrtable(self) -> Path:
        """get file path to a linked airrtable for spliting and reassembling"""
        return self.airr_table_inputs / "linked_airrtable_scfv.csv.gz"

    def get_json_as_dataframe(self) -> Path:
        """get a file path for the airr table that is in json.gz"""
        return self.airr_table_inputs / "airrtable.json.gz"

    def get_uncompressed_json_as_dataframe(self) -> Path:
        """get a file path for the airr table that is in json."""
        return _get_uncompressed_file(self.tmp_path, self.get_json_as_dataframe(), "gz")

    def get_busted_airrtable(self) -> Path:
        """get a file path for the airr table that is in csv. contains missing airrtables"""
        return self.airr_table_inputs / "busted.csv.gz"

    def get_catnap_heavy_airrtable(self) -> Path:
        """get a file path for the CATNAP heavy airr table results. should match fasta input"""
        return self.airr_table_inputs / "catnap_heavy_airrtable.feather"

    def get_catnap_light_airrtable(self) -> Path:
        """get a file path for the CATNAP light airr table results. Should match fasta input"""
        return self.airr_table_inputs / "catnap_light_airrtable.feather"

    def get_imgt_airrtable(self) -> Path:
        """get a file path for the OAS sequences ran through IMGT Hi-V"""
        return self.airr_table_inputs / "imgt_v_quest_airr.tsv.gz"


class SadieFixture(AirrSequences, AirrTables):
    def __init__(self, tmp_path_factory):
        self.tmp_path = tmp_path_factory.mktemp("sadie_fixture")
        self.base_datadir = Path("tests/data/fixtures/")
        self.reference_data = self.base_datadir / "reference/"
        super().__init__(self.tmp_path, self.base_datadir)

    def get_reference_dataset_csv(self) -> Path:
        """
        get a file path for the reference dataset csv
        this is a csv that contains all the reference data in a nice csv dataframe
        """
        return self.reference_data / "reference_object_dataframe.csv.gz"

    def get_card(self) -> Path:
        """get a png file path that has a card image. This is a nonsense file path to test unexpected files"""
        return self.base_datadir / "card.png"


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
