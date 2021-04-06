import os
from pathlib import Path

import pytest
from Bio.SeqIO import SeqRecord
from sadie.utility import SadieIO
from sadie.utility.exception import IOInferError


def get_file(file):
    """Helper method for test execution."""
    _file = os.path.join(os.path.abspath(os.path.dirname(__file__)), f"fixtures/{file}")
    if not os.path.exists(_file):
        raise FileNotFoundError(_file)
    return _file


def test_io_static_methods():
    assert SadieIO.guess_input_compression(get_file("multiple.fastq.gz")) == "gz"
    assert SadieIO.guess_input_compression(get_file("multiple.fastq")) is None
    assert SadieIO.guess_input_compression(get_file("multiple.fastq.bz2")) == "bz2"
    assert SadieIO.guess_input_compression(get_file("fastq_folder")) == "directory"

    assert SadieIO.guess_sequence_file_type(get_file("multiple.fasta")) == "fasta"
    assert SadieIO.guess_sequence_file_type(get_file("multiple.fastq.gz")) == "fastq"
    assert SadieIO.guess_sequence_file_type(get_file("ab1_files/file1.ab1")) == "abi"

    _dict = SadieIO.get_file_type_dict(get_file("fastq_folder"))
    assert isinstance(_dict, dict)
    assert all([x == "fastq" for x in _dict.values()])

    with pytest.raises(TypeError):
        SadieIO.get_file_type_dict(get_file("multiple.fasta"))


def test_io_single_files():
    # infer filetypes of non compressed
    input_file = get_file("multiple.fasta")
    io = SadieIO(input_file, out_format="csv")
    assert io.input == Path(input_file)
    assert io.infer_input
    assert not io.input_compressed
    assert io.input_file_type == "fasta"
    assert not io.isdir

    input_file = get_file("multiple.fastq")
    io = SadieIO(input_file, out_format="csv")
    assert io.input == Path(input_file)
    assert io.infer_input
    assert not io.input_compressed
    assert io.input_file_type == "fastq"
    assert not io.isdir

    input_file = get_file("ab1_files/file1.ab1")
    io = SadieIO(input_file, out_format="csv")
    assert io.input == Path(input_file)
    assert io.infer_input
    assert not io.input_compressed
    assert io.input_file_type == "abi"
    assert not io.isdir

    # infer filetypes of non compressed
    input_file = get_file("multiple.fasta.gz")
    io = SadieIO(input_file, out_format="csv")
    assert io.input == Path(input_file)
    assert io.infer_input
    assert io.input_compressed == "gz"
    assert io.input_file_type == "fasta"
    assert not io.isdir

    input_file = get_file("multiple.fastq.gz")
    io = SadieIO(input_file, out_format="csv")
    assert io.input == Path(input_file)
    assert io.infer_input
    assert io.input_compressed == "gz"
    assert io.input_file_type == "fastq"
    assert not io.isdir

    input_file = get_file("ab1_files/file1.ab1.gz")
    io = SadieIO(input_file, out_format="csv")
    assert io.input == Path(input_file)
    assert io.infer_input
    assert io.input_compressed == "gz"
    assert io.input_file_type == "abi"
    assert not io.isdir


def test_io_folders():
    # infer filetypes of non compressed
    input_file = get_file("fasta_folder")
    io = SadieIO(input_file, out_format="csv")
    assert io.input == Path(input_file)
    assert io.infer_input
    assert io.input_compressed == "directory"
    assert io.input_file_type == "fasta"
    assert io.isdir
    assert all([isinstance(i, SeqRecord) for i in io.get_input_records()])

    input_file = get_file("fastq_folder")
    io = SadieIO(input_file, out_format="csv")
    assert io.input == Path(input_file)
    assert io.infer_input
    assert io.input_compressed == "directory"
    assert io.input_file_type == "fastq"
    assert io.isdir
    assert all([isinstance(i, SeqRecord) for i in io.get_input_records()])

    input_file = get_file("ab1_files")
    io = SadieIO(input_file, out_format="csv")
    assert io.input == Path(input_file)
    assert io.infer_input
    assert io.input_compressed == "directory"
    assert io.input_file_type == "abi"
    assert io.isdir
    assert all([isinstance(i, SeqRecord) for i in io.get_input_records()])


def test_output(tmpdir):
    # uses tmpdir pytest fixture
    input_file = get_file("ab1_files/file1.ab1.gz")

    # case 1 no output path is set, but infer is set by default
    with pytest.raises(IOInferError):
        io_output = SadieIO(input_file)

    with pytest.raises(ValueError):
        io_output = SadieIO(input_file, output_path="myfile.gsv")

    with pytest.warns(UserWarning):
        io_output = SadieIO(input_file, output_path="myfile.tsv.zip")

    # case 2 no output path is set, but infer is set by default
    for suffix in ["csv", "tsv", "feather", "json"]:
        io_output = SadieIO(input_file, out_format=f"{suffix}")
        assert io_output.output == Path(f"file1.{suffix}")
        assert io_output.output_compressed == "infer"
        assert io_output.output_format == suffix
        assert io_output.infered_output_format is None
        assert io_output.infered_output_compression is None

    # case 3 no output path is set, but infer is set by default and compression
    for suffix in ["csv", "tsv", "feather", "json"]:
        for compression in ["gz", "bz2"]:
            io_output = SadieIO(input_file, out_format=f"{suffix}", compressed=f"{compression}")
            assert io_output.output == Path(f"file1.{suffix}.{compression}")
            assert io_output.output_compressed == compression
            assert io_output.infered_output_compression is None
            assert io_output.output_format == suffix
            assert io_output.infered_output_format is None

    # case 4 output and input are set with infer type and compression
    for suffix in ["csv", "tsv", "feather", "json"]:
        for compression in ["gz", "bz2"]:
            output_file_dummy = Path(f"somedummy_file.{suffix}.{compression}")
            io_output = SadieIO(input_path=input_file, output_path=output_file_dummy)
            assert io_output.output == output_file_dummy
            assert io_output.output_compressed == "infer"
            assert io_output.infered_output_compression == compression
            assert io_output.output_format == "infer"
            assert io_output.infered_output_format == suffix
