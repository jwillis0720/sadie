# All of these ignores are because pylance can't figure out fixtures
# pyright: reportUnknownMemberType=false
# pyright: reportUnknownParameterType=false
# pyright: reportUnknownVariableType=false
# pyright: reportMissingParameterType=false

# from gzip import GzipFile
# from pathlib import Path

from io import TextIOWrapper
from pathlib import Path
import shutil
from typing import List

import pytest

# third part
from Bio.SeqIO.AbiIO import AbiIterator
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio.SeqRecord import SeqRecord

# package
from sadie.utility.io import (
    SadieInputDir,
    SadieInputFile,
    SadieOutput,
    get_file_buffer,
    get_sequence_file_iter,
    get_sequence_file_type,
    guess_input_compression,
)
from sadie.utility.exception import (
    DirectoryExistsError,
    NotAValidCompression,
    NotAValidSequenceFile,
)


def test_io_methods(fixture_setup):
    """Text methods in IO that are used mostly in objects"""

    # get fastq files
    multi_fastq: Path = fixture_setup.get_multiple_fastq()
    multi_fastq_gz: Path = fixture_setup.get_multiple_fastq_compressed("gz")
    multi_fastq_bz: Path = fixture_setup.get_multiple_fastq_compressed("bz2")
    fastq_dir_path: Path = fixture_setup.fastq_inputs  # this is a directory
    # this is not implemented
    multi_fastq_xz: Path = fixture_setup.get_multiple_fastq_compressed("xz")

    # Step 1 is to explicitly get buffer
    # "rt" mode
    multi_fastq_fb = get_file_buffer(multi_fastq)
    multi_fastq_gz_fb = get_file_buffer(multi_fastq_gz, "gz")
    multi_fastq_bz_fb = get_file_buffer(multi_fastq_bz, "bz2")
    assert all(isinstance(x, TextIOWrapper) for x in [multi_fastq_fb, multi_fastq_gz_fb, multi_fastq_bz_fb])

    with pytest.raises(TypeError):
        get_file_buffer(multi_fastq_bz, "zip")

    # Step 2 is to guess compression type
    assert not guess_input_compression(multi_fastq)
    assert guess_input_compression(multi_fastq_gz) == "gz"
    assert guess_input_compression(multi_fastq_bz) == "bz2"
    assert guess_input_compression(fastq_dir_path) == "directory"
    with pytest.raises(NotImplementedError):
        guess_input_compression(multi_fastq_xz)

    # need to grab fasta files here
    fasta_file = fixture_setup.get_pg9_heavy_fasta()
    fasta_file_gz = fixture_setup.get_pg9_heavy_fasta_compressed("gz")
    fasta_file_bz2 = fixture_setup.get_pg9_heavy_fasta_compressed("bz2")
    assert all(get_sequence_file_type(x) == "fasta" for x in [fasta_file, fasta_file_gz, fasta_file_bz2])
    assert all(get_sequence_file_type(x) == "fastq" for x in [multi_fastq, multi_fastq_bz, multi_fastq_gz])
    assert all(get_sequence_file_type(x) == "abi" for x in fixture_setup.get_abi_files())

    with pytest.raises(NotAValidSequenceFile):
        # don't accept phylip format
        get_sequence_file_type(fixture_setup.get_phy_file())

    # get explicit iterators
    assert isinstance(get_sequence_file_iter(fixture_setup.get_abi_files()[0], "abi"), AbiIterator)
    assert isinstance(get_sequence_file_iter(fixture_setup.get_abi_files()[0], "abi-trim"), AbiIterator)
    assert isinstance(get_sequence_file_iter(fasta_file, "fasta"), FastaIterator)
    assert isinstance(get_sequence_file_iter(multi_fastq, "fastq"), FastqPhredIterator)

    with pytest.raises(NotImplementedError):
        get_sequence_file_iter(fixture_setup.get_phy_file(), "phy")


def test_sadie_input(fixture_setup):
    multi_fastq: Path = fixture_setup.get_multiple_fastq()
    multi_fastq_gz: Path = fixture_setup.get_multiple_fastq_compressed("gz")
    multi_fastq_bz: Path = fixture_setup.get_multiple_fastq_compressed("bz2")
    abi_files: List[Path] = fixture_setup.get_abi_files()
    fastq_dir_path: Path = fixture_setup.fastq_inputs  # this is a directory

    # infer both compression and filetype
    SadieInputFile(multi_fastq)
    SadieInputFile(multi_fastq_gz)
    SadieInputFile(multi_fastq_bz)

    # check abi files
    for file in abi_files:
        SadieInputFile(file)
        SadieInputFile(file, "abi")
        SadieInputFile(file, "abi-trim")
    SadieInputFile(multi_fastq, "fastq", None)
    SadieInputFile(multi_fastq_gz, "fastq", "gz")
    SadieInputFile(multi_fastq_bz, "fastq", "bz2")
    # test __str__
    print(SadieInputFile(multi_fastq_bz, "fastq", "bz2"))
    print(SadieInputFile(multi_fastq_bz, "fastq", "bz2").__repr__())

    gzip_fastq = SadieInputFile(multi_fastq_gz, "fastq", "gz")
    bzip_fastq = SadieInputFile(multi_fastq_bz, "fastq", "bz2")

    # test equivalence of get_seq_records across compression
    assert [i.seq for i in gzip_fastq.get_seq_records()] == [i.seq for i in bzip_fastq.get_seq_records()]
    assert [i.id for i in gzip_fastq.get_seq_records()] == [i.id for i in bzip_fastq.get_seq_records()]

    # test __iter__ directly
    assert [i.seq for i in gzip_fastq] == [i.seq for i in bzip_fastq]

    with pytest.raises(TypeError):
        # no paths in input file
        SadieInputFile(fastq_dir_path)
    with pytest.raises(NotAValidCompression):
        # bad compression
        SadieInputFile(multi_fastq_bz, "fastq", "bz")

    with pytest.warns(UserWarning):
        # compression does not match declared
        SadieInputFile(multi_fastq_gz, "fastq", "bz2")
    phy_file = fixture_setup.get_phy_file()
    with pytest.raises(NotAValidSequenceFile):
        SadieInputFile(phy_file, "phy")

    for seq in SadieInputFile(multi_fastq_gz):
        assert isinstance(seq, SeqRecord)


def test_sadie_directory(fixture_setup, tmpdir_factory):

    with pytest.raises(TypeError):
        # raise on a single file
        SadieInputDir(fixture_setup.get_abi_files()[0], "abi")

    abi_inputs_dir: Path = fixture_setup.abi_inputs
    fasta_inputs_dir: Path = fixture_setup.fasta_inputs
    fastq_inputs_dir: Path = fixture_setup.fastq_inputs

    abi_inputs_dir_obj = SadieInputDir(abi_inputs_dir)
    fasta_inputs_dir_obj = SadieInputDir(fasta_inputs_dir)
    fastq_inputs_dir_obj = SadieInputDir(fastq_inputs_dir)

    assert all(
        [x.directory_file_format_inferred for x in [abi_inputs_dir_obj, fasta_inputs_dir_obj, fastq_inputs_dir_obj]]
    )

    abi_inputs_dir_obj = SadieInputDir(abi_inputs_dir, "abi")
    fasta_inputs_dir_obj = SadieInputDir(fasta_inputs_dir, "fasta")
    fastq_inputs_dir_obj = SadieInputDir(fastq_inputs_dir, "fastq")
    with pytest.raises(NotImplementedError):
        SadieInputDir(fastq_inputs_dir, "fastz")

    with pytest.raises(ValueError):
        # if give fasta and specified fastq, will complain when try to get seq records
        SadieInputDir(fasta_inputs_dir, "fastq").get_combined_seq_records()
    fasta_inputs_dir_obj.__repr__()
    fastq_inputs_dir_obj.__repr__()
    abi_inputs_dir_obj.__repr__()
    mix_dir_path = tmpdir_factory.mktemp("test_sadie_directory")
    for directory in [abi_inputs_dir, fasta_inputs_dir, fastq_inputs_dir]:
        for file in directory.glob("*"):
            shutil.copy(file, mix_dir_path + "/.")
    mixed_dir_object = SadieInputDir(mix_dir_path)
    mixed_dir_object.__repr__()
    mixed_dir_object.get_combined_seq_records()

    # now add a bad file to the mix
    shutil.copy(fixture_setup.get_phy_file(), mix_dir_path + "/.")
    with pytest.warns(UserWarning):
        # warn that we have a bum file in the directory
        mixed_dir_object = SadieInputDir(mix_dir_path)

    with pytest.raises(NotAValidSequenceFile):
        mixed_dir_object = SadieInputDir(mix_dir_path, ignore_bad_seq_files=False)


def test_sadie_output(fixture_setup, tmpdir_factory):
    output_dir = tmpdir_factory.mktemp("test_sadie_output")
    abi_output_dir = fixture_setup.abi_inputs
    fasta_file = fixture_setup.get_pg9_heavy_fasta()
    shutil.copytree(abi_output_dir, output_dir + "/existing_dir")
    shutil.copy(fasta_file, output_dir + "/fasta_file.feather")

    with pytest.raises(DirectoryExistsError):
        SadieOutput(output_dir.join("/existing_dir"))
    with pytest.raises(FileExistsError):
        SadieOutput(output_dir.join("/fasta_file.feather"), overwrite=False)
    with pytest.warns(UserWarning):
        SadieOutput(output_dir.join("/fasta_file.feather"))
    with pytest.raises(ValueError):
        SadieOutput(output_dir.join("/output.feather.xs"))
    with pytest.raises(ValueError):
        SadieOutput(output_dir.join("/output.blah.gz"))

    with pytest.raises(ValueError):
        SadieOutput(output_dir.join("/output"))
    with pytest.raises(ValueError):
        SadieOutput(output_dir.join("/output.feather.gz.gz"))

    output_object = SadieOutput(output_dir.join("/output.json"))
    assert output_object.base_path == output_dir.join("/output")
    assert not output_object.compression_format
    assert output_object.output_format == "json"

    output_object = SadieOutput(output_dir.join("/output.json.gz"))
    assert str(output_object.base_path) == str(output_dir.join("/output"))
    assert output_object.compression_format == "gz"
    assert output_object.output_format == "json"

    output_object = SadieOutput(output_dir.join("/output.feather.bz2"))
    assert str(output_object.base_path) == str(output_dir.join("/output"))
    assert output_object.compression_format == "bz2"
    assert output_object.output_format == "feather"

    output_object = SadieOutput(output_dir.join("/output.csv.bz2"))
    assert str(output_object.base_path) == str(output_dir.join("/output"))
    assert output_object.compression_format == "bz2"
    assert output_object.output_format == "csv"
    print(output_object.__repr__())
    print(output_object)
