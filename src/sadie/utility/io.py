import glob
from pathlib import Path
from pprint import pformat
import warnings
from filetype import guess
from filetype.types.base import Type
from filetype.types.archive import Gz, Bz2
import gzip
import bz2
from typing import IO, Dict, Iterator, List, TextIO, Union, Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.AbiIO import AbiIterator
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO.QualityIO import FastqPhredIterator

# exceptions handled in own file
from sadie.utility.exception import NotAValidSequenceFile, NotAValidCompression, DirectoryExistsError


def guess_input_compression(input_path: Union[str, Path]) -> Union[str, None]:
    """Given a path, guess if it's compressed gz,bz2 or directory.
    If not compressed returned None

    Parameters
    ----------
    input_path : str
        input is compressed

    Returns
    -------
    Union[str,None]
        returns gz, bz2 or None

    Raises
    ------
    NotImplementedError
        If the file type compression is unsupported
    """

    # check if directory first
    if Path(input_path).is_dir():
        return "directory"

    # if not directory, check if it's bz2 or gz2 archive
    _filetype: Union[Type, None] = guess(str(input_path))
    if not _filetype:
        return None
    extension: str = _filetype.extension
    if isinstance(_filetype, (Gz, Bz2)):
        return extension
    else:
        raise NotImplementedError(f"filetype {extension} not supproted")


def get_file_buffer(file: Path, compression: Union[str, None] = None, mode: str = "rt") -> Union[TextIO, IO[Any]]:
    """Open a file with the correct file buffer

    Parameters
    ----------
    file : Path
        file path
    compression : str
        compression type, 'bz2', 'gz', None

    Returns
    -------
    TextIO
        Usually a TextIOWrapper

    Raises
    ------
    TypeError
        Can't determine file type
    """
    # get file buffer for compression
    if compression is None:
        return open(file, mode)
    elif compression == "gz":
        return gzip.open(file, "rt")
    elif compression == "bz2":
        return bz2.open(file, "rt")
    else:
        raise TypeError(f"{file} can't determine file encoding")


def get_sequence_file_type(file: Union[str, Path]) -> str:
    """Get file type of file

    Parameters
    ----------
    file : Path
        A file type

    Returns
    -------
    str
        the sequence file type, either 'abi','fasta' or 'fastq'
    Raises
    ------
    TypeError
        Can't determine filetype
    """
    file = Path(file)
    compression_extension = guess_input_compression(file)

    # The parse will consume the buffer so have to open it every time
    try:
        file_buffer = get_file_buffer(file, compression_extension)
        SeqIO.parse(file_buffer, "fasta").__next__()
        return "fasta"
    except (StopIteration, ValueError):
        pass

    try:
        file_buffer = get_file_buffer(file, compression_extension)
        SeqIO.parse(file_buffer, "fastq").__next__()
        return "fastq"
    except (StopIteration, ValueError):
        pass
    try:
        file_buffer = get_file_buffer(file, compression_extension, "rb")
        SeqIO.parse(file_buffer, "abi").__next__()
        return "abi"
    except (StopIteration, ValueError, OSError):
        raise NotAValidSequenceFile(f"can't determine sequence file type of {file}")


def get_sequence_file_iter(
    file: Union[str, Path, TextIO, IO[Any]], file_type: str
) -> Union[AbiIterator, FastaIterator, FastqPhredIterator]:
    """Get a sequence file iterator from SeqIO module"""
    if file_type not in ["fasta", "fastq", "abi", "abi-trim"]:
        raise NotImplementedError(f"{file_type} is not a supported sequence file type")
    return SeqIO.parse(file, file_type)


class SadieInputFile:
    """Sadie Input File will handle input files. Use SadieInputDir for input directories"""

    def __init__(
        self,
        input_path: Union[Path, str],
        input_format: str = "infer",
        compression_format: Union[str, None] = "infer",
    ):
        self.input = input_path
        self.compression_format_inferred = False
        self.file_type_inferred = False

        # infer compression no matter what
        inferred_compression_format = guess_input_compression(self.input)

        # check if directory
        if inferred_compression_format == "directory":
            raise TypeError(f"{self.input} is a directory, use SadieInputDir instead")

        # if we are to infer, we will get a filetype back
        if compression_format == "infer":
            self.compression_format_inferred = True
            self.compression_format = inferred_compression_format
        else:
            # set what the user put in
            if compression_format not in ["gz", "bz2", None]:
                raise NotAValidCompression(f"{compression_format} is not a valid compression type, need gz or bz2")

            # only warn if we got it wrong
            if compression_format != inferred_compression_format:
                warnings.warn(
                    f"{self.input} is detected to be {inferred_compression_format} but you specified {compression_format}",
                    UserWarning,
                )
            self.compression_format = compression_format

        # handle the input format next
        self.input_format = input_format

        # if inferred, try to guess type
        if self.input_format == "infer":
            self.input_format_inferred = True
            self.input_format = get_sequence_file_type(self.input)
        else:
            # else its explicit, but check that its implemented
            if self.input_format not in ["fasta", "fastq", "abi", "abi-trim"]:
                raise NotAValidSequenceFile(
                    f"{self.input_format} is not a supported sequence file type, only fasta, fastq, abi, abi-trim"
                )

        # get open file
        self.input = Path(self.input)
        if self.input_format in ["abi", "abi-trim"]:
            # if abi, force rb mode, no matter what compression format is
            mode = "rb"
        else:
            mode = "rt"
        self.open_input = get_file_buffer(self.input, self.compression_format, mode)

        # finally get the generator
        self.sequence_generator = get_sequence_file_iter(self.open_input, self.input_format)

    def get_seq_records(self) -> List[SeqRecord]:
        """Return all seqeunces as a list of sequence records

        Returns
        -------
        List[SeqRecord]
            A list of SeqRecord objects from file input
        """
        return [i for i in self]

    def __iter__(self) -> Iterator[SeqRecord]:
        for seq in self.sequence_generator:
            yield seq

    def __repr__(self) -> str:
        property_dict: Dict[str, Union[str, bool, None]] = {
            "input_path": str(self.input),
            "input_file_type": self.input_format,
            "input_inferred": self.file_type_inferred,
            "input_compression": self.compression_format,
            "compression_inferred": self.compression_format_inferred,
        }
        return "SadieInput" + pformat(property_dict, indent=4)

    def __str__(self) -> str:
        return str(self.input)


class SadieInputDir:
    """Sadie Input Dir will handle input directories. Use SadieInputFile for input files"""

    def __init__(
        self,
        directory_path: Union[Path, str],
        direcotry_file_format: str = "infer",
        compression_format: Union[str, None] = "infer",
        recurse: bool = False,
        ignore_bad_seq_files: bool = True,
    ):
        if Path(directory_path).is_file():
            raise TypeError(f"{directory_path} is a file, use SadieInputFile instead")
        self.directory_input = directory_path
        self.directory_file_format = direcotry_file_format
        self.directory_file_format_inferred = False
        self.compression_file_format = compression_format
        self.recurse = recurse
        self.ignore_bad_seq_files = ignore_bad_seq_files

        if self.directory_file_format == "infer":
            self.directory_file_format_inferred = True
        else:
            # else its explicit, but check that its implemented
            if self.directory_file_format not in ["fasta", "fastq", "abi", "abi-trim"]:
                raise NotImplementedError(
                    f"{self.directory_input} is not a supported sequence file type, only fasta, fastq, abi, abi-trim"
                )

        self.sequence_files: List[SadieInputFile] = self._get_parse()
        self.sequence_files_dict: Dict[Union[Path, str], str] = {i.input: i.input_format for i in self.sequence_files}

    def get_combined_seq_records(self) -> List[SeqRecord]:
        """For all sequence files in list, combine them into a single list of sequence records

        Returns
        -------
        List[SeqRecord]
            Combined list of all sequence records from all files in directory
        """
        _records = []
        for x in self.sequence_files:
            _records += x.get_seq_records()
        return _records

    def _get_parse(self) -> List[SadieInputFile]:
        """Private method to handle parsing directory"""
        _files: List[SadieInputFile] = []
        for x in glob.glob(str(self.directory_input) + "/*", recursive=self.recurse):
            try:
                _files.append(SadieInputFile(x, self.directory_file_format, self.compression_file_format))
            except NotAValidSequenceFile as e:
                if self.ignore_bad_seq_files:
                    warnings.warn(f"{x} is not a valid sequence file, skipping", UserWarning)
                else:
                    raise e
        return _files

    def __repr__(self) -> str:
        return "\nSadieDir:\n" + pformat(self.sequence_files_dict, indent=4)


class SadieOutput:
    """Sadie Output will handle output files, infer compression and filetype"""

    def __init__(
        self,
        output_path: Union[Path, str],
        overwrite: bool = True,
    ):
        """Make an output object from filepath

        Parameters
        ----------
        output_path : Union[Path, str]
            given a filepath
        overwrite : bool, optional
            is it okay to overwrite

        Raises
        ------
        DirectoryExistsError
            if output_path is a directory
        FileExistsError
            if output_path exists and overwrite is False
        ValueError
            if output_path is not a valid filepath type or compression isnot a valide filepath type
        """
        if Path(output_path).is_dir() and Path(output_path).exists():
            raise DirectoryExistsError(f"{output_path}  exists and is directory directory instead")
        elif Path(output_path).exists():
            if overwrite:
                warnings.warn(f"{output_path} exists, overwriting", UserWarning)
            else:
                raise FileExistsError(f"{output_path} is exists and is directory directory instead")

        # set accepted formats
        self._accepted_output_format = ["json", "csv", "tsv", "feather", "stdout"]
        self._accepted_output_compression = ["gz", "bz2"]
        self.output_path = Path(output_path)
        self.suffixes = list(map(lambda x: x.lstrip("."), self.output_path.suffixes))

        if not self.suffixes:
            raise ValueError(f"{self.output_path} has no suffix")
        else:
            # base path is path without suffix
            self.base_path = Path(self.output_path.__str__().split("." + self.suffixes[0])[0])

        if len(self.suffixes) <= 1:
            self.compression_format = None
        elif len(self.suffixes) == 2:
            self.compression_format = self.suffixes[-1]
        else:
            raise ValueError(f"{self.output_path} has too many suffixes; {self.suffixes}")

        # collect output format
        self.output_format = self.suffixes[0]
        if self.output_format not in self._accepted_output_format:
            raise ValueError(f"{self.output_format} is not a supported output format, {self._accepted_output_format}")
        if self.compression_format not in self._accepted_output_compression and self.compression_format is not None:
            raise ValueError(f"{self.compression_format} is not a supported output format, {self.compression_format}")

        self.overwrite = overwrite

        # make property dict that can be sliced later
        self.property_dict = {
            "output_path": str(self.output_path),
            "base_path": str(self.base_path),
            "output_format": self.output_format,
            "compression_format": str(self.compression_format),
            "overwrite": self.overwrite,
        }

    def __repr__(self) -> str:
        return pformat(self.property_dict, indent=4)

    def __str__(self) -> str:
        return str(self.output_path)
