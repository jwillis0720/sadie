from pathlib import Path
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio.SeqIO.AbiIO import AbiIterator
from io import TextIOWrapper
from filetype import guess
from typing import Union
from pprint import pformat
import glob
import gzip
import bz2
import itertools


class SadieIO:
    def __init__(self, input_path: Path, output_path=None, in_format="infer", out_format="infer", compressed=None):
        """Constructor for SadieIO to aid in input output operations

        Parameters
        ----------
        input_path : Path
            The input Path
        output_path : Tuple[Path,None], default=None
            The output Path, can be inferred from input if none
        in_format : str , default='infer'
            an explict input format, can infer, ex: 'fasta','fastq','abi'
        out_format : str, default='infer'
            an explicit ouptut format
        compressed : str, defualt=None
            the compression type of the output, "None","gzip","bzip2"
        """
        # input path
        self.input = input_path

        # we assume input is not a directorys
        self._isdir = False

        # check if input is a directory:
        if self.input.is_dir():
            self._isdir = True
            self.all_files_in_dir = self.get_file_type_dict(self.input)

        # check if input is compressed - gzip or bzip2
        self.input_compressed = self.guess_input_compression(self.input)

        # check if input format was specified
        self.infer_input = False
        if in_format == "infer":
            self.infer_input = True

        # infer input
        if self.infer_input:
            self.input_file_type = self.guess_sequence_file_type(self.input)
        # input type was set by user
        else:
            self.input_file_type = in_format

        # output
        self.output = output_path

    @property
    def input(self) -> Path:
        return self._input

    @input.setter
    def input(self, input_path: Union[str, Path]):
        if isinstance(input_path, str):
            input_path = Path(input_path)
        if not isinstance(input_path, Path):
            raise TypeError(f"{input_path} needs to be str or path")
        if not input_path.exists:
            raise FileNotFoundError(f"{input_path} not found")
        self._input = input_path

    @property
    def isdir(self):
        return self._isdir

    @property
    def input_compressed(self):
        return self._input_compressed

    @input_compressed.setter
    def input_compressed(self, compressed_format: str):
        if compressed_format not in ["bzip2", "gzip", "directory", None]:
            raise TypeError(f"{compressed_format} needs to be bzip2, gzip or None")
        self._input_compressed = compressed_format

    @property
    def input_file_type(self) -> str:
        return self._input_file_type

    @input_file_type.setter
    def input_file_type(self, input_file: str):
        if input_file not in ["fasta", "abi", "fastq"]:
            raise NameError(f"{input_file} needs to be fasta,abi,fastq")
        self._input_file_type = input_file

    def _get_parse(self, format, mode="rt"):
        return itertools.chain(
            *[
                SeqIO.parse(SadieIO.get_file_buffer(x, self.guess_input_compression(x), mode), format)
                for x in self.all_files_in_dir
            ]
        )

    def get_input_records(self) -> Union[FastaIterator, FastqPhredIterator, AbiIterator, itertools.chain]:
        if not self.isdir:
            if self.input_file_type == "fasta":
                return SeqIO.parse(SadieIO.get_file_buffer(self.input, self.input_compressed), "fasta")
            if self.input_file_type == "fastq":
                return SeqIO.parse(SadieIO.get_file_buffer(self.input, self.input_compressed), "fastq")
            if self.input_file_type == "abi":
                # A single instnace
                return SeqIO.parse(SadieIO.get_file_buffer(self.input, self.input_compressed), "abi")
            raise TypeError(f"Requested bad file type {self.input}")
        else:

            # these are directories which will chain together iterators
            if self.input_file_type == "fasta":
                return self._get_parse("fasta")
            if self.input_file_type == "fastq":
                return self._get_parse("fastq")
            if self.input_file_type == "abi":
                # has to be in rb
                return self._get_parse("abi", "rb")

    @staticmethod
    def get_file_type_dict(directory_path: Union[Path, str]) -> dict:
        """[summary]

        Parameters
        ----------
        directory_path : Union[Path, str]
            [description]

        Returns
        -------
        dict
            [description]

        Raises
        ------
        TypeError
            [description]
        """
        if isinstance(directory_path, str):
            directory_path = Path(directory_path)
        if not isinstance(directory_path, Path):
            raise TypeError(f"{directory_path} must be of type str or Path not {type(directory_path)}")
        if not directory_path.is_dir():
            raise TypeError(f"{directory_path} must be directory")
        if not directory_path.exists:
            return FileNotFoundError(f"{directory_path} does not exist")
        return {
            file: SadieIO.get_file_type(file, SadieIO.guess_input_compression(file))
            for file in glob.glob(str(directory_path) + "/*", recursive=True)
        }

    @staticmethod
    def guess_sequence_file_type(input_path: Union[Path, str]) -> str:
        """Guess Sequence File type.

        Returns
        -------
        filetype: str
            fasta, fastq or abi file format

        Raises
        ------
        NotImplementedError
            If a direcotry is passed as input and it contains heterogenous file types, we can not yet parse heterogenous file types
        """

        if isinstance(input_path, str):
            input_path = Path(input_path)
        if not isinstance(input_path, Path):
            raise TypeError(f"{input_path} is needs to be str of path or Path")

        if not input_path.is_dir():
            return SadieIO.get_file_type(input_path, SadieIO.guess_input_compression(input_path))

        # if dir
        else:
            _all_files_in_dir = SadieIO.get_file_type_dict(input_path)
            if len(set(_all_files_in_dir.values())) != 1:
                raise NotImplementedError(
                    f"all files in {input_path} directory are not of the same type {_all_files_in_dir}"
                )
            else:
                return list(_all_files_in_dir.values())[0]

    @staticmethod
    def get_file_buffer(file: Path, compression: str, mode="rt") -> Union[TextIOWrapper, gzip.GzipFile, bz2.BZ2File]:
        """Get an open file buffer

        Parameters
        ----------
        file : Path
            file path
        compression : str
            compression type, 'bzip2','gzip',None'
        mode : str, optional
            open mode, by default "rt"

        Returns
        -------
        Union[TextIOWrapper, gzip.GzipFile, bzip.BZ2File]
            Buffer of textfile

        Raises
        ------
        TypeError
            Can't determine file type
        """
        # get file buffer for compression
        if compression is None:
            file_buffer = open(file, mode)
        elif compression == "gzip":
            file_buffer = gzip.open(file, mode)
        elif compression == "bzip2":
            file_buffer = bz2.open(file, mode)
        else:
            raise TypeError(f"{file} can't determine file encoding")
        return file_buffer

    @staticmethod
    def get_file_type(file: Path, compression: str) -> str:
        """Get file type of file

        Parameters
        ----------
        file : Path
            A file type
        compression : str
            the compression type gzip, bzip2 or NOne

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

        if compression not in ["gzip", "bzip2", None]:
            raise TypeError("compression must be gzip,bzip2 or none")
        # now find bio type
        try:
            next(SeqIO.parse(SadieIO.get_file_buffer(file, compression), "fasta"))
            return "fasta"
        except Exception:
            pass
        try:
            next(SeqIO.parse(SadieIO.get_file_buffer(file, compression), "fastq"))
            return "fastq"
        except Exception:
            pass
        try:
            next(SeqIO.parse(SadieIO.get_file_buffer(file, compression, mode="rb"), "abi"))
            return "abi"
        except Exception:
            ValueError(f"can't determine file type of {file}")

    @staticmethod
    def guess_input_compression(input_path: str) -> Union[str, None]:
        """Guess if input compressed

        Parameters
        ----------
        input_path : str
            input is compressed

        Returns
        -------
        Tuple[str,None]
            returns filetype compression or none if uncompressed
        """
        if Path(input_path).is_dir():
            return "directory"

        _suffix = {"gz": "gzip", "bz2": "bzip2"}

        # filetype needs string
        _filetype = guess(str(input_path))
        if not _filetype:
            return None
        if _filetype.extension in _suffix.keys():
            return _suffix[_filetype.extension]
        raise ValueError(f"can't determine compression of {input_path}")

    def __repr__(self):
        return pformat(self.__dict__, indent=4)

    def __str__(self):
        return self.__repr__()
