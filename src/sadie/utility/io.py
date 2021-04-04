from pathlib import Path
from Bio import SeqIO
from filetype import guess
from typing import Tuple, List
from pprint import pformat
import glob
import gzip
import bz2
from abifpy import Trace


class Abi:
    def __init__(self, abi_file: Path):
        # Use abipy library to parse abi_file
        self.abi_trace_object = Trace(abi_file)

        # Resolve name
        self.abi_name = self.abi_trace_object.name

        # Resolve sequecnce base call
        self.abi_seq = self.abi_trace_object.seq

        # Resolve phred call
        self.abi_quality = self.abi_trace_object.qual_val


class AbiFiles:
    def __init__(self, abi_files: List):
        if not all([isinstance(x, Abi) for x in abi_files]):
            raise TypeError(f"{abi_files} must all be Abi")
        self.all_abi_files = abi_files

    def __add__(self, other: "AbiFiles"):
        if not isinstance(other, AbiFiles):
            raise TypeError(f"{other} must be AbiFiles")
        self.all_abi_files += other.all_abi_files

    def __iter__(self):
        for x in self.all_abi_files:
            yield x

    def append(self, abi: Abi):
        if not isinstance(abi, Abi):
            raise TypeError(f"{abi} must be Abi")
        self.all_abi_files.append(abi)


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
            the compression type of the output, "None","gzip","bzip"
        """

        # input path
        self.input = Path(input_path)
        if not self.input.exists:
            raise FileNotFoundError(f"{self.input} not found")

        # we assume input is not a directorys
        self.isdir = False

        # check if input is a directory:
        if self.input.is_dir():
            self.isdir = True

        # check if input is compressed - gzip or bzip2
        self.input_compressed = self.guess_input_compression(self.input)

        # check if input format was specified
        self.infer_input = False
        if in_format == "infer":
            self.infer_input = True
        if self.infer_input:
            self.input_file_type = self.guess_file_type()
        else:
            self.input_file_type = in_format

        # output
        self.output = output_path

    def guess_file_type(self):
        """Guess File Type of Input

        Returns
        -------
        filetype
            fasta, fastq or abi

        Raises
        ------
        NotImplementedError
            If a direcotry is passed as input, we can not yet parse heterogenous file types
        """
        if not self.isdir:
            return self._get_file_type(self.input, self.input_compressed)

        # if dir
        else:
            self.all_files_in_dir = {
                file: self._get_file_type(file, self.guess_input_compression(file))
                for file in glob.glob(str(self.input) + "/*", recursive=True)
            }
            if len(set(self.all_files_in_dir.values())) != 1:
                raise NotImplementedError(
                    f"all files in {self.input} directory are not of the same type {self.all_files_in_dir}"
                )
            else:
                return list(self.all_files_in_dir.values())[0]

    @staticmethod
    def _get_file_type(file: Path, compression: str) -> str:
        """Get file type of only file

        Parameters
        ----------
        file : Path
            A file type
        compression : str
            the compression type gzip, bzip2 or NOne

        Returns
        -------
        str
            the sequence file type
        Raises
        ------
        TypeError
            Can't determine file or compression
        """
        file = Path(file)

        def _get_file_buffer(mode="rt"):
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

        # now find bio type
        try:
            next(SeqIO.parse(_get_file_buffer(), "fasta"))
            return "fasta"
        except Exception:
            pass
        try:
            next(SeqIO.parse(_get_file_buffer(), "fastq"))
            return "fastq"
        except Exception:
            pass
        try:
            next(SeqIO.parse(_get_file_buffer(mode="rb"), "abi"))
            return "abi"
        except Exception:
            ValueError(f"can't determine file type of {file}")

    @staticmethod
    def guess_input_compression(input_path: str) -> Tuple[str, None]:
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
            return None

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
