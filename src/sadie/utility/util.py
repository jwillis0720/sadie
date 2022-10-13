from __future__ import annotations

import bz2
import gzip
import logging
import os
import warnings
from pathlib import Path
from typing import Union

import pandas as pd
from Bio.SeqIO import MultipleSeqAlignment

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio import pairwise2, SeqIO

from sadie.utility.io import SadieInputFile

# global
blosum_matrix = matlist.blosum62
logger = logging.getLogger("Utilility")


def getVerbosityLevel(verbosity_count: int) -> int:
    """Set verbosity level by how many --vvv were passed

    Arguments:
        verboisty_count {int} -- how many v's were passed

    50 - critical
    40 - error
    30 - warning
    20 - info
    10 -debug
    0 - notset
    """

    # If 5, we want 10 debug level logging
    if verbosity_count == 5:
        return 10
    # If 4, we want 20 info level logging
    if verbosity_count == 4:
        return 20
    # If 3, we want 30 warming level logging
    if verbosity_count == 3:
        return 30
    # If 2, we want 40 error level logging
    if verbosity_count == 2:
        return 40

    # always return critical
    return 50


def get_project_root() -> Path:
    """Get the pakage root direcotry"""
    return Path(__file__).parent.parent


def format_alignment(
    alignment: MultipleSeqAlignment, max_line_length: int = 80, ljust: int = 12, max_id: int = 30
) -> str:
    """Format an alignment object

    Parameters
    ----------
    alignment : Bio.MultipleSequenceAlignment
        Biopython multiple sequence alignment
    max_line_length : int, optional
        how many characters until wrap, by default 80
    ljust : int, optional
        what is the spacing between id and sequence, by default 12
    max_id : int, optional
        maximum characters in id of alignment, by default 30

    Returns
    -------
    str
        Formatted string alignment

    Raises
    ------
    ValueError
        Must have at least one sequence
    ValueError
        If sequence is empty
    """
    if len(alignment) == 0:
        raise ValueError("Must have at least one sequence")
    if alignment.get_alignment_length() == 0:
        # This doubles as a check for an alignment object
        raise ValueError("Non-empty sequences are required")

    output = ""
    cur_char = 0
    max_length = len(alignment[0])

    if max_length <= 0:
        raise ValueError("Non-empty sequences are required")

    # keep displaying sequences until we reach the end
    while cur_char != max_length:
        # calculate the number of sequences to show, which will
        # be less if we are at the end of the sequence
        if (cur_char + max_line_length) > max_length:
            show_num = max_length - cur_char
        else:
            show_num = max_line_length

        # go through all of the records and print out the sequences
        # when we output, we do a nice 80 column output, although this
        # may result in truncation of the ids.
        for record in alignment:
            # Make sure we don't get any spaces in the record
            # identifier when output in the file by replacing
            # them with underscores:
            line = record.id[0:max_id].replace(" ", "_").ljust(ljust)
            line += str(record.seq[cur_char : (cur_char + show_num)])
            output += line + "\n"
        output += "\n"
        cur_char += show_num
    return output


def split_fasta(
    parent_file: Union[str, Path], how_many: int, outdir: Union[str, Path] = ".", filetype: str = "fasta"
) -> None:
    # files_suffix = {"fasta": ".fasta", "fastq": ".fastq"}
    # get file_counter and base name of fasta_file
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    counter: int = 1
    sadie_handle: SadieInputFile = SadieInputFile(parent_file)
    compression_format: str | None = sadie_handle.compression_format
    if not compression_format:
        compression_format = "Uncompressed"
    # click.echo(f"Detected {encoding} filetype")
    # our first file name
    if compression_format == "gzip":
        parent_file_base_name = ".".join(os.path.basename(parent_file).split(".")[0:-2])
        suffix = f".{filetype}.gz"
    elif compression_format == "bzip":
        parent_file_base_name = ".".join(os.path.basename(parent_file).split(".")[0:-2])
        suffix = f".{filetype}.bz2"
    else:
        parent_file_base_name = os.path.basename(parent_file).split(f".{filetype}")[0]
        suffix = f".{filetype}"

    file: str = parent_file_base_name + "_" + str(counter) + suffix
    file = os.path.join(outdir, file)
    # click.echo(f"Detected {encoding} filetype")
    # carries all of our records to be written
    joiner = []

    # _open will . handle all
    with sadie_handle.open_input as parent_file_handle:
        for num, record in enumerate(SeqIO.parse(parent_file_handle, f"{filetype}"), start=1):

            # append records to our list holder
            joiner.append(">" + record.id + "\n" + str(record.seq))

            # if we have reached the maximum numbers to be in that file, write to a file, and then clear
            # record holder
            if num % how_many == 0:
                joiner.append("")
                if compression_format == "gzip":
                    with gzip.open(file, "wb") as f:
                        # print(file)
                        f.write("\n".join(joiner).encode("utf-8"))
                    # print(file)
                elif compression_format == "bzip":
                    with bz2.open(file, "wb") as f:
                        f.write("\n".join(joiner).encode("utf-8"))
                else:
                    with open(file, "w") as f:
                        f.write("\n".join(joiner))

                # click.echo(f"wrote to {file}")
                # change file name,clear record holder, and change the file count
                counter += 1
                file = parent_file_base_name + "_" + str(counter) + suffix
                file = os.path.join(outdir, file)
                joiner = []
        if joiner:
            # click.echo(f"writing final file")
            joiner.append("")
            if compression_format == "gzip":
                with gzip.open(file, "wb") as f:
                    f.write("\n".join(joiner).encode("utf-8"))
            elif compression_format == "bzip":
                with bz2.open(file, "wb") as f:
                    f.write("\n".join(joiner).encode("utf-8"))
            else:
                with open(file, "w") as f:
                    f.write("\n".join(joiner))


def get_verbosity_level(verbosity_count: int) -> int:
    """Get verbosity level by how many --vvv were passed

    Arguments:
        verboisty_count {int} -- how many v's were passed

    50 - critical
    40 - error
    30 - warning
    20 - info
    10 -debug
    0 - notset
    """

    # If 5, we want 10 debug level logging
    if verbosity_count >= 5:
        return 10
    # If 4, we want 20 info level logging
    elif verbosity_count == 4:
        return 20
    # If 3, we want 30 warming level logging
    elif verbosity_count == 3:
        return 30
    # If 2, we want 40 error level logging
    elif verbosity_count == 2:
        return 40

    # always return critical
    return 50


def guess_compress(filename: Union[str, Path]) -> Union[str, bool]:
    """Guess compression type

    Arguments:
        filename {str} -- filepath

    Returns:
        str -- file suffix '.gz.'
    """
    magic_dict = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        "\x50\x4b\x03\x04": "zip",
    }
    max_len = max(len(x) for x in magic_dict)
    with open(filename, "r") as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return False


def is_tool(name: str) -> bool:
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which

    return which(name) is not None


def correct_alignment(X: pd.Series, field_1: str, field_2: str) -> pd.Series:  # type: ignore
    alignment_aa_1 = X[field_1]
    alignment_aa_2 = X[field_2]
    if any([isinstance(alignment_aa_1, float), isinstance(alignment_aa_2, float)]):
        return pd.Series({field_1: alignment_aa_1, field_2: alignment_aa_2})

    # get alignments
    try:
        alignments = pairwise2.align.globalds(
            alignment_aa_1,
            alignment_aa_2,
            blosum_matrix,
            -12,
            -4,
            penalize_extend_when_opening=True,
            penalize_end_gaps=False,
        )
    except SystemError:
        logger.debug(f"System error most likely due to * in germline alignment...falling back {alignment_aa_2}")
        alignments = pairwise2.align.globalms(
            alignment_aa_1,
            alignment_aa_2,
            4,
            -1,
            -12,
            -4,
            penalize_extend_when_opening=True,
            penalize_end_gaps=False,
        )
    alignment_1_aa_corrected = alignments[0][0]
    alignment_2_aa_corrected = alignments[0][1]
    return pd.Series({field_1: alignment_1_aa_corrected, field_2: alignment_2_aa_corrected})


def get_consensus_of_paired_end_abi(abi_file_1: Union[str, Path], abi_file_2: Union[str, Path]) -> str:
    """Get consensus of paired end ABI files. Get base of higher quality read.

    Parameters
    ----------
    abi_file_1 : Union[str, Path]
        Path to first ABI file. Forward read.

    abi_file_2 : [type]
        Path to second ABI file. Reverse read. Must be reverse compliment of abi_file_1.

    Returns
    -------
    str
        The consensus sequence

    """
    # get the reads
    reads_1 = SeqIO.parse(abi_file_1, "abi-trim")
    reads_2 = SeqIO.parse(abi_file_2, "abi-trim")
    # get the read
    read_1 = next(reads_1)
    read_2 = next(reads_2)
    # get the read's seq
    seq_1 = read_1.seq
    seq_2 = read_2.seq
    # get the first read's phred
    phred_1 = read_1.letter_annotations["phred_quality"]
    phred_2 = read_2.letter_annotations["phred_quality"][::-1]

    # get the alignment of the two seqs with high gap penalties
    alignments = pairwise2.align.globalms(
        seq_1,
        seq_2.reverse_complement(),
        4,
        -1,
        -12,
        -4,
        penalize_extend_when_opening=True,
        penalize_end_gaps=False,
    )

    # get top scoring
    top_alignment = alignments[0]

    # destructure
    seq_1_aligned, seq_2_aligned = top_alignment[0], top_alignment[1]

    # set indexes to zero
    phred_1_indexer, phred_2_indexer = 0, 0
    consensus_seq = ""
    for i in range(len(seq_1_aligned)):
        seq_1_position = seq_1_aligned[i]
        seq_2_position = seq_2_aligned[i]
        if seq_1_position == "-":
            consensus_seq += seq_2_position
            phred_2_indexer += 1
        elif seq_2_position == "-":
            consensus_seq += seq_1_position
            phred_1_indexer += 1
        elif seq_1_position == seq_2_position:
            consensus_seq += seq_1_aligned[i]
            phred_1_indexer += 1
            phred_2_indexer += 1
        elif seq_1_position != seq_2_position:
            if seq_1_position == "N" and seq_2_position != "N":
                consensus_seq += seq_2_position
            elif seq_2_position == "N" and seq_1_position != "N":
                consensus_seq += seq_1_position
            else:
                if phred_1[phred_1_indexer] > phred_2[phred_2_indexer]:
                    consensus_seq += seq_1_position
                else:
                    consensus_seq += seq_2_position
            phred_1_indexer += 1
            phred_2_indexer += 1

        else:
            raise ValueError(f"Should not have reached this block for positon {i}")
    return consensus_seq
