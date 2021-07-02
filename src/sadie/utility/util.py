import bz2
import gzip
import os
import logging
import warnings
from functools import partial
from mimetypes import guess_type
from pathlib import Path
from typing import Union

import pandas as pd

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio import pairwise2, SeqIO

blosum_matrix = matlist.blosum62
logger = logging.getLogger("Utilility")


def get_project_root() -> Path:
    return Path(__file__).parent.parent


def determine_encoding(parent_file):
    guess_compress(parent_file)
    encoding = guess_type(parent_file)
    if encoding[1] == "gzip":
        return "gzip", partial(gzip.open, mode="rt")
    if encoding[0] == "application/x-bzip" or encoding[1] == "bzip2":
        return "bzip", partial(bz2.open, mode="rt")
    return None, open


def split_fasta(parent_file, how_many, outdir=".", filetype="fasta"):
    # files_suffix = {"fasta": ".fasta", "fastq": ".fastq"}
    # get file_counter and base name of fasta_file
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    counter = 1
    encoding, _open = determine_encoding(parent_file)
    if not encoding:
        encoding = "Uncompressed"
    # click.echo(f"Detected {encoding} filetype")
    # our first file name
    if encoding == "gzip":
        parent_file_base_name = ".".join(os.path.basename(parent_file).split(".")[0:-2])
        suffix = f".{filetype}.gz"
    elif encoding == "bzip":
        parent_file_base_name = ".".join(os.path.basename(parent_file).split(".")[0:-2])
        suffix = f".{filetype}.bz2"
    else:
        parent_file_base_name = os.path.basename(parent_file).split(f".{filetype}")[0]
        suffix = f".{filetype}"

    file = parent_file_base_name + "_" + str(counter) + suffix
    file = os.path.join(outdir, file)
    # click.echo(f"Detected {encoding} filetype")
    # carries all of our records to be written
    joiner = []

    # _open will . handle all
    with _open(parent_file) as parent_file_handle:
        for num, record in enumerate(SeqIO.parse(parent_file_handle, f"{filetype}"), start=1):

            # append records to our list holder
            joiner.append(">" + record.id + "\n" + str(record.seq))

            # if we have reached the maximum numbers to be in that file, write to a file, and then clear
            # record holder
            if num % how_many == 0:
                joiner.append("")
                if encoding == "gzip":
                    with gzip.open(file, "wb") as f:
                        # print(file)
                        f.write("\n".join(joiner).encode("utf-8"))
                    # print(file)
                elif encoding == "bzip":
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
            if encoding == "gzip":
                with gzip.open(file, "wb") as f:
                    f.write("\n".join(joiner).encode("utf-8"))
            elif encoding == "bzip":
                with bz2.open(file, "wb") as f:
                    f.write("\n".join(joiner).encode("utf-8"))
            else:
                with open(file, "w") as f:
                    f.write("\n".join(joiner))


def get_verbosity_level(verbosity_count):
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


def guess_compress(filename):
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


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which

    return which(name) is not None


def correct_alignment(X: pd.Series, field_1: str, field_2: str) -> pd.Series:
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
