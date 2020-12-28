import bz2
import gzip
import os
from functools import partial
from mimetypes import guess_type
from pathlib import Path

from Bio import SeqIO as so


def get_project_root() -> Path:
    return Path(__file__).parent.parent


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


def determine_encoding(parent_file):
    encoding = guess_type(parent_file)
    if encoding[1] == "gzip":
        return "gzip", partial(gzip.open, mode="rt")
    if encoding[0] == "application/x-bzip" or encoding[1] == "bzip2":
        return "bzip", partial(bz2.open, mode="rt")
    return None, open


def split_fasta(parent_file, how_many, outdir="."):
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
        suffix = ".fasta.gz"
    elif encoding == "bzip":
        parent_file_base_name = ".".join(os.path.basename(parent_file).split(".")[0:-2])
        suffix = ".fasta.bz2"
    else:
        parent_file_base_name = os.path.basename(parent_file).split(".fasta")[0]
        suffix = ".fasta"

    file = parent_file_base_name + "_" + str(counter) + suffix
    file = os.path.join(outdir, file)
    # click.echo(f"Detected {encoding} filetype")
    # carries all of our records to be written
    joiner = []

    # _open will . handle all
    with _open(parent_file) as parent_file_handle:
        for num, record in enumerate(so.parse(parent_file_handle, "fasta"), start=1):

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
