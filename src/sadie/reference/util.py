import logging
import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import List, Union

import pandas as pd
from Bio.SeqRecord import SeqRecord

# reference logger
logger = logging.getLogger("Reference")


def make_blast_db_for_internal(df: pd.DataFrame, dboutput: Union[str, Path]) -> None:
    """Make a blast database from dataframe"""
    out_fasta = Path(dboutput).with_suffix(".fasta")
    logger.debug("Writing fasta to {}".format(out_fasta))
    with open(out_fasta, "w") as f:
        for id_, seq in zip(df["gene"], df["sequence"]):
            f.write(">{}\n{}\n".format(id_, seq))
    # get basename
    write_blast_db(out_fasta, dboutput)


def write_blast_db(
    filename: Union[Path, str], output_db: Union[Path, str], make_blast_db_bin: Union[Path, None, str] = None
) -> None:
    """Write input fasta to a blast database

    Parameters
    ----------
    filename : Union[Path, str]
        The fasta file
    output_db : Union[Path, str]
        The output blast database path
    make_blast_db_bin : Union[Path, str], optional
        if there is a custom makeblastdb binary to use

    Raises
    ------
    ValueError
        If makeblastdb is not found
    Exception
    """

    system = platform.system().lower()
    if isinstance(filename, str):
        filename = Path(filename)
    if not make_blast_db_bin:
        make_blast_db_bin = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"bin/{system}/makeblastdb")
    if not shutil.which(make_blast_db_bin):
        raise ValueError(f"Make Blast DB {make_blast_db_bin} cant be found or is not executable")
    make_blast_db = subprocess.run(
        [
            make_blast_db_bin,
            "-dbtype",
            "nucl",
            "-hash_index",
            "-parse_seqids",
            "-in",
            filename,
            "-out",
            output_db,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout = make_blast_db.stdout.decode("utf-8")
    stderr = make_blast_db.stderr.decode("utf-8")
    logger.info(stdout)
    if make_blast_db.returncode:
        logger.error("Problem with {}".format(filename))
        logger.error(stderr)
        raise RuntimeError(stderr)


def write_out_fasta(sequences: List[SeqRecord], outpath: Union[str, Path]) -> Path:
    """Conveinence function to write out a list of SeqRecords to a fasta file

    Parameters
    ----------
    sequences : List[SeqRecord]
        List of SeqRecord objects from Biopython
    outpath : Path
        Path to write fasta to

    Returns
    -------
    Path
        The fully qualified path
    """
    logger = logging.getLogger(__file__)
    if isinstance(outpath, str):
        outpath = Path(outpath)
    output_fasta = outpath.with_suffix(".fasta")
    logger.debug("output fasta {}".format(output_fasta))

    with open(output_fasta, "w") as f:
        for sequence in sequences:
            name = sequence.name
            seq = str(sequence.seq).replace(".", "")
            f.write(f">{name}\n{seq}\n")
    return output_fasta
