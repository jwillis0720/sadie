import logging
import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import Union

logger = logging.getLogger(__name__)


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
