import logging
import os
from pathlib import Path
import platform
import subprocess
import shutil
import tempfile

logger = logging.getLogger("reference")


def write_blast_db(filename: str, output_db: str) -> None:
    """Take Input Fasta File

    Arguments:
        filename {str} -- the input fasta file
        output_db {str} -- the output path

    Raises:
        Exception: If makeblastdb fails for any reason
    """
    logger = logging.getLogger(__name__)
    
    # Validate Blast makeblastdb executable
    system = platform.system().lower()
    make_blast_db_exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"bin/{system}/makeblastdb")
    if not shutil.which(make_blast_db_exe):
        raise Exception(f"Make Blast DB {make_blast_db_exe} cant be found or is not executable")
    
    # Validate inpath fasta file and outpath directory
    filename = pathing(filename)
    output_db = pathing(output_db, new=True, overwrite=True)
    
    # Blast does not allow spaces in the path; need to use a temp file/dir
    with tempfile.NamedTemporaryFile(suffix=filename.suffix) as tmp_filename_obj:
        with tempfile.NamedTemporaryFile(suffix=None) as tmp_output_db_obj:
            # Only care for the path to the file
            tmp_filename = Path(tmp_filename_obj.name)
            tmp_output_db = Path(tmp_output_db_obj.name)
            # Write the fasta file to the temp file to garuntee no spaces in the path
            shutil.copy(filename, tmp_filename)
            # Execute makeblastdb
            logger.info(f"{make_blast_db_exe} -in {tmp_filename} -dbtype nucl -out {tmp_output_db}")
            make_blast_db = subprocess.run(
                [
                    make_blast_db_exe,
                    "-dbtype",
                    "nucl",
                    "-hash_index",
                    "-parse_seqids",
                    "-in",
                    tmp_filename,
                    "-out",
                    tmp_output_db,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            stdout = make_blast_db.stdout.decode("utf-8")
            stderr = make_blast_db.stderr.decode("utf-8")
            logger.debug((stdout))
            if make_blast_db.returncode:
                logger.critical("Problem with {}".format(filename))
                raise Exception(stderr)
            # We use the tmp file not as an actual output but to guarantee a unqiue prefix
            logger.info(f"Copying {tmp_output_db}.* to {output_db}")
            for tmpfile in tmp_output_db.parent.glob(f"{tmp_output_db.name}.*"):                
                logger.info(f"Copying {tmpfile} to {output_db.with_suffix(tmpfile.suffix)}")
                shutil.copy(tmpfile, output_db.with_suffix(tmpfile.suffix))
                tmpfile.unlink()  # we use the tmp file name only as a reference


def pathing(path: str, new: bool = False, overwrite: bool = True) -> Path:
    """Guarantees correct expansion rules for pathing.

    :param Union[str, Path] path: path of folder or file you wish to expand.
    :param bool new: will check if distination exists if new  (will check parent path regardless).
    :return: A pathlib.Path object.

    >>> pathing('~/Desktop/folderofgoodstuffs/')
    /home/user/Desktop/folderofgoodstuffs
    """
    path = Path(path)
    # Expand shortened path
    if str(path)[0] == "~":
        path = path.expanduser()
    # Exand local path
    if str(path)[0] == ".":
        path = path.resolve()
    else:
        path = path.absolute()
    # Making sure new paths don't exist while also making sure existing paths actually exist.
    if new:
        if not path.parent.exists():
            raise ValueError(f"ERROR ::: Parent directory of {path} does not exist.")
        if path.exists() and not overwrite:
            raise ValueError(f"ERROR ::: {path} already exists!")
    else:
        if not path.exists():
            raise ValueError(f"ERROR ::: Path {path} does not exist.")
    return path
