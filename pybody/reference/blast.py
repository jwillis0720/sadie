import logging
import subprocess


def write_blast_db(filename, output_db):
    """Take Input Fasta File

    Arguments:
        filename {str} -- the input fasta file 
        output_db {str} -- the output path 

    Raises:
        Exception: If makeblastdb fails for any reason
    """
    logger = logging.getLogger(__name__)
    make_blast_db = subprocess.run(
        [
            "makeblastdb",
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
    logger.debug((stdout))
    if make_blast_db.returncode:
        logger.critical("Problem with {}".format(filename))
        raise Exception(stderr)
