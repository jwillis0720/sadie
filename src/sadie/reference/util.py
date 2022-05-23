import logging
from pathlib import Path
from typing import List, Union

from Bio.SeqRecord import SeqRecord

# reference logger
logger = logging.getLogger("Reference")


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
