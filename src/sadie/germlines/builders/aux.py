"""
Auxiliary File Builder
======================

Generates IgBLAST auxiliary files from gapped germline sequences.

IgBLAST auxiliary files for J genes contain (5 columns, tab-separated):
<gene_name>\t<reading_frame>\t<chain_type>\t<cdr3_end>\t<is_functional>

Example:
IGHJ1*01	0	JH	17	1
IGHJ2*01	1	JH	18	1
IGKJ1*01	1	JK	6	1
"""

import logging
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO

from sadie.germlines.builders.j_gene_data import get_j_gene_data

logger = logging.getLogger(__name__)


# Constants
CHAINS = ["H", "K", "L"]
# Only J segments are needed for IgBLAST aux files
SEGMENTS = ["J"]


class AuxFileBuilder:
    """
    Build IgBLAST auxiliary files from gapped sequences.

    Auxiliary files provide CDR/FWR boundaries for IgBLAST annotation.

    Examples
    --------
    >>> builder = AuxFileBuilder()
    >>> builder.build_for_species(
    ...     "human",
    ...     source_dir=Path("normalized/human/gapped"),
    ...     output_file=Path("igblast/aux_db/human_gl.aux")
    ... )
    """

    def build_for_species(self, species: str, source_dir: Path, output_file: Path) -> None:
        """
        Build auxiliary file for species.

        IgBLAST auxiliary files only contain J gene data.

        Parameters
        ----------
        species : str
            Species name
        source_dir : Path
            Directory with gapped FASTA files
        output_file : Path
            Output auxiliary file path
        """
        output_file.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Building auxiliary file for {species}")

        aux_lines = []

        # Process J segments only (IgBLAST aux files only need J genes)
        for chain in CHAINS:
            lines = self._process_segment(species, chain, "J", source_dir)  # Only J segments
            aux_lines.extend(lines)

        # Write auxiliary file
        if aux_lines:
            output_file.write_text("\n".join(aux_lines) + "\n")
            logger.info(f"Wrote {len(aux_lines)} entries to {output_file}")
        else:
            logger.warning(f"No auxiliary entries generated for {species}")

    def _process_segment(self, species: str, chain: str, segment: str, source_dir: Path) -> List[str]:
        """
        Process single segment to generate aux entries.

        Parameters
        ----------
        species : str
            Species name
        chain : str
            Chain type
        segment : str
            Segment type
        source_dir : Path
            Source directory

        Returns
        -------
        List[str]
            Auxiliary file lines for this segment
        """
        fasta_path = source_dir / f"IG{chain}{segment}.fasta"

        # Guard: file doesn't exist
        if not fasta_path.exists():
            logger.debug(f"No file: {fasta_path}")
            return []

        aux_lines = []

        try:
            records = list(SeqIO.parse(fasta_path, "fasta"))
        except Exception as e:
            logger.error(f"Failed to parse {fasta_path}: {e}")
            return []

        for record in records:
            aux_line = self._create_aux_entry(record, chain, segment)
            if aux_line:
                aux_lines.append(aux_line)

        logger.info(f"Generated {len(aux_lines)} aux entries from {fasta_path.name}")

        return aux_lines

    def _create_aux_entry(self, record, chain: str, segment: str) -> Optional[str]:
        """
        Create auxiliary file entry for a J gene sequence.

        IgBLAST aux format for J genes (5 columns, tab-separated):
        <gene_name>\t<reading_frame>\t<chain_type>\t<cdr3_end>\t<is_functional>

        Parameters
        ----------
        record : SeqRecord
            Sequence record
        chain : str
            Chain type (H, K, L)
        segment : str
            Segment type (must be "J")

        Returns
        -------
        str or None
            Auxiliary file line for J gene, None for other segments
        """
        if segment != "J":
            # Only J genes go in aux file
            return None

        gene_name = record.id

        # Get reference data for this J gene
        reading_frame, chain_type, cdr3_end, is_functional = get_j_gene_data(gene_name, chain)

        return f"{gene_name}\t{reading_frame}\t{chain_type}\t{cdr3_end}\t{is_functional}"

    def validate_aux_file(self, aux_file: Path) -> bool:
        """
        Validate auxiliary file format.

        Parameters
        ----------
        aux_file : Path
            Auxiliary file path

        Returns
        -------
        bool
            True if valid
        """
        if not aux_file.exists():
            logger.error(f"Aux file doesn't exist: {aux_file}")
            return False

        try:
            lines = aux_file.read_text().strip().split("\n")

            if not lines:
                logger.error("Aux file is empty")
                return False

            # Validate J gene format (5 columns)
            for line in lines:
                fields = line.split("\t")
                if len(fields) != 5:
                    logger.error(f"Invalid aux line (expected 5 columns): {line}")
                    return False

            logger.info(f"Aux file valid: {len(lines)} entries")
            return True

        except Exception as e:
            logger.error(f"Aux file validation failed: {e}")
            return False
