"""
Auxiliary File Builder
======================

Generates IgBLAST auxiliary files from gapped germline sequences.

IgBLAST auxiliary files contain two formats:

V genes (10 columns, tab-separated):
<gene_name>\t<reading_frame>\t<fwr1_end>\t<cdr1_start>\t<cdr1_end>\t<fwr2_start>\t<fwr2_end>\t<cdr2_start>\t<cdr2_end>\t<fwr3_end>

Where all positions are amino acid positions (1-based) in the ungapped sequence.
Truncated regions use 0 for missing end positions.

Example:
IGHV1-2*01	1	26	27	38	39	55	56	65	104

J genes (5 columns, tab-separated):
<gene_name>\t<reading_frame>\t<chain_type>\t<cdr3_end>\t<extra_bps>

Example:
IGHJ1*01	0	JH	17	1
"""

import logging
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO

from sadie.germlines.builders.imgt_positions import calculate_imgt_v_positions
from sadie.germlines.builders.j_gene_data import get_j_gene_data

logger = logging.getLogger(__name__)


# Constants
CHAINS = ["H", "K", "L"]
# V and J segments for aux files (V for region boundaries, J for CDR3 end)
V_SEGMENTS = ["V"]
J_SEGMENTS = ["J"]


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

        IgBLAST auxiliary files contain both V-gene region boundaries
        and J-gene CDR3 end positions.

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

        # Process V segments first (10-column format for region boundaries)
        for chain in CHAINS:
            lines = self._process_v_segment(species, chain, source_dir)
            aux_lines.extend(lines)

        # Process J segments (5-column format for CDR3 end)
        for chain in CHAINS:
            lines = self._process_j_segment(species, chain, source_dir)
            aux_lines.extend(lines)

        # Write auxiliary file
        if aux_lines:
            output_file.write_text("\n".join(aux_lines) + "\n")
            logger.info(f"Wrote {len(aux_lines)} entries to {output_file}")
        else:
            logger.warning(f"No auxiliary entries generated for {species}")

    def _process_v_segment(self, species: str, chain: str, source_dir: Path) -> List[str]:
        """
        Process V segment to generate aux entries with region boundaries.

        Parameters
        ----------
        species : str
            Species name
        chain : str
            Chain type (H, K, L)
        source_dir : Path
            Source directory with gapped FASTA files

        Returns
        -------
        List[str]
            Auxiliary file lines (10-column format) for V genes
        """
        fasta_path = source_dir / f"IG{chain}V.fasta"

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
            aux_line = self._create_v_aux_entry(record)
            if aux_line:
                aux_lines.append(aux_line)

        logger.info(f"Generated {len(aux_lines)} V-gene aux entries from {fasta_path.name}")
        return aux_lines

    def _process_j_segment(self, species: str, chain: str, source_dir: Path) -> List[str]:
        """
        Process J segment to generate aux entries with CDR3 end positions.

        Parameters
        ----------
        species : str
            Species name
        chain : str
            Chain type (H, K, L)
        source_dir : Path
            Source directory with gapped FASTA files

        Returns
        -------
        List[str]
            Auxiliary file lines (5-column format) for J genes
        """
        fasta_path = source_dir / f"IG{chain}J.fasta"

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
            aux_line = self._create_j_aux_entry(record, chain, species)
            if aux_line:
                aux_lines.append(aux_line)

        logger.info(f"Generated {len(aux_lines)} J-gene aux entries from {fasta_path.name}")
        return aux_lines

    def _create_v_aux_entry(self, record) -> Optional[str]:
        """
        Create auxiliary file entry for a V gene sequence.

        IgBLAST aux format for V genes (10 columns, tab-separated):
        <gene_name>\t<reading_frame>\t<fwr1_end>\t<cdr1_start>\t<cdr1_end>\t<fwr2_start>\t<fwr2_end>\t<cdr2_start>\t<cdr2_end>\t<fwr3_end>

        All positions are amino acid positions (1-based) in the ungapped sequence.
        Truncated regions use 0 for missing positions.

        Parameters
        ----------
        record : SeqRecord
            Gapped V gene sequence record

        Returns
        -------
        str or None
            Auxiliary file line for V gene, None if positions cannot be determined
        """
        gene_name = record.id
        gapped_sequence = str(record.seq)

        # Calculate nucleotide positions from gapped sequence (0-based)
        nt_positions = calculate_imgt_v_positions(gapped_sequence, zero_based=True)

        if not nt_positions:
            logger.debug(f"Could not calculate positions for {gene_name}")
            return None

        # Convert nucleotide positions to amino acid positions (1-based)
        # Formula: aa_pos = (nt_pos // 3) + 1
        def nt_to_aa(nt_pos: int) -> int:
            """Convert 0-based nucleotide position to 1-based amino acid position."""
            return (nt_pos // 3) + 1

        # Extract region boundaries, using 0 for missing regions
        fwr1_end = nt_to_aa(nt_positions.get("fwr1", (0, 0))[1]) if "fwr1" in nt_positions else 0
        cdr1_start = nt_to_aa(nt_positions.get("cdr1", (0, 0))[0]) if "cdr1" in nt_positions else 0
        cdr1_end = nt_to_aa(nt_positions.get("cdr1", (0, 0))[1]) if "cdr1" in nt_positions else 0
        fwr2_start = nt_to_aa(nt_positions.get("fwr2", (0, 0))[0]) if "fwr2" in nt_positions else 0
        fwr2_end = nt_to_aa(nt_positions.get("fwr2", (0, 0))[1]) if "fwr2" in nt_positions else 0
        cdr2_start = nt_to_aa(nt_positions.get("cdr2", (0, 0))[0]) if "cdr2" in nt_positions else 0
        cdr2_end = nt_to_aa(nt_positions.get("cdr2", (0, 0))[1]) if "cdr2" in nt_positions else 0
        fwr3_end = nt_to_aa(nt_positions.get("fwr3", (0, 0))[1]) if "fwr3" in nt_positions else 0

        # Reading frame is typically 1 for V genes (starts at first codon)
        reading_frame = 1

        return f"{gene_name}\t{reading_frame}\t{fwr1_end}\t{cdr1_start}\t{cdr1_end}\t{fwr2_start}\t{fwr2_end}\t{cdr2_start}\t{cdr2_end}\t{fwr3_end}"

    def _create_j_aux_entry(self, record, chain: str, species: str) -> Optional[str]:
        """
        Create auxiliary file entry for a J gene sequence.

        IgBLAST aux format for J genes (5 columns, tab-separated):
        <gene_name>\t<reading_frame>\t<chain_type>\t<cdr3_end>\t<extra_bps>

        Parameters
        ----------
        record : SeqRecord
            J gene sequence record
        chain : str
            Chain type (H, K, L)
        species : str
            Species name for motif lookup

        Returns
        -------
        str or None
            Auxiliary file line for J gene
        """
        gene_name = record.id
        sequence = str(record.seq).replace(".", "").replace("-", "").upper()

        # Get reference data for this J gene (uses species-specific motif patterns)
        reading_frame, chain_type, cdr3_end, extra_bps = get_j_gene_data(
            gene_name, chain, sequence, species=species
        )

        return f"{gene_name}\t{reading_frame}\t{chain_type}\t{cdr3_end}\t{extra_bps}"

    def validate_aux_file(self, aux_file: Path) -> bool:
        """
        Validate auxiliary file format.

        Aux files contain two formats:
        - V genes: 10 columns (region boundaries)
        - J genes: 5 columns (CDR3 end)

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

            v_count = 0
            j_count = 0

            for line in lines:
                fields = line.split("\t")
                gene_name = fields[0] if fields else ""

                # V genes have 10 columns, J genes have 5 columns
                if "V" in gene_name and len(fields) == 10:
                    v_count += 1
                elif "J" in gene_name and len(fields) == 5:
                    j_count += 1
                else:
                    logger.warning(f"Unexpected aux line format: {line}")

            logger.info(f"Aux file valid: {v_count} V-gene entries, {j_count} J-gene entries")
            return True

        except Exception as e:
            logger.error(f"Aux file validation failed: {e}")
            return False
