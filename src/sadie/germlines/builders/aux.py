"""
Auxiliary File Builder
======================

Generates IgBLAST auxiliary files from gapped germline sequences.

IgBLAST auxiliary files contain:
- CDR and FWR region boundaries
- Chain type annotations
- Sequence orientation

Format (tab-separated):
<gene_name>\t<FWR1_start>\t<FWR1_end>\t<CDR1_start>\t<CDR1_end>...
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional
from Bio import SeqIO

logger = logging.getLogger(__name__)


# Constants
CHAINS = ["H", "K", "L"]
SEGMENTS = ["V", "J"]

IMGT_REGIONS_V = {
    "FWR1": (1, 26),
    "CDR1": (27, 38),
    "FWR2": (39, 55),
    "CDR2": (56, 65),
    "FWR3": (66, 104),
}

IMGT_REGIONS_J = {
    "FWR4": (1, 13),
}


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

    def build_for_species(
        self,
        species: str,
        source_dir: Path,
        output_file: Path
    ) -> None:
        """
        Build auxiliary file for species.

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

        # Process V and J segments (D doesn't have CDR/FWR)
        for chain in CHAINS:
            for segment in SEGMENTS:
                lines = self._process_segment(
                    species,
                    chain,
                    segment,
                    source_dir
                )
                aux_lines.extend(lines)

        # Write auxiliary file
        if aux_lines:
            output_file.write_text("\n".join(aux_lines) + "\n")
            logger.info(
                f"Wrote {len(aux_lines)} entries to {output_file}"
            )
        else:
            logger.warning(f"No auxiliary entries generated for {species}")

    def _process_segment(
        self,
        species: str,
        chain: str,
        segment: str,
        source_dir: Path
    ) -> List[str]:
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

        logger.info(
            f"Generated {len(aux_lines)} aux entries from {fasta_path.name}"
        )

        return aux_lines

    def _create_aux_entry(
        self,
        record,
        chain: str,
        segment: str
    ) -> Optional[str]:
        """
        Create auxiliary file entry for a single sequence.

        IgBLAST aux format:
        <gene_name>\t<FWR1_start>\t<FWR1_end>\t<CDR1_start>\t...

        Parameters
        ----------
        record : SeqRecord
            Sequence record
        chain : str
            Chain type
        segment : str
            Segment type

        Returns
        -------
        str or None
            Auxiliary file line
        """
        gene_name = record.id
        gapped_seq = str(record.seq)

        regions = self._parse_imgt_regions(gapped_seq, segment)
        if not regions:
            logger.debug(f"Could not parse regions for {gene_name}")
            return None

        if segment == "V":
            fields = [gene_name]
            for region in ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3"]:
                if region in regions:
                    start, end = regions[region]
                    fields.extend([str(start), str(end)])
                else:
                    fields.extend(["0", "0"])
            return "\t".join(fields)
        elif segment == "J":
            if "FWR4" in regions:
                start, end = regions["FWR4"]
                return f"{gene_name}\t{start}\t{end}"
            return None

        return None

    def _parse_imgt_regions(
        self,
        gapped_sequence: str,
        segment: str
    ) -> Dict[str, tuple]:
        """
        Parse CDR/FWR regions from IMGT-gapped sequence.

        Converts IMGT positions to ungapped sequence positions
        by counting non-gap characters.

        Parameters
        ----------
        gapped_sequence : str
            IMGT-gapped sequence (with dots)
        segment : str
            Segment type (V or J)

        Returns
        -------
        Dict[str, tuple]
            Region boundaries {region_name: (start, end)}
        """
        if segment == "V":
            imgt_regions = IMGT_REGIONS_V
        elif segment == "J":
            imgt_regions = IMGT_REGIONS_J
        else:
            return {}

        position_map = self._build_position_map(gapped_sequence)
        seq_len = len(gapped_sequence.replace(".", "").replace("-", ""))

        regions = {}
        for region_name, (imgt_start, imgt_end) in imgt_regions.items():
            ungapped_start = position_map.get(imgt_start)
            ungapped_end = position_map.get(imgt_end)

            if ungapped_start is not None and ungapped_end is not None:
                if ungapped_start <= seq_len and ungapped_end <= seq_len:
                    regions[region_name] = (ungapped_start, ungapped_end)

        return regions

    def _build_position_map(self, gapped_sequence: str) -> Dict[int, int]:
        """
        Map IMGT gapped positions to ungapped positions.

        Parameters
        ----------
        gapped_sequence : str
            IMGT-gapped sequence

        Returns
        -------
        Dict[int, int]
            Mapping from IMGT position (1-based) to ungapped position
        """
        position_map = {}
        ungapped_pos = 0

        for gapped_pos, char in enumerate(gapped_sequence, start=1):
            if char not in (".", "-"):
                ungapped_pos += 1
                position_map[gapped_pos] = ungapped_pos

        return position_map

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

            # Basic validation
            for line in lines:
                fields = line.split("\t")
                if len(fields) < 4:
                    logger.error(f"Invalid aux line: {line}")
                    return False

            logger.info(f"Aux file valid: {len(lines)} entries")
            return True

        except Exception as e:
            logger.error(f"Aux file validation failed: {e}")
            return False
