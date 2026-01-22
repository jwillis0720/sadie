"""
HMM Builder
===========

Generates Stockholm alignment files for HMMER hmmbuild from gapped germline sequences.

Per FR-013a-c:
- Accepts gapped V/J sequences
- Generates Stockholm format alignment files
- Minimum 3 sequences per segment/chain, maximum 100
"""

import logging
from pathlib import Path
from typing import List, Tuple, Optional
from Bio import SeqIO

logger = logging.getLogger(__name__)


class HMMBuilder:
    """
    Build Stockholm alignment files for HMMER from gapped germline sequences.

    Examples
    --------
    >>> builder = HMMBuilder()
    >>> builder.build_for_species(
    ...     "human",
    ...     source_dir=Path("normalized/human/gapped"),
    ...     output_dir=Path("hmm/human")
    ... )
    """

    MIN_SEQUENCES = 3
    MAX_SEQUENCES = 100

    def build_for_species(
        self,
        species: str,
        source_dir: Path,
        output_dir: Path
    ) -> None:
        """
        Build Stockholm alignments for all segments.

        Parameters
        ----------
        species : str
            Species name
        source_dir : Path
            Directory with gapped FASTA files
        output_dir : Path
            Output directory for Stockholm files
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        for chain in ["H", "K", "L"]:
            for segment in ["V", "J"]:
                self._build_segment_alignment(
                    species, chain, segment, source_dir, output_dir
                )

    def _build_segment_alignment(
        self,
        species: str,
        chain: str,
        segment: str,
        source_dir: Path,
        output_dir: Path
    ) -> None:
        """
        Build Stockholm alignment for single segment.

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
        output_dir : Path
            Output directory
        """
        fasta_path = source_dir / f"IG{chain}{segment}.fasta"

        if not fasta_path.exists():
            logger.debug(f"No file: {fasta_path}")
            return

        sequences = self._load_sequences(fasta_path)

        if len(sequences) < self.MIN_SEQUENCES:
            logger.warning(
                f"Insufficient sequences for {chain}{segment}: "
                f"{len(sequences)} < {self.MIN_SEQUENCES}"
            )
            return

        if len(sequences) > self.MAX_SEQUENCES:
            sequences = sequences[:self.MAX_SEQUENCES]
            logger.info(f"Truncated to {self.MAX_SEQUENCES} sequences for {chain}{segment}")

        output_file = output_dir / f"{species}_{chain}{segment}.sto"
        self._write_stockholm(sequences, output_file, f"{species}_{chain}{segment}")

        logger.info(f"Wrote {len(sequences)} sequences to {output_file}")

    def _load_sequences(self, fasta_path: Path) -> List[Tuple[str, str]]:
        """
        Load sequences from FASTA file.

        Parameters
        ----------
        fasta_path : Path
            Path to FASTA file

        Returns
        -------
        List[Tuple[str, str]]
            List of (name, sequence) tuples
        """
        sequences = []
        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                name = record.id.split("|")[0]
                seq = str(record.seq).upper()
                sequences.append((name, seq))
        except Exception as e:
            logger.error(f"Failed to parse {fasta_path}: {e}")

        return sequences

    def _write_stockholm(
        self,
        sequences: List[Tuple[str, str]],
        output_file: Path,
        alignment_id: str
    ) -> None:
        """
        Write sequences in Stockholm format.

        Parameters
        ----------
        sequences : List[Tuple[str, str]]
            List of (name, sequence) tuples
        output_file : Path
            Output file path
        alignment_id : str
            Alignment identifier
        """
        lines = ["# STOCKHOLM 1.0", f"#=GF ID {alignment_id}", ""]

        max_name_len = max(len(name) for name, _ in sequences)

        for name, seq in sequences:
            sto_seq = seq.replace(".", "-")
            lines.append(f"{name.ljust(max_name_len)}  {sto_seq}")

        lines.append("//")

        output_file.write_text("\n".join(lines) + "\n")


def get_gapped_sequences(
    manager,
    species: str,
    segment: str
) -> List[Tuple[str, str]]:
    """
    Get gapped sequences from GermlineManager.

    Per FR-013b: Returns List[Tuple[gene_name, gapped_sequence]]

    Parameters
    ----------
    manager : GermlineManager
        Germline manager instance
    species : str
        Species name
    segment : str
        Segment type (V or J)

    Returns
    -------
    List[Tuple[str, str]]
        List of (gene_name, gapped_sequence) tuples
    """
    result = []

    for chain in ["H", "K", "L"]:
        genes = manager.get_genes(species, segment, chain, functional_only=True)
        for gene in genes:
            if gene.sequence_gapped:
                result.append((gene.name, gene.sequence_gapped))

    return result
