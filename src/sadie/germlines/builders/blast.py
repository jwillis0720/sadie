"""
BLAST Database Builder
======================

Builds BLAST databases for IgBLAST from germline FASTA files.

IgBLAST requires:
- Separate BLAST databases for V, D, J segments
- Standard BLAST database format (created by makeblastdb)
- Specific naming convention: <species>_<segment>

Design Principles:
- Use subprocess for makeblastdb (explicit over implicit)
- Validate inputs before building
- Provide clear error messages
"""

import logging
import subprocess
from pathlib import Path
from typing import Optional
from Bio import SeqIO

logger = logging.getLogger(__name__)


# Constants
SEGMENTS = ["V", "D", "J"]
CHAINS = ["H", "K", "L"]


class BlastDBBuilder:
    """
    Build BLAST databases for IgBLAST.

    Creates separate BLAST databases for each segment type (V, D, J)
    from normalized ungapped FASTA files.

    Examples
    --------
    >>> builder = BlastDBBuilder()
    >>> builder.build_for_species(
    ...     "human",
    ...     source_dir=Path("normalized/human/ungapped"),
    ...     output_dir=Path("igblast/database/human")
    ... )
    """

    # Maximum sequence ID length for BLAST databases
    MAX_SEQ_ID_LENGTH = 50

    def _sanitize_seq_id(self, seq_id: str) -> str:
        """
        Sanitize sequence ID for BLAST database compatibility.

        BLAST databases have a 50-character limit for local IDs when using
        -parse_seqids. VDJbase sequences often have very long IDs with variant
        annotations that exceed this limit.

        Parameters
        ----------
        seq_id : str
            Original sequence ID

        Returns
        -------
        str
            Sanitized sequence ID (max 50 chars)
        """
        # Remove source annotation if present (e.g., "source=vdjbase")
        if " source=" in seq_id:
            seq_id = seq_id.split(" source=")[0]

        # If still too long, truncate to max length
        if len(seq_id) > self.MAX_SEQ_ID_LENGTH:
            # Try to keep the base gene name (before variant annotations)
            # VDJbase format: IGHV1-18*04_g107c_a110t_...
            parts = seq_id.split("_")
            if len(parts) > 1:
                # Keep gene name and as many variants as fit
                base = parts[0]  # e.g., IGHV1-18*04
                truncated = base
                for part in parts[1:]:
                    if len(truncated) + 1 + len(part) <= self.MAX_SEQ_ID_LENGTH:
                        truncated = f"{truncated}_{part}"
                    else:
                        break
                seq_id = truncated
            else:
                # Just truncate if no underscores
                seq_id = seq_id[:self.MAX_SEQ_ID_LENGTH]

        return seq_id

    def build_for_species(
        self,
        species: str,
        source_dir: Path,
        output_dir: Path
    ) -> None:
        """
        Build all BLAST databases for a species.

        Parameters
        ----------
        species : str
            Species name
        source_dir : Path
            Directory containing ungapped FASTA files (normalized/)
        output_dir : Path
            Output directory for BLAST databases
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Building BLAST databases for {species}")

        # Combine all chains for each segment
        for segment in SEGMENTS:
            self._build_segment_database(
                species,
                segment,
                source_dir,
                output_dir
            )

    def _build_segment_database(
        self,
        species: str,
        segment: str,
        source_dir: Path,
        output_dir: Path
    ) -> None:
        """
        Build BLAST database for single segment.

        Combines all chains (H, K, L) for this segment.

        Parameters
        ----------
        species : str
            Species name
        segment : str
            Segment type (V, D, or J)
        source_dir : Path
            Source directory with FASTA files
        output_dir : Path
            Output directory
        """
        # Collect sequences from all chains
        combined_sequences = []
        seen_ids: dict[str, int] = {}  # Track duplicate IDs

        for chain in CHAINS:
            fasta_path = source_dir / f"IG{chain}{segment}.fasta"

            if not fasta_path.exists():
                logger.debug(f"No file: {fasta_path}")
                continue

            try:
                records = list(SeqIO.parse(fasta_path, "fasta"))
                # Sanitize sequence IDs for BLAST (max 50 chars)
                for record in records:
                    sanitized_id = self._sanitize_seq_id(record.id)
                    # Handle duplicates by appending counter
                    if sanitized_id in seen_ids:
                        seen_ids[sanitized_id] += 1
                        sanitized_id = f"{sanitized_id}_{seen_ids[sanitized_id]}"
                    else:
                        seen_ids[sanitized_id] = 0
                    record.id = sanitized_id
                    record.description = ""  # Clear description to avoid issues
                combined_sequences.extend(records)
                logger.info(
                    f"Added {len(records)} sequences from {fasta_path.name}"
                )
            except Exception as e:
                logger.error(f"Failed to read {fasta_path}: {e}")

        # Guard: no sequences found
        if not combined_sequences:
            logger.warning(f"No sequences for {species} {segment}")
            return

        # Write combined FASTA
        combined_fasta = output_dir / f"{species}_{segment}.fasta"
        SeqIO.write(combined_sequences, combined_fasta, "fasta")
        logger.info(
            f"Wrote {len(combined_sequences)} sequences to {combined_fasta}"
        )

        # Build BLAST database
        self._run_makeblastdb(combined_fasta, segment)

    def _run_makeblastdb(self, fasta_path: Path, segment: str) -> None:
        """
        Run makeblastdb to create BLAST database.

        Parameters
        ----------
        fasta_path : Path
            Input FASTA file
        segment : str
            Segment type (for logging)
        """
        db_name = fasta_path.stem

        try:
            cmd = [
                "makeblastdb",
                "-dbtype", "nucl",
                "-in", str(fasta_path),
                "-out", str(fasta_path.with_suffix("")),
                "-parse_seqids",
                "-hash_index",
            ]

            logger.debug(f"Running: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            logger.info(f"Built BLAST database: {db_name}")
            logger.debug(result.stdout)

        except subprocess.CalledProcessError as e:
            logger.error(
                f"makeblastdb failed for {db_name}: {e.stderr}"
            )
            raise

        except FileNotFoundError:
            logger.error(
                "makeblastdb not found. "
                "Ensure BLAST+ is installed and in PATH."
            )
            raise

    def validate_database(self, db_path: Path) -> bool:
        """
        Validate BLAST database was created successfully.

        Parameters
        ----------
        db_path : Path
            Path to BLAST database (without extension)

        Returns
        -------
        bool
            True if valid
        """
        # Check for required BLAST database files
        required_extensions = [".nhr", ".nin", ".nsq"]

        for ext in required_extensions:
            file_path = db_path.with_suffix(ext)
            if not file_path.exists():
                logger.error(f"Missing BLAST file: {file_path}")
                return False

        return True
