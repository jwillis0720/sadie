"""
Sequence Gapper - IMGT Numbering Alignment
==========================================

Gaps ungapped germline sequences to IMGT numbering using BioPython
pairwise alignment against IMGT-gapped templates.

Design Principles:
- Use BioPython PairwiseAligner for alignment (explicit tool choice)
- Per-gene template when available, per-segment consensus as fallback
- Translate to amino acid before alignment for better accuracy
- Map gaps back to nucleotide positions at codon boundaries
- Graceful failure: log warning and return ungapped on error

IMGT Numbering:
- V genes: positions 1-104 (CDR1: 27-38, CDR2: 56-65)
- J genes: positions 1-13 (minimal structure)
- D genes: no gapping (ungapped only)

Gap Character: "." (period) per IMGT convention
"""

import logging
from pathlib import Path
from typing import Optional, Dict, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

logger = logging.getLogger(__name__)


# IMGT gap character
GAP_CHAR = "."

# Standard genetic code for translation
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


class GapperService:
    """
    Service for gapping ungapped germline sequences to IMGT numbering.

    Uses BioPython pairwise alignment against IMGT-gapped templates.
    Supports V and J segments only (D segments remain ungapped).

    Examples
    --------
    >>> gapper = GapperService(template_dir=Path("sources/imgt/human"))
    >>> gapped = gapper.gap_sequence(
    ...     sequence="CAGGTGCAGCTGGTGGAG...",
    ...     segment="V",
    ...     chain="H",
    ...     gene_name="IGHV1-69*01"
    ... )
    """

    def __init__(self, template_dir: Optional[Path] = None):
        """
        Initialize gapper with template directory.

        Parameters
        ----------
        template_dir : Path, optional
            Directory containing IMGT-gapped FASTA templates.
            If not provided, consensus templates will be built on demand.
        """
        self.template_dir = template_dir
        self._template_cache: Dict[str, Dict[str, str]] = {}
        self._consensus_cache: Dict[str, str] = {}
        self._aligner = self._create_aligner()

    def _create_aligner(self) -> PairwiseAligner:
        """
        Create BioPython PairwiseAligner with appropriate settings.

        Returns
        -------
        PairwiseAligner
            Configured aligner for amino acid sequences
        """
        aligner = PairwiseAligner()

        # Global alignment for full sequence alignment
        aligner.mode = 'global'

        # Scoring parameters optimized for antibody sequences
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5

        # Allow end gaps without penalty (sequences may be partial)
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0

        return aligner

    def gap_sequence(
        self,
        sequence: str,
        segment: str,
        chain: str,
        gene_name: Optional[str] = None,
        species: str = "human"
    ) -> Optional[str]:
        """
        Gap an ungapped nucleotide sequence to IMGT numbering.

        Parameters
        ----------
        sequence : str
            Ungapped nucleotide sequence (ACGT only)
        segment : str
            Segment type: "V" or "J" (D segments return None)
        chain : str
            Chain type: "H", "K", or "L"
        gene_name : str, optional
            Gene name for per-gene template lookup
        species : str
            Species for template lookup (default: "human")

        Returns
        -------
        str or None
            IMGT-gapped nucleotide sequence, or None if gapping fails
            or segment is "D"
        """
        # D segments are not gapped
        if segment.upper() == "D":
            logger.debug(f"D segment '{gene_name}' - no gapping required")
            return None

        # Validate segment
        if segment.upper() not in ["V", "J"]:
            logger.warning(f"Unknown segment '{segment}' - skipping gapping")
            return None

        try:
            # Get template (per-gene or consensus)
            template = self._get_template(
                segment=segment,
                chain=chain,
                gene_name=gene_name,
                species=species
            )

            if template is None:
                logger.warning(
                    f"No template found for {gene_name or f'{chain}{segment}'} - "
                    f"returning ungapped"
                )
                return None

            # Translate to amino acid for alignment
            query_aa = self._translate(sequence)
            template_aa = self._translate(template.replace(GAP_CHAR, ""))

            if query_aa is None or template_aa is None:
                logger.warning(
                    f"Translation failed for {gene_name} - returning ungapped"
                )
                return None

            # Perform alignment
            gapped_sequence = self._align_and_gap(
                sequence=sequence,
                template=template,
                query_aa=query_aa,
                template_aa=template_aa
            )

            if gapped_sequence:
                logger.debug(f"Successfully gapped {gene_name}")

            return gapped_sequence

        except Exception as e:
            logger.warning(
                f"Failed to gap sequence {gene_name}: {e}"
            )
            return None

    def _get_template(
        self,
        segment: str,
        chain: str,
        gene_name: Optional[str] = None,
        species: str = "human"
    ) -> Optional[str]:
        """
        Get IMGT-gapped template for alignment.

        Strategy:
        1. Look for per-gene template (exact match)
        2. Fall back to per-segment consensus template

        Parameters
        ----------
        segment : str
            Segment type
        chain : str
            Chain type
        gene_name : str, optional
            Gene name for exact match
        species : str
            Species name

        Returns
        -------
        str or None
            IMGT-gapped template sequence
        """
        cache_key = f"{species}_{chain}{segment}"

        # Load templates if not cached
        if cache_key not in self._template_cache:
            self._load_templates(species, segment, chain)

        templates = self._template_cache.get(cache_key, {})

        # Try per-gene template first
        if gene_name and gene_name in templates:
            logger.debug(f"Using per-gene template for {gene_name}")
            return templates[gene_name]

        # Fall back to consensus
        consensus_key = f"{species}_{chain}{segment}_consensus"
        if consensus_key not in self._consensus_cache:
            self._build_consensus(species, segment, chain, templates)

        consensus = self._consensus_cache.get(consensus_key)

        if consensus:
            logger.debug(
                f"Using consensus template for {gene_name or cache_key}"
            )

        return consensus

    def _load_templates(
        self,
        species: str,
        segment: str,
        chain: str
    ) -> None:
        """
        Load IMGT-gapped templates from FASTA file.

        Supports multiple naming conventions:
        - {template_dir}/{species}/IG{chain}{segment}_gapped.fasta
        - {template_dir}/{species}/IG{chain}{segment}.fasta
        - {template_dir}/IG{chain}{segment}_gapped.fasta (if template_dir is species-specific)

        Parameters
        ----------
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type
        """
        cache_key = f"{species}_{chain}{segment}"
        self._template_cache[cache_key] = {}

        if self.template_dir is None:
            logger.debug("No template directory configured")
            return

        segment_name = f"IG{chain}{segment}"
        
        # Try multiple path patterns
        candidate_paths = [
            # Pattern 1: {template_dir}/{species}/{segment}_gapped.fasta
            self.template_dir / species / f"{segment_name}_gapped.fasta",
            # Pattern 2: {template_dir}/{species}/{segment}.fasta
            self.template_dir / species / f"{segment_name}.fasta",
            # Pattern 3: {template_dir}/{segment}_gapped.fasta (template_dir is species-specific)
            self.template_dir / f"{segment_name}_gapped.fasta",
            # Pattern 4: {template_dir}/{segment}.fasta
            self.template_dir / f"{segment_name}.fasta",
        ]
        
        fasta_path = None
        for path in candidate_paths:
            if path.exists():
                fasta_path = path
                break
        
        if fasta_path is None:
            logger.debug(f"No template file found for {segment_name} in {self.template_dir}")
            return

        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                sequence = str(record.seq).upper()

                # Only cache gapped sequences (containing . or -)
                if GAP_CHAR in sequence or "-" in sequence:
                    # Extract gene name from IMGT header format: accession|gene_name|...
                    header_parts = record.id.split("|")
                    if len(header_parts) >= 2:
                        gene_name = header_parts[1]  # Second field is gene name
                    else:
                        gene_name = header_parts[0]  # Fallback to first field
                    
                    # Normalize to use "." as gap character
                    sequence = sequence.replace("-", GAP_CHAR)
                    self._template_cache[cache_key][gene_name] = sequence

            logger.info(
                f"Loaded {len(self._template_cache[cache_key])} "
                f"gapped templates for {cache_key} from {fasta_path}"
            )

        except Exception as e:
            logger.error(f"Failed to load templates from {fasta_path}: {e}")

    def _build_consensus(
        self,
        species: str,
        segment: str,
        chain: str,
        templates: Dict[str, str]
    ) -> None:
        """
        Build consensus template from available gapped sequences.

        Uses the most common character at each position.

        Parameters
        ----------
        species : str
            Species name
        segment : str
            Segment type
        chain : str
            Chain type
        templates : Dict[str, str]
            Available gapped templates
        """
        consensus_key = f"{species}_{chain}{segment}_consensus"

        if not templates:
            self._consensus_cache[consensus_key] = None
            return

        # Find maximum length
        max_len = max(len(seq) for seq in templates.values())

        # Build consensus by position
        consensus = []
        for pos in range(max_len):
            chars = []
            for seq in templates.values():
                if pos < len(seq):
                    chars.append(seq[pos])

            if chars:
                # Most common character at this position
                consensus.append(max(set(chars), key=chars.count))

        self._consensus_cache[consensus_key] = "".join(consensus)

        logger.debug(
            f"Built consensus template for {consensus_key} "
            f"from {len(templates)} sequences"
        )

    def _translate(self, nucleotide_seq: str) -> Optional[str]:
        """
        Translate nucleotide sequence to amino acids.

        Handles sequences not divisible by 3 by truncating.

        Parameters
        ----------
        nucleotide_seq : str
            Nucleotide sequence (gaps removed)

        Returns
        -------
        str or None
            Amino acid sequence, or None if translation fails
        """
        # Remove any gap characters
        clean_seq = nucleotide_seq.replace(GAP_CHAR, "").replace("-", "")

        # Truncate to codon boundary
        truncated = clean_seq[:len(clean_seq) - (len(clean_seq) % 3)]

        if len(truncated) < 3:
            return None

        amino_acids = []
        for i in range(0, len(truncated), 3):
            codon = truncated[i:i+3].upper()
            aa = CODON_TABLE.get(codon, 'X')  # X for unknown
            amino_acids.append(aa)

        return "".join(amino_acids)

    def _align_and_gap(
        self,
        sequence: str,
        template: str,
        query_aa: str,
        template_aa: str
    ) -> Optional[str]:
        """
        Align amino acid sequences and map gaps to nucleotide.

        Parameters
        ----------
        sequence : str
            Ungapped nucleotide sequence
        template : str
            IMGT-gapped nucleotide template
        query_aa : str
            Query amino acid sequence
        template_aa : str
            Template amino acid sequence (ungapped)

        Returns
        -------
        str or None
            Gapped nucleotide sequence
        """
        try:
            # Perform amino acid alignment
            alignments = list(self._aligner.align(template_aa, query_aa))

            if not alignments:
                logger.debug("No alignments found")
                return None

            # Use best alignment
            best = alignments[0]

            # Extract gap positions from alignment
            gap_positions = self._extract_gap_positions(template)

            # Apply gaps to query sequence at codon boundaries
            gapped = self._apply_gaps(sequence, gap_positions)

            return gapped

        except Exception as e:
            logger.debug(f"Alignment failed: {e}")
            return None

    def _extract_gap_positions(self, template: str) -> List[int]:
        """
        Extract positions of gaps in template sequence.

        Parameters
        ----------
        template : str
            IMGT-gapped template

        Returns
        -------
        List[int]
            Positions where gaps occur (0-indexed)
        """
        positions = []
        for i, char in enumerate(template):
            if char == GAP_CHAR or char == "-":
                positions.append(i)
        return positions

    def _apply_gaps(
        self,
        sequence: str,
        gap_positions: List[int]
    ) -> str:
        """
        Apply gaps to nucleotide sequence.

        The gap_positions are already in nucleotide coordinates from the template.
        We need to map these template positions to query positions, accounting
        for the fact that gap positions in the template include the gaps themselves.

        Parameters
        ----------
        sequence : str
            Ungapped nucleotide sequence
        gap_positions : List[int]
            Positions where gaps exist in the template (nucleotide coordinates)

        Returns
        -------
        str
            Gapped nucleotide sequence
        """
        result = list(sequence)
        
        # The gap_positions are positions in the GAPPED template where '.' appears.
        # To apply to the ungapped query, we need to:
        # 1. Convert template gapped positions to ungapped positions
        # 2. Insert gaps at those positions in the query
        
        # Sort and insert from end to preserve indices
        for i, pos in enumerate(sorted(gap_positions, reverse=True)):
            # The position in the ungapped sequence is:
            # gapped_pos - (number of gaps before this position)
            # Count gaps before this position
            gaps_before = sum(1 for p in gap_positions if p < pos)
            ungapped_pos = pos - gaps_before
            
            # Only insert if within bounds
            if ungapped_pos <= len(result):
                result.insert(ungapped_pos, GAP_CHAR)

        return "".join(result)


def gap_sequences_batch(
    sequences: List[Tuple[str, str, str, str]],
    template_dir: Path,
    species: str = "human"
) -> Dict[str, Optional[str]]:
    """
    Gap multiple sequences in batch.

    Parameters
    ----------
    sequences : List[Tuple[str, str, str, str]]
        List of (gene_name, sequence, segment, chain) tuples
    template_dir : Path
        Directory containing IMGT templates
    species : str
        Species name

    Returns
    -------
    Dict[str, Optional[str]]
        Mapping of gene_name to gapped sequence (None if failed)
    """
    gapper = GapperService(template_dir=template_dir)
    results = {}

    for gene_name, sequence, segment, chain in sequences:
        results[gene_name] = gapper.gap_sequence(
            sequence=sequence,
            segment=segment,
            chain=chain,
            gene_name=gene_name,
            species=species
        )

    return results
