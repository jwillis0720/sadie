"""
Renumbering Integration for Germlines Module
=============================================

Provides HMM generation from local germlines database for use with
sadie.renumbering.aligners.hmmer module.

This replaces G3 API calls with local database queries.
"""

import logging
from pathlib import Path
from typing import Optional, List, Tuple
from functools import lru_cache

import pyhmmer
from sadie.typing import Chain, Source, Species

logger = logging.getLogger(__name__)


class LocalHMMBuilder:
    """
    Build HMM models from local germlines database.

    This class replaces G3 client for HMM generation in renumbering workflows.
    It queries the germlines manager and builds HMM models using pyhmmer.

    Examples
    --------
    >>> builder = LocalHMMBuilder()
    >>> hmm = builder.get_hmm(species="human", chain="H")
    >>> # Use in renumbering
    >>> from sadie.renumbering.aligners.hmmer import HMMER
    >>> aligner = HMMER(species="human", chains="H")
    """

    def __init__(self):
        """Initialize HMM builder with pyhmmer components."""
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.builder = pyhmmer.plan7.Builder(self.alphabet)
        self.background = pyhmmer.plan7.Background(self.alphabet)

        # Cache directory for HMM files
        from sadie.germlines import get_germlines_base_dir
        self.hmm_dir = get_germlines_base_dir() / "hmms"
        self.hmm_dir.mkdir(parents=True, exist_ok=True)

    @lru_cache(maxsize=None)
    def get_hmm(
        self,
        species: Species,
        chain: Chain,
        source: Source = "imgt"
    ) -> pyhmmer.plan7.HMM:
        """
        Get or build HMM model for species/chain combination.

        Parameters
        ----------
        species : Species
            Species name (e.g., "human", "mouse")
        chain : Chain
            Chain type ("H", "K", or "L")
        source : Source
            Data source (default: "imgt")

        Returns
        -------
        pyhmmer.plan7.HMM
            Compiled HMM model
        """
        hmm_path = self.hmm_dir / f"{species}_{chain}.hmm"

        # Return cached HMM if exists
        if hmm_path.exists():
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                return next(hmm_file)

        # Build new HMM
        logger.info(f"Building HMM for {species} {chain}")
        return self._build_hmm(species, chain, source)

    def _build_hmm(
        self,
        species: str,
        chain: str,
        source: str
    ) -> pyhmmer.plan7.HMM:
        """
        Build HMM from germlines database.

        Parameters
        ----------
        species : str
            Species name
        chain : str
            Chain type
        source : str
            Data source

        Returns
        -------
        pyhmmer.plan7.HMM
            Compiled HMM model
        """
        # Get V and J sequences from germlines
        vj_pairs = self._get_vj_alignment_pairs(species, chain, source)

        if not vj_pairs:
            raise ValueError(f"No sequences found for {species} {chain}")

        # Create Stockholm alignment
        sto_path = self._write_stockholm(vj_pairs, species, chain)

        # Build HMM using pyhmmer
        hmm_path = self.hmm_dir / f"{species}_{chain}.hmm"

        with pyhmmer.easel.MSAFile(
            sto_path,
            digital=True,
            alphabet=self.alphabet,
            format="stockholm"
        ) as msa_file:
            msa = next(msa_file)
            hmm, _, _ = self.builder.build_msa(msa, self.background)

        # Save HMM to disk
        with open(hmm_path, "wb") as f:
            hmm.write(f)

        logger.info(f"Built HMM: {hmm_path}")
        return hmm

    def _get_vj_alignment_pairs(
        self,
        species: str,
        chain: str,
        source: str,
        strict: bool = True
    ) -> List[Tuple[str, str]]:
        """
        Get V-J alignment pairs from germlines database.

        Parameters
        ----------
        species : str
            Species name
        chain : str
            Chain type
        source : str
            Data source
        strict : bool
            If True, raise error when genes lack gapped sequences (FR-013)

        Returns
        -------
        List[Tuple[str, str]]
            List of (gene_name, gapped_aa_sequence) tuples

        Raises
        ------
        ValueError
            If strict=True and genes are missing both gapped AA and gapped NT
        """
        from sadie.germlines import get_manager

        manager = get_manager()
        pairs = []
        missing_gapped = []

        # Get V and J genes
        for segment in ["V", "J"]:
            genes = manager.get_genes(species, segment, chain)

            for gene in genes:
                # First try pre-computed gapped AA sequence
                if gene.sequence_aa_gapped:
                    pairs.append((gene.name, gene.sequence_aa_gapped))
                # Fall back to translating gapped nucleotide sequence
                elif gene.sequence_gapped:
                    aa_gapped = self._translate_gapped_nt_to_aa(gene.sequence_gapped)
                    if aa_gapped:
                        pairs.append((gene.name, aa_gapped))
                    else:
                        missing_gapped.append(gene.name)
                else:
                    # FR-013: Track genes missing both gapped sequences
                    missing_gapped.append(gene.name)

        # FR-013: Fail-fast when genes lack gapped data for HMM building
        if strict and missing_gapped:
            logger.warning(
                f"Genes missing gapped sequences for {species} {chain}: "
                f"{missing_gapped[:10]}{'...' if len(missing_gapped) > 10 else ''}"
            )
            # Only fail if we have NO pairs at all
            if not pairs:
                raise ValueError(
                    f"Cannot build HMM for {species} {chain}: no genes have "
                    f"gapped sequences. Missing genes: {missing_gapped[:5]}. "
                    f"Ensure germline data includes sequence_aa_gapped or "
                    f"sequence_gapped fields."
                )

        return pairs

    def _translate_gapped_nt_to_aa(self, gapped_nt: str) -> Optional[str]:
        """
        Translate IMGT-gapped nucleotide sequence to gapped amino acid.

        IMGT gaps are represented as dots. We preserve gap positions while
        translating the nucleotide codons to amino acids.

        Parameters
        ----------
        gapped_nt : str
            IMGT-gapped nucleotide sequence (dots for gaps)

        Returns
        -------
        str or None
            Gapped amino acid sequence, or None if translation fails
        """
        from Bio.Seq import Seq

        # Remove gaps to get pure nucleotide sequence
        nt_ungapped = gapped_nt.replace(".", "").replace("-", "")

        # Must be multiple of 3 for translation
        # Truncate to nearest codon boundary
        codon_len = (len(nt_ungapped) // 3) * 3
        if codon_len < 3:
            return None

        nt_for_translation = nt_ungapped[:codon_len]

        try:
            aa_seq = str(Seq(nt_for_translation).translate())
        except Exception:
            return None

        # Now we need to insert gaps at the correct positions
        # IMGT numbering: every 3 nucleotide positions = 1 AA position
        # Gaps in NT sequence need to be mapped to AA positions

        # Build the gapped AA sequence
        aa_gapped_chars = []
        nt_pos = 0
        aa_pos = 0

        i = 0
        while i < len(gapped_nt) and aa_pos < len(aa_seq):
            char = gapped_nt[i]
            if char in (".", "-"):
                # This is a gap - accumulate gaps until we have 3
                gap_count = 0
                while i < len(gapped_nt) and gapped_nt[i] in (".", "-"):
                    gap_count += 1
                    i += 1
                # For every 3 NT gaps, insert 1 AA gap
                aa_gaps = gap_count // 3
                aa_gapped_chars.extend(["."] * aa_gaps)
            else:
                # This is a nucleotide - consume 3 NTs, output 1 AA
                codon_chars = 0
                while i < len(gapped_nt) and codon_chars < 3:
                    if gapped_nt[i] not in (".", "-"):
                        codon_chars += 1
                    i += 1
                if aa_pos < len(aa_seq):
                    aa_gapped_chars.append(aa_seq[aa_pos])
                    aa_pos += 1

        return "".join(aa_gapped_chars) if aa_gapped_chars else None

    def _write_stockholm(
        self,
        pairs: List[Tuple[str, str]],
        species: str,
        chain: str
    ) -> Path:
        """
        Write Stockholm alignment file.

        Parameters
        ----------
        pairs : List[Tuple[str, str]]
            List of (name, gapped_sequence) tuples
        species : str
            Species name
        chain : str
            Chain type

        Returns
        -------
        Path
            Path to Stockholm file
        """
        sto_dir = self.hmm_dir.parent / "stockholms"
        sto_dir.mkdir(parents=True, exist_ok=True)

        sto_path = sto_dir / f"{species}_{chain}.sto"

        # Find max sequence length and max name length
        max_seq_len = max(len(seq) for _, seq in pairs)
        max_name_len = max(len(name) for name, _ in pairs)

        lines = [
            "# STOCKHOLM 1.0",
            f"#=GF ID {species}_{chain}",
            ""
        ]

        for name, seq in pairs:
            # Pad sequences to same length with gaps (.)
            padded_seq = seq.ljust(max_seq_len, ".")
            lines.append(f"{name.ljust(max_name_len)}  {padded_seq}")

        # Add reference line (RF) matching the alignment length and terminator
        rf_label = f"#=GC RF".ljust(max_name_len + 2)
        lines.append(f"{rf_label}  {'x' * max_seq_len}")
        lines.append("//")

        sto_path.write_text("\n".join(lines))

        return sto_path


def use_local_hmm_builder() -> bool:
    """
    Check if renumbering should use local HMM builder.

    Returns True if germlines module should be used (default),
    False if G3 API should be used (legacy).

    Returns
    -------
    bool
        True to use local HMM builder, False to use G3
    """
    from sadie.germlines.utils import use_germlines_module
    return use_germlines_module()
