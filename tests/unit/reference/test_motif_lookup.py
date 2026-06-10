"""Empirical validation of ``sadie.reference.settings.MOTIF_LOOKUP``.

Background (GitHub issue #267)
------------------------------
A user asked three things about the J-gene motifs in ``MOTIF_LOOKUP``:

1. Can these motifs identify the FWR4 start on germline J alleles?
2. Where do the non-canonical variations (e.g. ``WGQG`` platypus, ``W[GD].G``
   horse, ``[WL]G[TQK][VG]`` alpaca, human IGKJ ``FG``) come from?
3. How reliable are they for *de-novo* J alleles not in IMGT/OGRDB?

``MOTIF_LOOKUP`` is a curated reference of the conserved J framework-4 anchor
(the IMGT J-TRP 118 / J-PHE 118 invariant -> ``W-G-x-G`` for heavy, ``F-G-x-G``
for light), empirically fitted to the IMGT/OGRDB J alleles per organism. It is
**not** wired into SADIE's current FWR4 pipeline (that is driven by the IgBLAST
``.aux`` files + numbering schemes / G3) -- see the docstring above the constant.

This module pins that empirical reality: for every *functional* Ig germline J
allele SADIE ships, the species/locus motif must locate the FR4 anchor at the
C-terminus of the translated J segment. It guards ``MOTIF_LOOKUP`` against silent
drift and documents the known, deliberate exceptions (pseudogenes, ``ignore``
entries, and one genuinely divergent allele).
"""
import re
from pathlib import Path
from typing import List, Optional, Tuple

import pytest
from Bio.Seq import Seq

from sadie.airr import __file__ as _sadie_airr_file
from sadie.reference.settings import MOTIF_LOOKUP

# Install-location-independent path to the shipped germline blast databases,
# resolved the same way as tests/unit/airr/test_igblast.py.
BLASTDB_DIR = Path(_sadie_airr_file).parent / "data" / "germlines" / "Ig" / "blastdb"

IG_J_LOCI = ("IGHJ", "IGKJ", "IGLJ")

# Functional alleles whose motif legitimately does NOT match. Each is a
# deliberate, documented exception rather than a motif defect:
#   - ("macaque", "IGLJ7*02"): divergent allele translating to ...L-G-R-G-T...
#     (no Phe/Cys before the G-x-G-T), so the macaque IGL motif F[GC].GT cannot
#     anchor it. A real limitation, not drift.
KNOWN_MOTIF_MISSES = {
    ("macaque", "IGLJ7*02"),
}


def _read_fasta(path: Path) -> List[Tuple[str, str]]:
    """Return ``[(allele_name, nt_sequence), ...]`` from a germline J FASTA."""
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    chunks: List[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:]
            chunks = []
        elif line:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def _translate(nt: str, frame: int) -> str:
    trimmed = nt[frame:]
    trimmed = trimmed[: len(trimmed) - (len(trimmed) % 3)]
    return str(Seq(trimmed).translate())


def locates_fr4_anchor(nt: str, motif: str) -> bool:
    """True if ``motif`` finds the conserved FR4 anchor in any reading frame.

    The germline J FASTA headers carry no reading frame, so all three frames are
    tried. A hit only counts as the FR4 anchor when the match sits at the
    C-terminus -- there is no stop codon between the anchor and the end of the
    segment, and the anchor-to-end distance is within the J framework-4 length
    range (W/F-G-x-G ... V-T-V-S-S is ~10-13 residues). This rejects spurious
    internal matches in out-of-frame translations.
    """
    for frame in (0, 1, 2):
        aa = _translate(nt, frame)
        for match in re.finditer(motif, aa):
            tail = aa[match.start() :]
            if "*" in tail:
                continue
            if 5 <= len(tail) <= 16:
                return True
    return False


def _ignore_set(species: str) -> set:
    """Normalize the heterogeneous ``ignore`` field into a set of gene/allele ids."""
    raw = MOTIF_LOOKUP[species].get("ignore", "")
    if isinstance(raw, str):
        raw = [raw]
    return {x for x in raw if x}


def _is_pseudogene(allele: str) -> bool:
    # IMGT pseudogene gene names carry a 'P' after the locus+number, e.g. IGHJ3P, IGLJ2P.
    return re.match(r"IG[HKL]J\d+P", allele) is not None


def _is_functional(species: str, allele: str) -> bool:
    if _is_pseudogene(allele):
        return False
    ignore = _ignore_set(species)
    gene = allele.split("*", 1)[0]
    return allele not in ignore and gene not in ignore


def _testable_species_loci() -> List[Tuple[str, str]]:
    """(species, locus) pairs that have BOTH a motif and shipped germline data."""
    pairs: List[Tuple[str, str]] = []
    if not BLASTDB_DIR.exists():
        return pairs
    for species_dir in sorted(BLASTDB_DIR.iterdir()):
        species = species_dir.name
        if species not in MOTIF_LOOKUP:
            continue
        if not (species_dir / f"{species}_J.fasta").exists():
            continue
        for locus in IG_J_LOCI:
            if locus in MOTIF_LOOKUP[species]:
                pairs.append((species, locus))
    return pairs


@pytest.mark.parametrize("species,locus", _testable_species_loci())
def test_motif_locates_fr4_on_functional_alleles(species: str, locus: str) -> None:
    """Every functional germline J allele's FR4 start is found by its motif.

    This is the empirical answer to issue #267 Q1/Q3: on the curated in-database
    alleles the motifs are essentially 100% accurate at locating FR4; the only
    misses are pseudogenes, deliberately ``ignore``-d alleles, or the documented
    divergent allele(s) in ``KNOWN_MOTIF_MISSES``.
    """
    motif = MOTIF_LOOKUP[species][locus]
    alleles = [
        (name, nt)
        for name, nt in _read_fasta(BLASTDB_DIR / species / f"{species}_J.fasta")
        if name[:4] == locus and _is_functional(species, name)
    ]
    if not alleles:
        pytest.skip(f"no functional {locus} alleles shipped for {species}")

    unexpected_misses = [
        name
        for name, nt in alleles
        if not locates_fr4_anchor(nt, motif) and (species, name) not in KNOWN_MOTIF_MISSES
    ]
    assert not unexpected_misses, (
        f"{species} {locus} motif {motif!r} failed to locate FR4 on functional "
        f"alleles {unexpected_misses}. Either MOTIF_LOOKUP drifted from the shipped "
        f"germline data, or a new exception must be documented in KNOWN_MOTIF_MISSES."
    )


def test_known_misses_are_still_missing() -> None:
    """Pin the documented exceptions so a fixed allele forces us to revisit the list."""
    for species, allele in KNOWN_MOTIF_MISSES:
        fasta = BLASTDB_DIR / species / f"{species}_J.fasta"
        if not fasta.exists():
            pytest.skip(f"{species} germline data not shipped")
        locus = allele[:4]
        motif = MOTIF_LOOKUP[species][locus]
        seqs = dict(_read_fasta(fasta))
        assert allele in seqs, f"{allele} no longer shipped; prune KNOWN_MOTIF_MISSES"
        assert not locates_fr4_anchor(seqs[allele], motif), (
            f"{species} {allele} now matches {motif!r}; remove it from KNOWN_MOTIF_MISSES."
        )


def test_human_igkj_motif_is_anchored_against_spurious_fg() -> None:
    """human IGKJ motif must require the full F-G-x-G anchor, not a bare 'FG'.

    The historical motif was ``FG`` (matches the FIRST Phe-Gly dipeptide anywhere
    in the translation). The old motif consumer took ``re.findall(motif)[0]``, so a
    spurious upstream ``FG`` on a de-novo allele would mis-call FR4 too early. The
    anchored ``FG.G`` (the human IGKJ FR4 consensus FGQG/FGPG/FGGG) avoids that.
    """
    motif = MOTIF_LOOKUP["human"]["IGKJ"]
    assert motif == r"FG.G", "human IGKJ motif regressed to an unanchored pattern"

    # Synthetic CDR3 tail containing an early spurious 'FG' before the true FR4.
    spurious = "CARDFGYTLDYWGFGQGTKVEIK"  # internal 'FG' at index 4, real FR4 'FGQG' later
    bare_fg = re.search(r"FG", spurious)
    anchored = re.search(r"FG.G", spurious)
    assert bare_fg is not None and anchored is not None
    assert bare_fg.start() < anchored.start()  # bare 'FG' fires early; anchored motif does not
