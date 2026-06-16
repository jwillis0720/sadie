"""Regression guard for MOTIF_LOOKUP against SADIE's IgBLAST aux data."""

import re
import warnings
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from sadie.reference.settings import MOTIF_LOOKUP

J_FASTA_DIR = Path(__file__).parents[3] / "src/sadie/airr/data/germlines/Ig/blastdb"
AUX_DIR = Path(__file__).parents[3] / "src/sadie/airr/data/germlines/aux_db/imgt"
SHIPPED_SPECIES = ["human", "mouse", "rat", "macaque", "dog", "rabbit"]
MARKER_TO_LOCUS = {"JH": "IGHJ", "JK": "IGKJ", "JL": "IGLJ"}


def _ignore_set(entry):
    ig = entry.get("ignore")
    if not ig:
        return set()
    if isinstance(ig, str):
        return {ig}
    return {x for x in ig if x}


def _translate_frame(seq, reading_frame):
    s = str(seq)
    sub = s[reading_frame:]
    sub = sub[: len(sub) - len(sub) % 3]
    return str(Seq(sub).translate())


def _aux_rows(path):
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 5:
            continue
        gene, reading_frame, marker, fr4_start, _left_over = cols[:5]
        yield gene, int(reading_frame), marker, int(fr4_start)


@pytest.mark.parametrize("species", SHIPPED_SPECIES)
def test_motif_matches_igblast_aux_fr4_start(species):
    fasta = J_FASTA_DIR / species / f"{species}_J.fasta"
    aux = AUX_DIR / f"{species}_gl.aux"
    assert fasta.exists(), f"missing shipped J data for {species}"
    assert aux.exists(), f"missing IgBLAST aux data for {species}"

    records = {rec.id: rec for rec in SeqIO.parse(fasta, "fasta")}
    motifs = MOTIF_LOOKUP[species]
    ignore = _ignore_set(motifs)

    checked = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Biopython partial-codon warning on J segments
        for gene, reading_frame, marker, fr4_start in _aux_rows(aux):
            rec = records.get(gene)
            if rec is None:
                continue
            locus = MARKER_TO_LOCUS.get(marker)
            pattern = motifs.get(locus)
            if pattern is None:
                continue
            if gene in ignore or gene.split("*")[0] in ignore:
                continue

            aa = _translate_frame(rec.seq, reading_frame)
            motif_starts = [match.start() * 3 + reading_frame - 1 for match in re.finditer(pattern, aa)]
            assert motif_starts == [
                fr4_start
            ], f"{species} {gene}: motif {pattern!r} starts at {motif_starts}, aux FR4 starts at {fr4_start}"
            checked += 1

    assert checked > 0, f"no {species} aux alleles exercised a motif"
