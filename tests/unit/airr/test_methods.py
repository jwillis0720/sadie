import logging

import pandas as pd
import pytest
from numpy import nan

from sadie.airr import AirrTable
from sadie.airr.methods import codon_table, find_best_codon, get_igl_aa

LOGGER = logging.getLogger("AirrMethod")


@pytest.mark.parametrize(
    "codon, aa, expected",
    [
        ("", "M", "ATG"),
        ("ATG", "M", "ATG"),
        ("FAKE", "M", "ATG"),
    ]
    + [("".join(codons[::-1]), aa, codons[0]) for aa, codons in codon_table.items()],
)
def test_find_best_codon(codon: str, aa: str, expected: str) -> None:
    assert find_best_codon(codon, aa) == expected


def test_get_igl_aa(fixture_setup, caplog) -> None:
    df = pd.read_feather(fixture_setup.get_bum_igl_assignment())
    table = AirrTable(df)  # init and verify
    row = table.iloc[0]
    # germline does not exist
    _row = row.copy()
    _row.v_germline_alignment_aa = nan
    assert get_igl_aa(_row) is nan
    # mature != germline
    _row = row.copy()
    _row.germline_alignment_aa = _row.germline_alignment_aa[:10]
    assert get_igl_aa(_row) is nan
    # insertion in germline
    _row = row.copy()
    _row.germline_alignment_aa = _row.germline_alignment_aa[:1] + "-" * 10 + _row.germline_alignment_aa[11:]
    with caplog.at_level(logging.WARNING):
        get_igl_aa(_row)
    assert "has a insertion in it at CDR3" in caplog.text
    # X or * in aa position for both germline and mature
    _row = row.copy()
    _row.germline_alignment_aa = _row.germline_alignment_aa[:1] + "*" * 10 + _row.germline_alignment_aa[11:]
    _row.sequence_alignment_aa = _row.sequence_alignment_aa[:1] + "*" * 10 + _row.germline_alignment_aa[11:]
    with pytest.raises(ValueError):
        get_igl_aa(_row)
    # full igl has a stop codon
    _row = row.copy()
    _row.v_germline_alignment_aa = _row.v_germline_alignment_aa[:1] + "*" * 10 + _row.v_germline_alignment_aa[11:]
    with pytest.raises(ValueError):
        get_igl_aa(_row)
