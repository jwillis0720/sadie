import logging

import pandas as pd
import pytest
from numpy import nan

from sadie.airr import AirrTable
from sadie.airr.methods import (
    codon_table,
    find_best_codon,
    get_igl_aa,
    get_igl_nt,
    run_five_prime_buffer,
    run_igl_assignment,
    run_mutational_analysis,
    run_termini_buffers,
    run_three_prime_buffer,
)
from sadie.reference.reference import References

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
    _row.germline_alignment_aa = _row.germline_alignment_aa[:1] + "*" * 100 + _row.germline_alignment_aa[101:]
    _row.sequence_alignment_aa = _row.germline_alignment_aa
    with pytest.raises(ValueError):
        get_igl_aa(_row)
    # full igl has a stop codon - should return nan instead of raising error
    _row = row.copy()
    _row.v_germline_alignment_aa = _row.v_germline_alignment_aa[:1] + "*" * 10 + _row.v_germline_alignment_aa[11:]
    # Stop codons now return nan instead of raising ValueError (more robust behavior)
    assert get_igl_aa(_row) is nan


def test_get_igl_nt(fixture_setup, caplog):
    df = pd.read_feather(fixture_setup.get_bum_igl_assignment())
    table = AirrTable(df)  # init and verify
    row = table.iloc[0]
    # germline does not exist
    _row = row.copy()
    _row.germline_alignment_aa = nan
    assert get_igl_nt(_row) is nan
    # seq alignment does not exist
    _row = row.copy()
    _row.sequence_alignment_aa = nan
    assert get_igl_nt(_row) is nan
    # seq alignment is shorter/longer than germline
    _row = row.copy()
    _row.germline_alignment_aa = _row.germline_alignment_aa[:-1]
    with pytest.raises(ValueError):
        get_igl_nt(_row)
    _row = row.copy()
    _row.germline_alignment_aa = _row.germline_alignment_aa[:1] + "-" + _row.germline_alignment_aa[2:]
    _row.sequence_alignment_aa = _row.germline_alignment_aa
    with pytest.raises(ValueError):
        get_igl_nt(_row)
    # germline has a insertion at the end & germline does not equal complete vdj
    _row = row.copy()
    _row["iGL_aa"] = get_igl_aa(_row)
    _row.germline_alignment_aa = _row.germline_alignment_aa.replace("-", "A")
    _row.sequence_alignment_aa = _row.germline_alignment_aa[:]
    _row.germline_alignment_aa = _row.germline_alignment_aa[:-1] + "-"
    with caplog.at_level(logging.WARNING):
        with pytest.warns(UserWarning):
            get_igl_nt(_row)


def test_run_igl_assignment() -> None:
    with pytest.raises(TypeError):
        run_igl_assignment(None)


def test_run_mutational_analysis(fixture_setup, caplog) -> None:
    with pytest.raises(TypeError):
        run_mutational_analysis(airrtable=None, scheme="imgt")


def test_run_five_prime_buffer(fixture_setup, caplog) -> None:
    df = pd.read_feather(fixture_setup.get_bum_igl_assignment())
    table = AirrTable(df)  # init and verify
    run_five_prime_buffer(table, references=References())
    with pytest.raises(TypeError):
        run_five_prime_buffer(None)
    _table = table.copy()
    table.at[0, "reference_name"] = "test"
    table.at[1, "reference_name"] = "test"
    with pytest.raises(ValueError):
        run_five_prime_buffer(table)


def test_run_three_prime_buffer(fixture_setup, caplog) -> None:
    df = pd.read_feather(fixture_setup.get_bum_igl_assignment())
    table = AirrTable(df)  # init and verify
    run_three_prime_buffer(table, references=References())
    with pytest.raises(TypeError):
        run_three_prime_buffer(None)
    _table = table.copy()
    table.at[0, "reference_name"] = "test"
    table.at[1, "reference_name"] = "test"
    with pytest.raises(ValueError):
        run_three_prime_buffer(table)


def test_run_termini_buffers(fixture_setup, caplog) -> None:
    with pytest.raises(TypeError):
        run_termini_buffers(None)
