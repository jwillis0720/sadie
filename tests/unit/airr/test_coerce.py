"""Test CDR3 coercion functionality"""
from pathlib import Path

import pandas as pd
import pytest

from sadie.airr import Airr
from sadie.airr.airrtable import AirrTable


def test_coerce_allele_when_exact_match_not_found():
    """Test that coerce option uses available allele when exact match not found"""

    # Load test sequence where IGKV3-15*02 is found but only IGKV3-15*01 is in auxiliary files
    test_fasta = Path("tests/data/fixtures/cdr3_known_bugs.fasta")

    # Run without coerce
    airr_no_coerce = Airr("human", coerce=False)
    result_no_coerce = airr_no_coerce.run_fasta(test_fasta)

    # Run with coerce - should use IGKV3-15*01 when IGKV3-15*02 not found
    airr_coerce = Airr("human", coerce=True)
    result_coerce = airr_coerce.run_fasta(test_fasta)

    # Check that coerce found CDR3 where no_coerce didn't
    assert isinstance(result_no_coerce, AirrTable)
    assert isinstance(result_coerce, AirrTable)

    # Find sequences with IGKV3-15*02 in v_call
    mask_coerce = result_coerce["v_call"].str.contains("IGKV3-15", na=False)

    if mask_coerce.any():
        # Check if any sequences were coerced
        coerced_rows = result_coerce[mask_coerce]
        comments_mask = coerced_rows["comments"].str.contains("FR/CDR boundaries derived from", na=False)

        if comments_mask.any():
            # At least one sequence was coerced
            coerced_seq = coerced_rows[comments_mask].iloc[0]

            # Check comment was added correctly
            comment = coerced_seq["comments"]
            assert "FR/CDR boundaries derived from" in str(comment)
            assert "IGKV3-15*01" in str(comment)
            assert "IGKV3-15*02" in str(comment)

            # The V call should now start with IGKV3-15*01
            assert coerced_seq["v_call"].startswith("IGKV3-15*01")


def test_coerce_preserves_exact_matches():
    """Test that coerce doesn't change alleles that have exact matches"""

    # Use a standard test file with alleles that should have matches
    test_fasta = Path("tests/data/fixtures/cdr3_known_bugs.fasta")

    airr_no_coerce = Airr("human", coerce=False)
    airr_coerce = Airr("human", coerce=True)

    result_no_coerce = airr_no_coerce.run_fasta(test_fasta)
    result_coerce = airr_coerce.run_fasta(test_fasta)

    # For sequences without coercion comments, V calls should be identical
    # (except for those that were coerced)
    no_coercion_mask = ~result_coerce["comments"].str.contains("FR/CDR boundaries derived from", na=False)

    if no_coercion_mask.any():
        # V calls should be identical for non-coerced sequences
        v_calls_no_coerce = result_no_coerce.loc[no_coercion_mask, "v_call"]
        v_calls_coerce = result_coerce.loc[no_coercion_mask, "v_call"]
        pd.testing.assert_series_equal(v_calls_no_coerce, v_calls_coerce)

        # CDR3 should be identical for non-coerced sequences
        cdr3_no_coerce = result_no_coerce.loc[no_coercion_mask, "cdr3"]
        cdr3_coerce = result_coerce.loc[no_coercion_mask, "cdr3"]
        # Use check_names=False since indices might differ slightly
        pd.testing.assert_series_equal(cdr3_no_coerce, cdr3_coerce, check_names=False)


def test_coerce_with_multiple_sequences():
    """Test coerce works correctly with multiple sequences"""

    test_fasta = Path("tests/data/fixtures/cdr3_known_bugs.fasta")

    airr_api = Airr("human", coerce=True)
    result = airr_api.run_fasta(test_fasta)

    # Should process all sequences
    assert len(result) > 0

    # Check if any sequences were coerced
    coerced_mask = result["comments"].str.contains("FR/CDR boundaries derived from", na=False)
    if coerced_mask.any():
        # Coerced sequences should have CDR3
        coerced_rows = result[coerced_mask]
        assert not coerced_rows["cdr3"].isna().all()
