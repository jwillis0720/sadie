"""Parity test validating G3 and Germlines produce identical AIRR output.

This test ensures that both backends produce equivalent results when built
from the same reference.g3.yml configuration. Any difference (except source
tracking columns) indicates a regression or bug.
"""
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from sadie.airr import Airr


# Columns that are expected to differ between backends (source tracking)
EXCLUDED_COLUMNS = frozenset([
    "v_call_source",
    "d_call_source",
    "j_call_source",
    "c_call_source",
])


def values_equal(v1: Any, v2: Any) -> bool:
    """Compare values treating NaN == NaN as True.
    
    Args:
        v1: First value to compare.
        v2: Second value to compare.
        
    Returns:
        True if values are equal (NaN == NaN is True).
    """
    if pd.isna(v1) and pd.isna(v2):
        return True
    if pd.isna(v1) or pd.isna(v2):
        return False
    return v1 == v2


def compare_airr_outputs(
    g3_df: pd.DataFrame,
    germlines_df: pd.DataFrame,
    fasta_name: str,
) -> None:
    """Compare two AIRR DataFrames, fail fast on first mismatch.
    
    Args:
        g3_df: DataFrame from G3 backend annotation.
        germlines_df: DataFrame from Germlines backend annotation.
        fasta_name: Name of the FASTA file being tested (for error messages).
        
    Raises:
        pytest.fail: On first detected mismatch with detailed report.
    """
    # 1. Check column presence (excluding source columns)
    g3_cols = set(g3_df.columns) - EXCLUDED_COLUMNS
    germlines_cols = set(germlines_df.columns) - EXCLUDED_COLUMNS

    if g3_cols != germlines_cols:
        g3_only = g3_cols - germlines_cols
        germlines_only = germlines_cols - g3_cols
        pytest.fail(
            f"\nCOLUMN MISMATCH DETECTED\n"
            f"\nSequence file: {fasta_name}\n"
            f"G3-only columns: {sorted(g3_only)}\n"
            f"Germlines-only columns: {sorted(germlines_only)}"
        )

    # 2. Sort and align rows by sequence_id
    g3_df = g3_df.sort_values("sequence_id").reset_index(drop=True)
    germlines_df = germlines_df.sort_values("sequence_id").reset_index(drop=True)

    # 3. Check row count
    if len(g3_df) != len(germlines_df):
        pytest.fail(
            f"\nROW COUNT MISMATCH\n"
            f"\nSequence file: {fasta_name}\n"
            f"G3 rows: {len(g3_df)}\n"
            f"Germlines rows: {len(germlines_df)}"
        )

    # 4. Compare cell by cell, excluding source columns
    compare_cols = sorted(g3_cols)

    for row_idx in range(len(g3_df)):
        for col in compare_cols:
            g3_val = g3_df.loc[row_idx, col]
            germlines_val = germlines_df.loc[row_idx, col]

            if not values_equal(g3_val, germlines_val):
                seq_id = g3_df.loc[row_idx, "sequence_id"]
                pytest.fail(
                    f"\nPARITY MISMATCH DETECTED\n"
                    f"\nSequence file: {fasta_name}\n"
                    f"Row: {row_idx}\n"
                    f"Column: {col}\n"
                    f"Sequence ID: {seq_id}\n"
                    f"\nG3 value:        {g3_val!r}\n"
                    f"Germlines value: {germlines_val!r}"
                )


# Test data: human FASTA files
FIXTURES_DIR = Path(__file__).parent.parent / "data" / "fixtures" / "fasta_inputs"

HUMAN_TEST_FASTAS = [
    FIXTURES_DIR / "PG9_H.fasta",
    FIXTURES_DIR / "PG9_L.fasta",
    FIXTURES_DIR / "OAS_subsample_1000.fasta",
    FIXTURES_DIR / "catnap_nt_heavy.fasta",
    FIXTURES_DIR / "catnap_nt_light.fasta",
]


@pytest.mark.parametrize(
    "fasta_file",
    HUMAN_TEST_FASTAS,
    ids=lambda p: p.name,
)
def test_airr_parity(
    g3_database: Path,
    germlines_database: Path,
    fasta_file: Path,
) -> None:
    """Test that G3 and Germlines backends produce identical AIRR output.
    
    This test:
    1. Runs AIRR annotation with G3-built database
    2. Runs AIRR annotation with Germlines-built database
    3. Compares all columns except source tracking columns
    4. Fails immediately on first mismatch with detailed report
    
    Args:
        g3_database: Path to G3-built database (session fixture).
        germlines_database: Path to Germlines-built database (session fixture).
        fasta_file: Path to FASTA file to annotate.
    """
    # Skip if fixture file doesn't exist
    if not fasta_file.exists():
        pytest.skip(f"Test file not found: {fasta_file}")

    # Run annotation with G3 database
    g3_api = Airr("human", database=g3_database)
    g3_result = g3_api.run_fasta(fasta_file)

    # Run annotation with Germlines database
    germlines_api = Airr("human", database=germlines_database)
    germlines_result = germlines_api.run_fasta(fasta_file)

    # Compare outputs (AirrTable has DataFrame interface)
    compare_airr_outputs(
        pd.DataFrame(g3_result),
        pd.DataFrame(germlines_result),
        fasta_file.name,
    )
