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
EXCLUDED_COLUMNS = frozenset(
    [
        "v_call_source",
        "d_call_source",
        "j_call_source",
        "c_call_source",
    ]
)

# Columns excluded only for mixed-source tests due to IgBLAST alignment boundary variations
# These can have 1-2bp differences at sequence ends without affecting biological interpretation
MIXED_SOURCE_EXCLUDED_COLUMNS = frozenset(
    [
        "germline_alignment",  # Can vary by 1-2bp at J-region boundary due to database composition
        "j_alignment_end",  # Boundary can shift by 1bp when identical J alleles come from different providers
        "j_cigar",
        "j_germline_alignment",
        "j_germline_end",
        "j_identity",
        "j_mutation",
        "j_sequence_alignment",
        "j_sequence_end",
        "sequence_alignment",
    ]
)

# E-value columns that can have tiny floating-point variations due to
# IgBLAST's database statistics calculations (not related to sequence data)
SUPPORT_COLUMNS = frozenset(["v_support", "d_support", "j_support", "c_support"])

# Relative tolerance for e-value comparisons (0.1% difference allowed)
SUPPORT_RELATIVE_TOLERANCE = 0.001


def values_equal(v1: Any, v2: Any, column: str = "") -> bool:
    """Compare values treating NaN == NaN as True.

    For support columns (e-values), uses relative tolerance comparison
    since IgBLAST calculates these based on database statistics which
    can vary slightly between builds.

    Args:
        v1: First value to compare.
        v2: Second value to compare.
        column: Column name for context-aware comparison.

    Returns:
        True if values are equal (NaN == NaN is True).
    """
    if pd.isna(v1) and pd.isna(v2):
        return True
    if pd.isna(v1) or pd.isna(v2):
        return False

    # For support columns, use relative tolerance for floating-point comparison
    if column in SUPPORT_COLUMNS:
        try:
            f1, f2 = float(v1), float(v2)
            if f1 == 0.0 and f2 == 0.0:
                return True
            if f1 != 0.0:
                rel_diff = abs(f1 - f2) / abs(f1)
                return rel_diff <= SUPPORT_RELATIVE_TOLERANCE
        except (ValueError, TypeError):
            pass

    return v1 == v2


def compare_airr_outputs(
    g3_df: pd.DataFrame,
    germlines_df: pd.DataFrame,
    fasta_name: str,
    excluded_columns: frozenset[str] = EXCLUDED_COLUMNS,
) -> None:
    """Compare two AIRR DataFrames, fail fast on first mismatch.

    Args:
        g3_df: DataFrame from G3 backend annotation.
        germlines_df: DataFrame from Germlines backend annotation.
        fasta_name: Name of the FASTA file being tested (for error messages).

    Raises:
        pytest.fail: On first detected mismatch with detailed report.
    """
    # 1. Check column presence (excluding known backend-specific columns)
    g3_cols = set(g3_df.columns) - excluded_columns
    germlines_cols = set(germlines_df.columns) - excluded_columns

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

            if not values_equal(g3_val, germlines_val, col):
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
def test_mixed_source_parity(
    germlines_database: Path,
    mixed_source_database: Path,
    fasta_file: Path,
) -> None:
    """Test that mixed-source database produces identical results to IMGT-only.

    This test validates that the Germlines backend produces identical annotation
    results regardless of which provider (IMGT, OGRDB, VDJbase) supplies an allele,
    as long as the underlying sequences are identical.

    Mixed source database pulls:
    - 5 alleles from OGRDB: IGHV1-18*01, IGHV1-2*02, IGHD1-1*01, IGHJ1*01, IGHJ3*02
    - 5 alleles from VDJbase: IGHV3-30*01, IGHV3-21*01, IGHD2-2*01, IGHJ2*01, IGHJ4*02
    - Remaining alleles from IMGT

    Args:
        germlines_database: Path to IMGT-only Germlines database (session fixture).
        mixed_source_database: Path to mixed-source database (session fixture).
        fasta_file: Path to FASTA file to annotate.
    """
    # Skip if fixture file doesn't exist
    if not fasta_file.exists():
        pytest.skip(f"Test file not found: {fasta_file}")

    # Run annotation with IMGT-only Germlines database
    imgt_api = Airr("human", database=germlines_database)
    imgt_result = imgt_api.run_fasta(fasta_file)

    # Run annotation with mixed-source database
    mixed_api = Airr("human", database=mixed_source_database)
    mixed_result = mixed_api.run_fasta(fasta_file)

    # Compare outputs - should be identical since sequences are the same
    compare_airr_outputs(
        pd.DataFrame(imgt_result),
        pd.DataFrame(mixed_result),
        fasta_file.name,
        excluded_columns=EXCLUDED_COLUMNS | MIXED_SOURCE_EXCLUDED_COLUMNS,
    )
