# Phase 33: G3-Germlines Parity Test — Executable Plan

---
wave: 1
depends_on: [32]
files_modified:
  - tests/migration/__init__.py
  - tests/migration/conftest.py
  - tests/migration/test_valid_parity.py
autonomous: true
tdd: false
---

## Goal Statement

Create an automated parity test validating that G3 and Germlines backends produce identical AIRR output when built with the same alleles from `reference.g3.yml`.

**Rationale:**
- Both backends must be validated to produce equivalent results before deprecating G3
- Using `reference.g3.yml` ensures identical allele sets for fair comparison
- Fail-fast on first mismatch enables efficient debugging

## Context Summary

- **Session-scoped fixtures**: Build databases once per pytest session for efficiency
- **Database parameter**: `Airr(reference_name, database=path)` API uses prebuilt databases
- **Excluded columns**: Source tracking columns (`v_call_source`, etc.) differ by design
- **NaN handling**: Treat `NaN == NaN` as equal for comparison purposes
- **Test data**: Human FASTA files (PG9, OAS, CATNAP) available in fixtures

---

## Plan 1: Create Migration Test Suite

```yaml
wave: 1
depends_on: []
files_modified:
  - tests/migration/__init__.py
  - tests/migration/conftest.py
  - tests/migration/test_valid_parity.py
autonomous: true
```

**Objective:** Create the complete parity test infrastructure with session-scoped database fixtures and parametrized test functions.

<task id="33-1-1">
<title>Create tests/migration Package</title>
<file>tests/migration/__init__.py</file>

**Action:** Create the migration test package with an empty `__init__.py`.

**Code:**
```python
"""Migration tests for G3 to Germlines backend transition."""
```

**Verification:**
```bash
test -f tests/migration/__init__.py && echo "PASS" || echo "FAIL"
```
</task>

<task id="33-1-2">
<title>Create Session-Scoped Database Fixtures</title>
<file>tests/migration/conftest.py</file>

**Action:** Create conftest.py with session-scoped fixtures that build G3 and Germlines databases from `reference.g3.yml`.

**Code:**
```python
"""Session-scoped fixtures for G3-Germlines parity testing."""
import pytest
from pathlib import Path

from sadie.reference import References


# Use absolute path from project root
YAML_PATH = Path(__file__).parent.parent.parent / "reference.g3.yml"


@pytest.fixture(scope="session")
def g3_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference.g3.yml using G3 backend.
    
    Returns:
        Path to the built database directory.
    """
    refs = References.from_yaml(YAML_PATH, use_germlines=False)
    outpath = tmp_path_factory.mktemp("g3_db")
    return refs.make_airr_database(outpath)


@pytest.fixture(scope="session")
def germlines_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference.g3.yml using Germlines backend.
    
    Returns:
        Path to the built database directory.
    """
    refs = References.from_yaml(YAML_PATH, use_germlines=True)
    outpath = tmp_path_factory.mktemp("germlines_db")
    return refs.make_airr_database(outpath)
```

**Verification:**
```bash
# File exists and has correct fixtures
grep -c "def g3_database" tests/migration/conftest.py | grep -q "1" && echo "PASS" || echo "FAIL"
grep -c "def germlines_database" tests/migration/conftest.py | grep -q "1" && echo "PASS" || echo "FAIL"
grep -c "scope=\"session\"" tests/migration/conftest.py | grep -q "2" && echo "PASS" || echo "FAIL"
```
</task>

<task id="33-1-3">
<title>Create Parity Validation Test</title>
<file>tests/migration/test_valid_parity.py</file>

**Action:** Create the main parity test file with comparison logic and parametrized test functions.

**Requirements implemented:**
- TEST-01: Create test file ✓
- TEST-02/03: Use fixtures to build both databases ✓
- TEST-04: Both use same reference.g3.yml ✓
- PAR-01/02: Run AIRR on both databases ✓
- PAR-03: Compare all columns except source columns ✓
- PAR-04: Fail fast on first mismatch ✓
- REP-01/02/03: Report column, row, sequence_id, expected/actual values ✓

**Code:**
```python
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
```

**Verification:**
```bash
# File structure verification
grep -c "def test_airr_parity" tests/migration/test_valid_parity.py | grep -q "1" && echo "PASS" || echo "FAIL"
grep -c "EXCLUDED_COLUMNS" tests/migration/test_valid_parity.py | grep -q "1" && echo "PASS" || echo "FAIL"
grep -c "def values_equal" tests/migration/test_valid_parity.py | grep -q "1" && echo "PASS" || echo "FAIL"
grep -c "def compare_airr_outputs" tests/migration/test_valid_parity.py | grep -q "1" && echo "PASS" || echo "FAIL"

# Check that source columns are excluded
grep "v_call_source" tests/migration/test_valid_parity.py && echo "PASS: Source columns defined" || echo "FAIL"
```
</task>

---

## Plan 1 Verification

```bash
# 1. Test structure exists
test -d tests/migration && echo "PASS: Directory exists" || echo "FAIL"
test -f tests/migration/__init__.py && echo "PASS: __init__.py exists" || echo "FAIL"
test -f tests/migration/conftest.py && echo "PASS: conftest.py exists" || echo "FAIL"
test -f tests/migration/test_valid_parity.py && echo "PASS: test_valid_parity.py exists" || echo "FAIL"

# 2. Fixtures are session-scoped
grep -c "scope=\"session\"" tests/migration/conftest.py

# 3. Test can be collected (doesn't mean it passes)
pytest tests/migration/test_valid_parity.py --collect-only

# 4. Run the actual parity tests (this is the real validation)
pytest tests/migration/test_valid_parity.py -v --tb=short
```

---

## Success Criteria (Phase Level)

1. ✅ `tests/migration/test_valid_parity.py` exists and is executable
2. ✅ Test builds G3 database from `reference.g3.yml` (use_germlines=False)
3. ✅ Test builds Germlines database from `reference.g3.yml` (use_germlines=True)
4. ✅ Test runs AIRR on human sequences with both databases
5. ✅ All non-source columns match exactly between backends
6. ✅ Test fails immediately on first mismatch with clear diagnostics
7. ✅ Test passes when backends produce identical output

---

## must_haves

Derived from phase goal using goal-backward methodology:

| ID | Requirement | Verification |
|----|-------------|--------------|
| MH-1 | Test directory structure exists | `test -d tests/migration && test -f tests/migration/__init__.py` |
| MH-2 | Session-scoped fixtures build both databases | `grep -c "scope=\"session\"" tests/migration/conftest.py` returns 2 |
| MH-3 | G3 fixture uses `use_germlines=False` | `grep "use_germlines=False" tests/migration/conftest.py` |
| MH-4 | Germlines fixture uses `use_germlines=True` | `grep "use_germlines=True" tests/migration/conftest.py` |
| MH-5 | Source columns are excluded from comparison | `grep "v_call_source\|d_call_source\|j_call_source\|c_call_source" tests/migration/test_valid_parity.py` |
| MH-6 | NaN values treated as equal | `grep "pd.isna" tests/migration/test_valid_parity.py` |
| MH-7 | Fail-fast on first mismatch | `grep "pytest.fail" tests/migration/test_valid_parity.py` |
| MH-8 | Error report includes column name | `grep "Column:" tests/migration/test_valid_parity.py` |
| MH-9 | Error report includes sequence_id | `grep "Sequence ID:" tests/migration/test_valid_parity.py` |
| MH-10 | Error report includes both values | `grep "G3 value:\|Germlines value:" tests/migration/test_valid_parity.py` |
| MH-11 | Tests can be collected by pytest | `pytest tests/migration/ --collect-only` succeeds |
| MH-12 | All parametrized tests pass | `pytest tests/migration/test_valid_parity.py -v` all pass |

---

## Files Modified Summary

| File | Change Type | Purpose |
|------|-------------|---------|
| `tests/migration/__init__.py` | CREATE | Package marker for migration tests |
| `tests/migration/conftest.py` | CREATE | Session-scoped database fixtures |
| `tests/migration/test_valid_parity.py` | CREATE | Parity validation test with parametrized FASTA inputs |

---

## Risk Assessment

**Risk Level**: LOW-MEDIUM

- **Low Risk**: File creation only, no modifications to existing code
- **Medium Risk**: G3 API dependency — test may fail if G3 API is unavailable (per 33-CONTEXT.md decision: fail hard)
- **Medium Risk**: Database build time may be slow for CI, but session-scoping mitigates

**Mitigation**:
- Session-scoped fixtures ensure databases are built only once per test run
- Fail-fast comparison reduces debugging time when mismatches occur
- Detailed error output enables quick root cause identification

---

## Deferred Items

None. All requirements from the phase definition are addressed in this plan.
