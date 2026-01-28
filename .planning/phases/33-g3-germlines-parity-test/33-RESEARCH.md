# Phase 33 Research: G3-Germlines Parity Test

## Summary

This phase creates an automated parity test validating that G3 and Germlines backends produce identical AIRR output when built with the same alleles. The implementation leverages existing pytest patterns in the codebase: session-scoped fixtures for database builds, pandas DataFrame comparison utilities, and the established `References.from_yaml()` / `Airr(database=path)` APIs.

**Primary Recommendation:** Use pytest session-scoped fixtures with `tmp_path_factory` for database builds, `pd.testing.assert_frame_equal` for initial column validation, then iterate cell-by-cell for fail-fast mismatch reporting.

---

## Standard Stack

| Component | Library/API | Version | Usage |
|-----------|-------------|---------|-------|
| Test Framework | pytest | (project version) | Test runner with session fixtures |
| DataFrame Comparison | pandas.testing | (project version) | `assert_frame_equal` for column presence checks |
| Database Build | `References.from_yaml()` | SADIE internal | Build G3 vs Germlines databases |
| Annotation | `Airr(database=path)` | SADIE internal | Run annotation with prebuilt database |
| Temp Directories | `tmp_path_factory` | pytest built-in | Session-scoped temporary directories |

**Confidence Level:** HIGH (all patterns verified in existing codebase)

---

## Architecture Patterns

### 1. Session-Scoped Fixture Pattern

The codebase uses session-scoped fixtures for expensive operations. From `tests/conftest.py`:

```python
@pytest.fixture(scope="session", autouse=True)
def fixture_setup(tmp_path_factory: pytest.TempPathFactory):
    return SadieFixture(tmp_path_factory)
```

**For this phase:**

```python
# tests/migration/conftest.py
import pytest
from pathlib import Path
from sadie.reference import References

@pytest.fixture(scope="session")
def g3_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference.g3.yml using G3 backend."""
    yaml_path = Path("reference.g3.yml")
    refs = References.from_yaml(yaml_path, use_germlines=False)
    outpath = tmp_path_factory.mktemp("g3_db")
    return refs.make_airr_database(outpath)

@pytest.fixture(scope="session")
def germlines_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference.g3.yml using Germlines backend."""
    yaml_path = Path("reference.g3.yml")
    refs = References.from_yaml(yaml_path, use_germlines=True)
    outpath = tmp_path_factory.mktemp("germlines_db")
    return refs.make_airr_database(outpath)
```

**Confidence Level:** HIGH (exact pattern used in `tests/conftest.py` and `test_reference_integration.py`)

### 2. Database Parameter for Airr

The `Airr` class accepts a `database` parameter to use prebuilt databases:

```python
# From src/sadie/airr/airr.py:
database: Optional[Path | str] = None,
    """Path to prebuilt database from `sadie reference build`.
    When provided, uses database directly without germlines/G3 lookup.
    Expected structure: Ig/blastdb/, Ig/internal_data/, aux_db/."""
```

**Usage:**
```python
from sadie.airr import Airr

# Use prebuilt database
airr_api = Airr("human", database=db_path)
result = airr_api.run_fasta(fasta_path)
```

**Confidence Level:** HIGH (documented in source code)

### 3. Excluded Columns Pattern

The `_add_source_columns` method adds source tracking columns:

```python
EXCLUDED_COLUMNS = [
    "v_call_source",
    "d_call_source", 
    "j_call_source",
    "c_call_source",
]
```

These columns differ by design (G3 vs Germlines) and must be excluded from comparison.

**Confidence Level:** HIGH (verified in source code)

---

## Don't Hand-Roll

| Problem | Existing Solution | Why |
|---------|------------------|-----|
| Temporary directories | `tmp_path_factory.mktemp()` | Pytest manages cleanup; session-scoped |
| DataFrame column presence | `set(df.columns)` comparison | Simple, reliable |
| NaN equality | `pd.isna()` check before comparison | Built-in pandas support |
| Database build | `References.from_yaml(...).make_airr_database()` | Existing tested API |
| Annotation | `Airr(database=path).run_fasta()` | Existing tested API |

---

## Common Pitfalls

### 1. NaN Handling in Comparisons

**Problem:** Direct `==` comparison fails for NaN values.

**Solution:**
```python
def values_equal(v1, v2):
    """Compare values treating NaN == NaN as True."""
    if pd.isna(v1) and pd.isna(v2):
        return True
    if pd.isna(v1) or pd.isna(v2):
        return False
    return v1 == v2
```

**Confidence Level:** HIGH (decision from 33-CONTEXT.md)

### 2. Column Order May Differ

**Problem:** DataFrames may have same columns but different order.

**Solution:** Compare column sets before comparing values:
```python
g3_cols = set(g3_df.columns)
germlines_cols = set(germlines_df.columns)
assert g3_cols == germlines_cols, f"Column mismatch: G3-only={g3_cols - germlines_cols}, Germlines-only={germlines_cols - g3_cols}"
```

**Confidence Level:** HIGH (standard pandas practice)

### 3. Row Alignment by sequence_id

**Problem:** Row order may differ between backends.

**Solution:** Sort by `sequence_id` and reset index before comparison:
```python
g3_df = g3_df.sort_values("sequence_id").reset_index(drop=True)
germlines_df = germlines_df.sort_values("sequence_id").reset_index(drop=True)
```

**Confidence Level:** MEDIUM (verify sequence_id is reliable key)

### 4. reference_name vs database Parameter

**Problem:** When using `database=` parameter, still need a reference name.

**Solution:** Use "human" as reference_name even with custom database:
```python
airr_api = Airr("human", database=db_path)
```

The reference_name is used for species-specific settings; the database path overrides germline sources.

**Confidence Level:** HIGH (verified in source code)

### 5. Path Type Handling

**Problem:** `make_airr_database()` returns Path; `Airr(database=)` accepts Path or str.

**Solution:** Use Path objects consistently:
```python
db_path: Path = refs.make_airr_database(outpath)
airr_api = Airr("human", database=db_path)
```

**Confidence Level:** HIGH (type hints in source)

---

## Code Examples

### Complete conftest.py

```python
# tests/migration/conftest.py
"""Session-scoped fixtures for G3-Germlines parity testing."""
import pytest
from pathlib import Path
from sadie.reference import References

YAML_PATH = Path("reference.g3.yml")


@pytest.fixture(scope="session")
def g3_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference.g3.yml using G3 backend."""
    refs = References.from_yaml(YAML_PATH, use_germlines=False)
    outpath = tmp_path_factory.mktemp("g3_db")
    return refs.make_airr_database(outpath)


@pytest.fixture(scope="session")
def germlines_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference.g3.yml using Germlines backend."""
    refs = References.from_yaml(YAML_PATH, use_germlines=True)
    outpath = tmp_path_factory.mktemp("germlines_db")
    return refs.make_airr_database(outpath)
```

### Parity Test Structure

```python
# tests/migration/valid_parity.py
"""Parity test validating G3 and Germlines produce identical AIRR output."""
from pathlib import Path
from typing import List

import pandas as pd
import pytest

from sadie.airr import Airr

EXCLUDED_COLUMNS = [
    "v_call_source",
    "d_call_source",
    "j_call_source",
    "c_call_source",
]


def values_equal(v1, v2) -> bool:
    """Compare values treating NaN == NaN as True."""
    if pd.isna(v1) and pd.isna(v2):
        return True
    if pd.isna(v1) or pd.isna(v2):
        return False
    return v1 == v2


def compare_airr_outputs(g3_df: pd.DataFrame, germlines_df: pd.DataFrame, fasta_name: str) -> None:
    """Compare two AIRR DataFrames, fail fast on first mismatch."""
    # 1. Check column presence (excluding source columns)
    g3_cols = set(g3_df.columns) - set(EXCLUDED_COLUMNS)
    germlines_cols = set(germlines_df.columns) - set(EXCLUDED_COLUMNS)
    
    if g3_cols != germlines_cols:
        g3_only = g3_cols - germlines_cols
        germlines_only = germlines_cols - g3_cols
        pytest.fail(f"Column mismatch in {fasta_name}: G3-only={g3_only}, Germlines-only={germlines_only}")
    
    # 2. Sort and align rows
    g3_df = g3_df.sort_values("sequence_id").reset_index(drop=True)
    germlines_df = germlines_df.sort_values("sequence_id").reset_index(drop=True)
    
    if len(g3_df) != len(germlines_df):
        pytest.fail(f"Row count mismatch in {fasta_name}: G3={len(g3_df)}, Germlines={len(germlines_df)}")
    
    # 3. Compare cell by cell, excluding source columns
    compare_cols = sorted(g3_cols)
    
    for row_idx in range(len(g3_df)):
        for col in compare_cols:
            g3_val = g3_df.loc[row_idx, col]
            germlines_val = germlines_df.loc[row_idx, col]
            
            if not values_equal(g3_val, germlines_val):
                seq_id = g3_df.loc[row_idx, "sequence_id"]
                pytest.fail(
                    f"\nPARITY MISMATCH DETECTED\n"
                    f"\nSequence: {fasta_name}\n"
                    f"Row: {row_idx}\n"
                    f"Column: {col}\n"
                    f"Sequence ID: {seq_id}\n"
                    f"\nG3 value:      {g3_val!r}\n"
                    f"Germlines value: {germlines_val!r}"
                )


def get_human_test_fastas() -> List[Path]:
    """Return list of human FASTA files for testing."""
    fixtures_dir = Path("tests/data/fixtures/fasta_inputs")
    return [
        fixtures_dir / "PG9_H.fasta",
        fixtures_dir / "PG9_L.fasta",
        fixtures_dir / "OAS_subsample_1000.fasta",
        fixtures_dir / "catnap_nt_heavy.fasta",
        fixtures_dir / "catnap_nt_light.fasta",
    ]


@pytest.mark.parametrize("fasta_file", get_human_test_fastas(), ids=lambda p: p.name)
def test_airr_parity(g3_database: Path, germlines_database: Path, fasta_file: Path) -> None:
    """Test that G3 and Germlines backends produce identical AIRR output."""
    # Run annotation with G3 database
    g3_api = Airr("human", database=g3_database)
    g3_result = g3_api.run_fasta(fasta_file)
    
    # Run annotation with Germlines database
    germlines_api = Airr("human", database=germlines_database)
    germlines_result = germlines_api.run_fasta(fasta_file)
    
    # Compare outputs
    compare_airr_outputs(
        pd.DataFrame(g3_result),
        pd.DataFrame(germlines_result),
        fasta_file.name
    )
```

---

## Sources

| Source | Confidence | Notes |
|--------|------------|-------|
| `tests/conftest.py` | HIGH | Session-scoped fixture patterns |
| `tests/integration/reference/test_reference_integration.py` | HIGH | Database build patterns, `pd._testing.assert_series_equal` usage |
| `src/sadie/reference/reference.py` | HIGH | `References.from_yaml(path, use_germlines=bool)` API |
| `src/sadie/airr/airr.py` | HIGH | `Airr(reference_name, database=path)` API |
| `src/sadie/airr/airrtable/constants.py` | HIGH | Column definitions |
| `33-CONTEXT.md` | HIGH | Decision rationale for excluded columns, NaN handling, mismatch reporting |
| pytest documentation | HIGH | `tmp_path_factory` and session scope patterns |

---

## Verification Checklist

- [x] All domains investigated (database build, annotation, comparison)
- [x] Negative claims verified (source columns differ by design - verified in code)
- [x] Multiple sources for critical claims (fixtures pattern in conftest.py AND test_reference_integration.py)
- [x] Confidence levels assigned honestly
- [x] Section names match what plan-phase expects

---

*Created: 2026-01-28*
