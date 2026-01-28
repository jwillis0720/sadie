# Phase 33: G3-Germlines Parity Test — SUMMARY

## Status: COMPLETE (with findings)

## What Was Done

### Task 33-1-1: Create tests/migration Package
- Created `tests/migration/__init__.py` 
- Commit: `bb33048c`

### Task 33-1-2: Create Session-Scoped Database Fixtures
- Created `tests/migration/conftest.py` with:
  - `g3_database` fixture (scope="session", use_germlines=False)
  - `germlines_database` fixture (scope="session", use_germlines=True)
- Commit: `cfc15a0d`

### Task 33-1-3: Create Parity Validation Test
- Created `tests/migration/test_valid_parity.py` with:
  - `EXCLUDED_COLUMNS` frozenset for source tracking columns
  - `values_equal()` function for NaN-safe comparison
  - `compare_airr_outputs()` function with fail-fast on first mismatch
  - `test_airr_parity()` parametrized test for 5 human FASTA files
- Commit: `aef8aebc`

### Deviation: Human-Only Reference
- **Issue**: Original plan specified using `reference.g3.yml`, but G3 backend fails to build databases for multi-species references due to missing IMGT V-region position data for some genes.
- **Resolution**: Created `tests/migration/reference_parity_test.yml` with human-only genes from IMGT that both backends support.
- Commit: `b23e7fd3`

## Parity Test Findings

**Test infrastructure is working correctly.** The parity test detected a real difference between backends:

```
PARITY MISMATCH DETECTED

Sequence file: PG9_H.fasta
Row: 0
Column: j_cigar
Sequence ID: GU272045.1

G3 value:        '355S9N53M'
Germlines value: '355S9N53M1N'
```

The Germlines backend appends `1N` to the J gene CIGAR string that G3 does not include. This is not a test bug - the test correctly detected a difference between the two backends.

## Files Created/Modified

| File | Change | Commit |
|------|--------|--------|
| `tests/migration/__init__.py` | CREATE | bb33048c |
| `tests/migration/conftest.py` | CREATE | cfc15a0d, b23e7fd3 |
| `tests/migration/test_valid_parity.py` | CREATE | aef8aebc |
| `tests/migration/reference_parity_test.yml` | CREATE | b23e7fd3 |

## Verification Commands

```bash
# All pass:
test -d tests/migration                           # ✓
test -f tests/migration/__init__.py               # ✓
test -f tests/migration/conftest.py               # ✓
test -f tests/migration/test_valid_parity.py      # ✓
grep -c 'scope="session"' tests/migration/conftest.py  # Returns 2 ✓
pytest tests/migration/ --collect-only            # Collects 5 tests ✓

# Parity test runs and detects difference:
pytest tests/migration/test_valid_parity.py -v    # FAILS (by design - detected parity issue)
```

## must_haves Status

| ID | Requirement | Status |
|----|-------------|--------|
| MH-1 | Test directory structure exists | ✅ |
| MH-2 | Session-scoped fixtures build both databases | ✅ |
| MH-3 | G3 fixture uses `use_germlines=False` | ✅ |
| MH-4 | Germlines fixture uses `use_germlines=True` | ✅ |
| MH-5 | Source columns are excluded from comparison | ✅ |
| MH-6 | NaN values treated as equal | ✅ |
| MH-7 | Fail-fast on first mismatch | ✅ |
| MH-8 | Error report includes column name | ✅ |
| MH-9 | Error report includes sequence_id | ✅ |
| MH-10 | Error report includes both values | ✅ |
| MH-11 | Tests can be collected by pytest | ✅ |
| MH-12 | All parametrized tests pass | ❌ (j_cigar difference detected) |

## Recommended Next Steps

1. **Investigate j_cigar difference**: The Germlines backend adds `1N` suffix to J gene CIGAR strings. Determine if this is:
   - A bug in Germlines (remove the extra `1N`)
   - A bug in G3 (should include `1N`)
   - An expected difference to exclude (unlikely, as CIGAR format should be consistent)

2. **Check modified files**: The following uncommitted files may be related:
   - `src/sadie/germlines/builders/aux.py`
   - `src/sadie/germlines/builders/j_gene_data.py`
   - `src/sadie/germlines/g3_adapter.py`

3. **Once parity is achieved**: All 5 test cases should pass, confirming backends are equivalent.
