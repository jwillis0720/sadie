# Phase 29 Plan 01 Summary: Germline Source Tracking

## Status: ✅ Complete

## Commit
- Hash: `e3129873`
- Message: `feat(29-01): add germline source tracking columns to AIRR output`

## Tasks Completed

### Task 1: Add get_source_lookup() to GermlineData
- **File:** `src/sadie/airr/igblast/germline.py`
- **Changes:**
  - Added `lru_cache` import from functools
  - Added `get_source_lookup()` method with LRU caching
  - Method iterates over V/D/J/C segments and H/K/L chains
  - Uses GermlineManager to get genes and builds lookup dict mapping gene name → source

### Task 2: Add _lookup_source() to Airr
- **File:** `src/sadie/airr/airr.py`
- **Changes:**
  - Added `numpy` and `Dict` imports
  - Added `_lookup_source()` helper method
  - Handles edge cases: NaN calls → np.nan, comma-separated → first allele, missing → "unknown"

### Task 3: Add _add_source_columns() and integrate into run_fasta()
- **File:** `src/sadie/airr/airr.py`
- **Changes:**
  - Added `_add_source_columns()` method that adds v_call_source, d_call_source, j_call_source, c_call_source columns
  - Integrated before AirrTable() conversion in run_fasta()

### Task 4: Add source columns in _run_scfv()
- **File:** `src/sadie/airr/airr.py`
- **Changes:**
  - Added `_add_source_columns()` calls before creating AirrTables
  - Source columns are added before merge, so LinkedAirrTable gets _heavy/_light suffixed columns

## Verification Results

| Check | Status |
|-------|--------|
| get_source_lookup() returns 1730 genes | ✅ |
| _lookup_source() handles NaN, comma-separated, unknown | ✅ |
| _add_source_columns() method exists | ✅ |
| run_fasta() produces output with source columns | ✅ |
| _run_scfv() produces LinkedAirrTable with suffixed source columns | ✅ |
| Syntax check passed | ✅ |
| Imports work correctly | ✅ |

## Output Example

**Standard AIRR output columns:**
- `v_call_source`: vdjbase/imgt/ogrdb/custom/unknown
- `d_call_source`: vdjbase/imgt/ogrdb/custom/unknown
- `j_call_source`: vdjbase/imgt/ogrdb/custom/unknown
- `c_call_source`: vdjbase/imgt/ogrdb/custom/unknown

**scfv LinkedAirrTable columns:**
- `v_call_source_heavy`, `v_call_source_light`
- `d_call_source_heavy`, `d_call_source_light`
- `j_call_source_heavy`, `j_call_source_light`
- `c_call_source_heavy`, `c_call_source_light`

## Files Modified
- `src/sadie/airr/igblast/germline.py` (+30 lines)
- `src/sadie/airr/airr.py` (+60 lines)

## Deviations
None - plan executed as specified.

---
*Generated: 2026-01-25*
