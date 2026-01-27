# Summary: Plan 32-02 — Remove Phase 17 Workaround Code

## Status: ✅ Complete

## Overview

Removed the Phase 17 workaround code that post-processed IgBLAST results to recalculate `complete_vdj`. This workaround was necessary because IgBLAST's symlink-based internal_data structure couldn't calculate J gene lengths correctly. Plan 32-01 fixed this by creating combined VDJC files, making the workaround obsolete.

## Changes Made

### Task 1: Remove _recalculate_complete_vdj from Airr Class
**Commit:** cb79ad29

Removed from `src/sadie/airr/airr.py`:
- `_recalculate_complete_vdj()` method definition (~47 lines)
- Call in `run_fasta()` method
- Two calls in `_run_scfv()` method
- Import of `get_j_gene_length` (was inside method)

**Lines removed:** 55

### Task 2: Clean Up j_gene_data.py
**Commit:** 3a5d83af

Removed from `src/sadie/germlines/builders/j_gene_data.py`:
- `J_GENE_LENGTHS` dictionary (34 hardcoded human J allele lengths)
- `get_j_gene_length()` function

Preserved for aux file generation:
- `HUMAN_J_GENE_DATA` (reading frame, CDR3 end positions)
- `CHAIN_TYPE_MAP`
- `get_j_gene_data()`

**Lines removed:** 63

### Task 3: Verify complete_vdj Works for Multiple Species
**Commit:** 37738487

Created `tests/unit/airr/test_complete_vdj.py` with 8 tests:
- 4 tests verifying workaround code removed
- 2 tests verifying human complete_vdj works (regression)
- 2 tests verifying mouse complete_vdj works (new capability)

Updated `tests/data/fixtures/airr_tables/igl_out.feather`:
- 1 macaque sequence (ID: 51260) now correctly gets `complete_vdj=True`
- Changed from False→True because IgBLAST can now calculate J gene length

## Verification

```
$ pytest tests/unit/airr/test_complete_vdj.py tests/unit/airr/test_airr.py::test_hard_igl_seqs -v
9 passed
```

## Benefits

1. **Code Simplification:** Removed ~118 lines of workaround code
2. **Multi-Species Support:** `complete_vdj` now works for all 29+ species, not just human
3. **Accuracy:** IgBLAST calculates directly from internal_data, no post-processing needed
4. **Maintainability:** No need to maintain hardcoded J gene length dictionary

## Commit Hashes

| Task | Commit |
|------|--------|
| Task 1 | cb79ad29 |
| Task 2 | 3a5d83af |
| Task 3 | 37738487 |

## Files Modified

- `src/sadie/airr/airr.py` (-55 lines)
- `src/sadie/germlines/builders/j_gene_data.py` (-63 lines)
- `tests/unit/airr/test_complete_vdj.py` (new, +180 lines)
- `tests/data/fixtures/airr_tables/igl_out.feather` (updated)
