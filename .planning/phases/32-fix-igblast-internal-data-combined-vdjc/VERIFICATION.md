---
status: passed
verified_at: "2026-01-27"
plans_verified:
  - 32-01
  - 32-02
must_haves_verified: 16
must_haves_total: 16
gaps: []
human_verification_needed: []
---

# Phase 32 Verification: Fix IgBLAST internal_data Combined VDJC

## Verification Status: ✅ PASSED

All 16 must-haves verified against actual codebase implementation.

---

## Plan 32-01: Create Combined VDJC Files for internal_data

### Must-Have 1: build_internal_data.py creates combined VDJC FASTA (not symlinks)
**Status:** ✅ Verified

**Evidence:**
- File: `src/sadie/germlines/scripts/build_internal_data.py`
- Function `build_combined_fasta()` concatenates V+D+J+C segments
- Deduplication logic present (tracks `seen_ids` set)
- No symlink creation code found

### Must-Have 2: Combined file named `{species}_V.fasta` in `internal_data/{species}/`
**Status:** ✅ Verified

**Evidence:**
```python
# From build_internal_data.py:40
combined_fasta = internal_data_dir / f"{species}_V.fasta"
```
- File exists: `src/sadie/germlines/igblast/Ig/internal_data/human/human_V.fasta`
- Contains 881 unique sequences

### Must-Have 3: BLAST database built from combined file with all segments
**Status:** ✅ Verified

**Evidence:**
- BLAST database files exist in `internal_data/human/`:
  - `human_V.nhr`, `human_V.nin`, `human_V.nsq`, etc.
- `build_blast_db()` function calls makeblastdb with combined FASTA

### Must-Have 4: No symlinks exist in any internal_data/{species}/ directory
**Status:** ✅ Verified

**Evidence:**
```bash
$ find src/sadie/germlines/igblast/Ig/internal_data -type l | wc -l
0
```

### Must-Have 5: NDM file generation unchanged (V genes only)
**Status:** ✅ Verified

**Evidence:**
- `build_ndm_file()` only processes files matching `IG{chain}V_gapped.fasta`
- NDM file exists: `internal_data/human/human.ndm.imgt`

### Must-Have 6: Reference builder includes D/J/C in BLAST db
**Status:** ✅ Verified

**Evidence:**
```python
# From reference.py:_make_internal_annotaion_file()
all_segments = group_df.loc[group_df["gene_segment"].isin(["V", "D", "J", "C"])].copy()
# ...
make_blast_db_for_internal(all_segments, db_outpath)
```

### Must-Have 7: GermlineData points V/D/J/C dirs to `database/{species}/`
**Status:** ✅ Verified

**Evidence:**
```python
# From germline.py:GermlineData.__init__()
blast_prefix = database_species / f"{name}_"
self.v_gene_dir = Path(str(blast_prefix) + "V")
self.d_gene_dir = Path(str(blast_prefix) + "D")
self.j_gene_dir = Path(str(blast_prefix) + "J")
self.c_gene_dir = Path(str(blast_prefix) + "C")
```

### Must-Have 8: Tests verify new structure and no symlinks
**Status:** ✅ Verified

**Evidence:**
- Test file: `tests/unit/germlines/test_build_internal_data.py`
- 9 tests covering:
  - `test_no_symlinks_in_internal_data` 
  - `test_combined_fasta_exists_for_human`
  - `test_combined_fasta_contains_multiple_segments`
  - `test_blast_database_files_exist`
  - `test_germline_data_paths_point_to_database`
  - And 4 more

---

## Plan 32-02: Remove Phase 17 Workaround Code

### Must-Have 1: `_recalculate_complete_vdj()` method removed from Airr class
**Status:** ✅ Verified

**Evidence:**
- Grep search for `_recalculate_complete_vdj` in `src/` returned 0 matches
- Method not present in `src/sadie/airr/airr.py`
- Test confirms: `test_recalculate_complete_vdj_method_removed`

### Must-Have 2: All calls to `_recalculate_complete_vdj()` removed
**Status:** ✅ Verified

**Evidence:**
- No calls found in `run_fasta()` or `_run_scfv()` methods
- Full review of `airr.py` confirms removal

### Must-Have 3: `J_GENE_LENGTHS` dictionary removed from j_gene_data.py
**Status:** ✅ Verified

**Evidence:**
- Not present in `src/sadie/germlines/builders/j_gene_data.py`
- Test confirms: `test_j_gene_lengths_dict_removed`

### Must-Have 4: `get_j_gene_length()` function removed from j_gene_data.py
**Status:** ✅ Verified

**Evidence:**
- Not present in `src/sadie/germlines/builders/j_gene_data.py`
- Test confirms: `test_get_j_gene_length_function_removed`

### Must-Have 5: `HUMAN_J_GENE_DATA` and `get_j_gene_data()` preserved
**Status:** ✅ Verified

**Evidence:**
```python
# Present in j_gene_data.py:
HUMAN_J_GENE_DATA = {
    "IGHJ1*01": (0, 17, 1),
    # ... 34 more entries
}

def get_j_gene_data(allele_name: str, chain: str) -> tuple:
    # ...
```

### Must-Have 6: No import errors after cleanup
**Status:** ✅ Verified

**Evidence:**
- All 17 tests pass without import errors
- `from sadie.airr import Airr` works
- `from sadie.germlines.builders import j_gene_data` works

### Must-Have 7: `complete_vdj` works correctly for human (regression)
**Status:** ✅ Verified

**Evidence:**
- Tests pass: `TestCompleteVdjHuman` (2 tests)
- `complete_vdj` column populated by IgBLAST

### Must-Have 8: `complete_vdj` works for non-human species (new capability)
**Status:** ✅ Verified

**Evidence:**
- Tests pass: `TestCompleteVdjMouse` (2 tests)
- Mouse annotation works without workaround

---

## Test Execution Summary

```
$ pytest tests/unit/germlines/test_build_internal_data.py tests/unit/airr/test_complete_vdj.py -v
17 passed in 5.40s
```

### Tests Executed:
| Test Class | Tests | Status |
|------------|-------|--------|
| TestInternalDataStructure | 5 | ✅ All passed |
| TestGermlineDataPaths | 2 | ✅ All passed |
| TestBuildInternalDataScript | 2 | ✅ All passed |
| TestCompleteVdjWorkaroundRemoved | 4 | ✅ All passed |
| TestCompleteVdjHuman | 2 | ✅ All passed |
| TestCompleteVdjMouse | 2 | ✅ All passed |

---

## Observable Truths Verified

1. **Combined FASTA contains all segment types:**
   - V genes: IGHV, IGKV, IGLV ✓
   - D genes: IGHD ✓
   - J genes: IGHJ, IGKJ, IGLJ ✓
   - C genes: IGHA, IGHE, IGHG, IGHM, IGKC, IGLC ✓

2. **No symlinks in internal_data:**
   - `find ... -type l | wc -l` returns 0

3. **IgBLAST calculates complete_vdj correctly:**
   - Human sequences: working ✓
   - Mouse sequences: working ✓
   - No post-processing workaround needed ✓

---

## Conclusion

Phase 32 goal achieved: IgBLAST internal_data restructured to use combined VDJC files (not symlinks), and the Phase 17 workaround code has been completely removed. The system now relies on IgBLAST's native `complete_vdj` calculation for all species.
