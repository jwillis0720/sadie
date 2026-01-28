# UAT: Phase 32 — Fix IgBLAST internal_data to Use Combined VDJC File

**Status:** In Progress  
**Date:** 2026-01-27

## Test Summary

| # | Test | Status | Notes |
|---|------|--------|-------|
| 1 | No symlinks in internal_data directories | ⏳ | |
| 2 | Combined FASTA contains V+D+J+C segments | ⏳ | |
| 3 | complete_vdj=True for full-length human sequence | ⏳ | |
| 4 | complete_vdj works for mouse sequences | ⏳ | |
| 5 | Workaround code removed (_recalculate_complete_vdj) | ⏳ | |
| 6 | J_GENE_LENGTHS dictionary removed | ⏳ | |
| 7 | Unit tests pass (17 total: 9 + 8) | ⏳ | |

## Test Details

### Test 1: No symlinks in internal_data directories
**Expected:** Zero symlinks in `src/sadie/germlines/igblast/Ig/internal_data/`
**Result:** 

### Test 2: Combined FASTA contains V+D+J+C segments
**Expected:** `human_V.fasta` contains all segment types combined
**Result:** 

### Test 3: complete_vdj=True for full-length human sequence
**Expected:** Running AIRR annotation on a full-length human antibody returns complete_vdj=True
**Result:** 

### Test 4: complete_vdj works for mouse sequences
**Expected:** Mouse sequences also get correct complete_vdj without hardcoded dictionary
**Result:** 

### Test 5: Workaround code removed
**Expected:** `_recalculate_complete_vdj` method no longer exists in Airr class
**Result:** 

### Test 6: J_GENE_LENGTHS dictionary removed
**Expected:** `J_GENE_LENGTHS` no longer exists in j_gene_data.py
**Result:** 

### Test 7: Unit tests pass
**Expected:** All 17 new tests pass (9 from 32-01, 8 from 32-02)
**Result:** 

---
*Created: 2026-01-27*
