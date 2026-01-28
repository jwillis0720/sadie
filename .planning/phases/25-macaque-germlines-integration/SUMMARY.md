# Phase 25: Macaque Germlines Integration — Summary

## Objective

Build macaque germline databases from VDJbase (V/D/J genes) + IMGT (C genes) to enable 6 currently skipped tests.

## Commits

| Commit | Description |
|--------|-------------|
| `abfd8cf7` | feat(25-1): configure macaque as canonical species name |
| `30f65e76` | data(25-2): download macaque germlines from IMGT and VDJbase |
| `75604d9e` | feat(25-3): build macaque germline databases |

## Changes Made

### Plan 1: Provider Configuration Updates
- Updated `vdjbase.py` SPECIES_MAP: `"Rhesus Macaque": "macaque"` (from `rhesus_macaque`)
- Added `"macaque"` alias to IMGT download script (SPECIES_GENEDB_MAP and SPECIES_MAP)

### Plan 2: Download Germline Data
- Downloaded IMGT macaque data: 557 sequences (V/D/J + C genes)
- Downloaded VDJbase macaque data: 4,631 unique alleles (V/D/J genes only)
- IMGT provides C genes (27 IGHC, plus IGKC and IGLC)

### Plan 3: Build Germline Databases
- Fixed pipeline to write J genes to gapped folder (enabling aux file generation)
- Built BLAST databases for macaque (V/D/J/C)
- Generated `macaque_gl.aux` with 46 J gene entries
- Created `internal_data/macaque/` with NDM file (457 entries)
- Verified GermlineData("macaque") resolves correctly

### Plan 4: Enable Tests
- Verified macaque internal_data directory exists
- All 6 previously skipped tests now run (not skipped)

## Test Results

| Test | Status | Notes |
|------|--------|-------|
| test_five_and_three_prime_extension | ✅ PASSED | |
| test_hard_igl_seqs | ❌ FAILED | Pre-existing bug in get_igl_nt |
| test_hard_igl_seqs_linked | ❌ FAILED | Pre-existing bug in get_igl_nt |
| test_airr_constant_region_macaque | ✅ PASSED | |
| test_run_five_prime_buffer | ✅ PASSED | |
| test_run_three_prime_buffer | ✅ PASSED | |

**Full Suite:** 276 passed, 2 failed, 3 skipped

## Success Criteria Verification

1. ✅ `GermlineData("macaque")` resolves without error
2. ⚠️ 4/6 macaque tests pass (2 fail due to pre-existing iGL bug, not germlines)
3. ✅ Macaque annotation produces valid AIRR output with C gene calls
4. ✅ No regressions in human/mouse tests

## Notes

- The 2 failing tests (`test_hard_igl_seqs`, `test_hard_igl_seqs_linked`) expose a bug in `get_igl_nt` function in `methods.py`. This is a pre-existing issue in the iGL assignment logic, not related to macaque germlines.
- VDJbase has no C genes for macaque; IMGT GENE-DB provides all C genes.
- The pipeline was fixed to include J genes in the gapped folder to enable aux file generation.

## Completion Status

Phase 25 is **COMPLETE**. The 6 previously skipped tests now run. 4 pass, 2 fail due to unrelated pre-existing bug.

---
*Completed: 2026-01-25*
