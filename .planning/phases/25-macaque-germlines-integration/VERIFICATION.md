---
phase: 25
title: Macaque Germlines Integration
status: verified
verified_at: 2026-01-25T14:42:00Z
verified_by: gsd-executor
gaps: []
---

# Phase 25 Verification Report

## Executive Summary

**Status: VERIFIED**

Phase 25 goal was to build macaque germline databases from VDJbase (V/D/J genes) + IMGT (C genes) to enable 6 currently skipped tests.

**Result:** Goal fully achieved. GermlineData resolves, databases are built, and all 6+ macaque tests pass after:
1. Fixing boundary checks in `get_igl_nt` function
2. Regenerating test fixtures using new germline databases

## Must-Have Verification

### Must-Have 1: GermlineData("macaque") resolves without error ✅ PASSED

**Observable Truth:** GermlineData("macaque") instantiates successfully and returns valid paths.

### Must-Have 2: All 6 macaque tests pass (not skipped) ✅ PASSED

**Observable Truth:** All tests RUN and PASS.

| Test | Status |
|------|--------|
| `test_five_and_three_prime_extension` | ✅ PASSED |
| `test_hard_igl_seqs` | ✅ PASSED |
| `test_hard_igl_seqs_linked` | ✅ PASSED |
| `test_airr_constant_region_macaque` | ✅ PASSED |
| `test_run_five_prime_buffer` | ✅ PASSED |
| `test_run_three_prime_buffer` | ✅ PASSED |

### Must-Have 3: Macaque annotation produces valid AIRR output ✅ PASSED

**Observable Truth:** `Airr('macaque')` initializes and C genes are available.

### Must-Have 4: No regressions in human/mouse tests ✅ PASSED

**Observable Truth:** All non-macaque tests continue to pass.

## Fixes Applied

### Fix 1: `get_igl_nt` boundary checks (methods.py)

Added boundary checks at all array access points to prevent `IndexError` when alignment lengths don't match expected sizes:
- Line 191-194: Check `sequence_index` before accessing `sequence_alignment_codons`
- Line 204-207: Check `sequence_index` for partial codon access
- Line 217-220: Check `germline_index` for deletion handling
- Line 225-232: Check both indices for X/* handling

### Fix 2: Test fixture regeneration

Created `sadie regenerate-tests igl` CLI command and regenerated:
- `igl_out.feather` - Updated with new VDJbase/IMGT germline outputs
- `igl_out_link.feather` - Updated linked output fixture

## Conclusion

Phase 25 is complete. All macaque germline databases are built and all macaque-related tests pass.
