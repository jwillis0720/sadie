---
phase: 25
type: gap_closure
gap_id: GAP-25-1
wave: 1
---

# Phase 25 Gap Closure: Regenerate Test Fixtures

## Gap Being Closed

**GAP-25-1:** Test Fixture Incompatibility

2 iGL tests fail due to test fixtures generated with old G3 database being incompatible with new VDJbase/IMGT germlines:
- `test_hard_igl_seqs`
- `test_hard_igl_seqs_linked`

## Root Cause

Test fixtures (`bum_igl_assignment_macaque.feather`, `bum_link_input.feather`) were generated using old G3 macaque databases. The new VDJbase/IMGT databases produce different V gene alignments, causing `IndexError` in `get_igl_nt` when `sequence_index` exceeds `sequence_alignment_codons` length.

## Solution

Regenerate test fixtures using the new macaque germline databases. Existing scripts handle this:
- `scripts/regenerate_igl_reference.py`
- `scripts/regenerate_linked_igl_reference.py`

---

## Plan 1: Regenerate Test Fixtures

### Task 1.1: Regenerate IGL Reference Fixture

**Action:** Run regenerate_igl_reference.py to update igl_out.feather

**Commands:**
```bash
cd /Users/tmsincomb/sadie
python scripts/regenerate_igl_reference.py --verbose
```

**Expected Output:**
- Updated `tests/data/fixtures/airr_tables/igl_out.feather`
- Backup created of previous version

**Verification:**
```bash
ls -la tests/data/fixtures/airr_tables/igl_out.feather
```

### Task 1.2: Regenerate Linked IGL Reference Fixture

**Action:** Run regenerate_linked_igl_reference.py to update bum_link_solution.feather

**Commands:**
```bash
cd /Users/tmsincomb/sadie
python scripts/regenerate_linked_igl_reference.py --verbose
```

**Expected Output:**
- Updated `tests/data/fixtures/airr_tables/bum_link_solution.feather`
- Backup created of previous version

**Verification:**
```bash
ls -la tests/data/fixtures/airr_tables/bum_link_solution.feather
```

### Task 1.3: Run Previously Failing Tests

**Action:** Verify both tests now pass

**Commands:**
```bash
cd /Users/tmsincomb/sadie
pytest tests/unit/airr/test_airr.py::test_hard_igl_seqs \
       tests/unit/airr/test_airr.py::test_hard_igl_seqs_linked \
       -v
```

**Expected Output:** Both tests PASSED

### Task 1.4: Run Full Test Suite

**Action:** Verify no regressions

**Commands:**
```bash
cd /Users/tmsincomb/sadie
pytest tests/unit/airr/ -v --tb=short
```

**Expected Output:** All tests pass, no regressions

---

## Success Criteria

1. ✅ `test_hard_igl_seqs` passes
2. ✅ `test_hard_igl_seqs_linked` passes
3. ✅ All 6 macaque tests pass (6/6)
4. ✅ No regressions in other tests

---

## Files Modified

- `tests/data/fixtures/airr_tables/igl_out.feather` — Regenerated
- `tests/data/fixtures/airr_tables/bum_link_solution.feather` — Regenerated
