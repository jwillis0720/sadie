# Phase 27: Remove Deprecated G3 Tests — Summary

**Status:** ✅ Complete
**Date:** 2026-01-25

## Tasks Completed

### Task 27-1: Remove TestGermlineDataPaths class ✅
**Commit:** `f8f6962b`

Removed the `TestGermlineDataPaths` class from `tests/unit/germlines/test_germline_data_legacy.py`:
- Deleted 47 lines (the entire class with 2 skipped tests)
- Kept `TestGermlineDataLegacyAPI` class (5 tests for feature flag behavior)

### Task 27-2: Create G3 deprecation documentation ✅
**Commit:** `fae2fd28`

Created `docs/G3-Deprecation.md` with:
- Deprecation timeline (removal after 2026-06-01)
- Migration guidance for users and developers
- `SADIE_USE_GERMLINES_MODULE` environment variable documentation
- Benefits of germlines module over legacy G3 API
- List of what will be removed

### Task 27-3: Update CONCERNS.md ✅
**Commit:** `e666e80d`

Removed the row referencing `test_germline_data_legacy.py:74` from `.planning/codebase/CONCERNS.md` skipped tests table.

## Verification

### Must-Haves Verified

| Requirement | Status |
|-------------|--------|
| No tests reference deprecated G3 API | ✅ `TestGermlineDataPaths` removed |
| All test functionality covered | ✅ Verified equivalent tests exist and pass |
| Clear documentation of G3 deprecation | ✅ `docs/G3-Deprecation.md` created |
| Documentation updated | ✅ CONCERNS.md updated |
| No regressions | ✅ All 5 `TestGermlineDataLegacyAPI` tests pass |

### Test Results

```
tests/unit/germlines/test_germline_data_legacy.py: 5 passed
tests/unit/germlines/test_airr_integration.py::TestGermlineDataPaths::test_germlines_path_enabled: PASSED
tests/unit/airr/test_airr.py::test_germline_init: PASSED
```

## Files Modified

| File | Change |
|------|--------|
| `tests/unit/germlines/test_germline_data_legacy.py` | Removed `TestGermlineDataPaths` class |
| `docs/G3-Deprecation.md` | Created (new) |
| `.planning/codebase/CONCERNS.md` | Removed obsolete test reference |

## Commits

1. `f8f6962b` — refactor(27-1): remove deprecated TestGermlineDataPaths class
2. `fae2fd28` — docs(27-2): add G3 API deprecation documentation
3. `e666e80d` — docs(27-3): remove deleted test reference from CONCERNS.md

---

*Phase completed: 2026-01-25*
