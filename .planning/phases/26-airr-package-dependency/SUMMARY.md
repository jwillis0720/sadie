# Phase 26 Summary: Add AIRR Package Dependency

## Overview

Moved the `airr` package from dev dependencies to main dependencies, enabling AIRR table validation with the official airr package at runtime.

## Changes Made

### Task 26-1-1: Move airr from dev to main dependencies
**File:** `pyproject.toml`
- Added `airr = "^1.5.0"` to `[tool.poetry.dependencies]` section (after `rich = "^14.1.0"`)
- Removed `airr = "^1.5.0"` from `[tool.poetry.group.dev.dependencies]` section

### Task 26-1-2: Update test to use regular import  
**File:** `tests/unit/airr/test_airr.py`
- Added `import airr` to imports section at top of file
- Removed `airr = pytest.importorskip("airr", reason="airr package not installed")` from `test_write_and_check_airr`

## Verification Results

```
✓ airr in main dependencies: pyproject.toml line 68
✓ airr NOT in dev dependencies: confirmed removed
✓ import airr at top of test file: confirmed
✓ importorskip removed: grep count = 0
✓ test_write_and_check_airr: PASSED (not skipped)
```

## Commit

```
19b0a855 feat(26-1): add airr to main dependencies
```

## Impact

- `pip install sadie` now includes the airr package as a runtime dependency
- `test_write_and_check_airr` runs and passes (no longer conditionally skipped)
- AIRR table validation using the official airr package works out of the box

## Status

✅ **Complete** — All must_haves satisfied
