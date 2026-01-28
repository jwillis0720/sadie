# Phase 26: Add AIRR Package Dependency

## Overview

Move the `airr` package from dev dependencies to main dependencies and update test to use regular import.

## Goal

Enable AIRR table validation with the official airr package by making it a runtime dependency.

---

## Plan 26-1: Add AIRR as Main Dependency

```yaml
wave: 1
depends_on: []
files_modified:
  - pyproject.toml
  - tests/unit/airr/test_airr.py
autonomous: true
```

<task id="26-1-1" title="Move airr from dev to main dependencies">
**File:** `pyproject.toml`

**Current state:**
- `airr = "^1.5.0"` is in `[tool.poetry.group.dev.dependencies]` section

**Action:**
1. Add `airr = "^1.5.0"` to `[tool.poetry.dependencies]` section (after `rich = "^14.1.0"`)
2. Remove `airr = "^1.5.0"` from `[tool.poetry.group.dev.dependencies]` section

**Result:** airr becomes a runtime dependency installed with `pip install sadie`
</task>

<task id="26-1-2" title="Update test to use regular import">
**File:** `tests/unit/airr/test_airr.py`

**Current state:**
- Line 608: `airr = pytest.importorskip("airr", reason="airr package not installed")`
- This makes the test skip when airr is not installed

**Action:**
1. Add `import airr` to the imports section at the top of the file (after line 10: `import pytest`)
2. Remove line 608: `airr = pytest.importorskip("airr", reason="airr package not installed")`

**Note:** The internal `sadie.airr` module is distinct from the PyPI `airr` package - no naming conflict.

**Result:** Test uses the always-available airr package directly
</task>

### Verification

```bash
# Verify dependency moved correctly
grep -A 5 '^\[tool.poetry.dependencies\]' pyproject.toml | grep airr
grep -A 5 '^\[tool.poetry.group.dev.dependencies\]' pyproject.toml | grep -v airr

# Verify import updated
grep "^import airr" tests/unit/airr/test_airr.py
grep -c "importorskip.*airr" tests/unit/airr/test_airr.py  # should be 0

# Run the specific test
cd /Users/tmsincomb/sadie && python -m pytest tests/unit/airr/test_airr.py::test_write_and_check_airr -v
```

---

## must_haves

Derived from phase goal "Add airr package to dependencies so AIRR validation test runs":

1. **airr in main dependencies**: `pyproject.toml` has `airr = "^1.5.0"` under `[tool.poetry.dependencies]`
2. **airr NOT in dev dependencies**: `airr` line removed from `[tool.poetry.group.dev.dependencies]`
3. **Test uses regular import**: `import airr` at top of test file
4. **No importorskip**: `pytest.importorskip("airr"...)` removed from test
5. **Test passes**: `test_write_and_check_airr` runs and passes (not skipped)

## Success Criteria

1. `pip install sadie` includes airr package
2. `test_write_and_check_airr` passes (not skipped)
3. AIRR table validation works with official airr package
