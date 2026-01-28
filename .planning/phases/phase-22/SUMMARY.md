# Phase 22: Runtime Usage — SUMMARY

**Status:** Complete
**Date:** 2026-01-23

## Tasks Executed

### Task 1: Add Database Structure Validation ✓
**File:** `src/sadie/airr/igblast/germline.py`

- Added `validate_prebuilt_database(database_path, name)` function
- Validates: Ig/, blastdb/{name}/, internal_data/{name}/, aux_db/
- Returns path dict or raises `FileNotFoundError` with clear message

### Task 2: Update GermlineData for Prebuilt Support ✓
**File:** `src/sadie/airr/igblast/germline.py`

- Added `prebuilt=True` parameter to `GermlineData.__init__()`
- When `prebuilt=True`: validates and uses paths directly
- Skips germlines module and G3 lookup entirely

### Task 3: Add `database` Parameter to Airr ✓
**File:** `src/sadie/airr/airr.py`
**Commit:** `b6d0c7fd`

- Added `database: Optional[Path | str] = None` parameter
- Stores `_database_path` for recursive penalty adaptation calls
- Validates path exists before use
- Creates `GermlineData(..., prebuilt=True)` when database provided

## Requirements Satisfied

| Requirement | Description | Status |
|-------------|-------------|--------|
| RUN-01 | Add `Airr(database=<path>)` parameter | ✓ |
| RUN-02 | Skip germlines/G3 lookup when database provided | ✓ |
| RUN-03 | Validate database structure on load | ✓ |

## Success Criteria Verification

1. ✓ `Airr(database="./db")` uses prebuilt database instead of default
2. ✓ No network calls or germlines lookups when database path provided
3. ✓ Clear error if database path missing required structure
4. ✓ Annotation results identical whether using prebuilt or runtime-built
5. ⚠ Performance not measured (infrastructure for startup timing not added)

## Test Results

```python
# Invalid path raises clear error
>>> Airr("human", database="/nonexistent")
FileNotFoundError: Database path not found: /nonexistent

# Missing structure raises validation error
>>> Airr("foo", database="/tmp/empty_dir")
FileNotFoundError: Invalid prebuilt database at /tmp/empty_dir:
  - Missing Ig/ directory...
  - Missing internal_data/foo/...

# Valid prebuilt works
>>> airr = Airr("short", database="/tmp/test_db")
>>> airr.germline_data.base_dir
/tmp/test_db
```

## Files Modified

| File | Change |
|------|--------|
| `src/sadie/airr/igblast/germline.py` | Add `validate_prebuilt_database()`, `prebuilt` param |
| `src/sadie/airr/airr.py` | Add `database` param, pass to recursive calls |
