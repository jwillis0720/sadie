# Phase 19: Source Validation

**Goal:** Expand source validation to accept all germline database providers

**Requirements:**
- SRC-01: Expand VALID_SOURCES to include `ogrdb`, `vdjbase`
- SRC-02: Validate source exists in germlines before processing (deferred to Phase 20)

## Tasks

### Task 1: Add VALID_SOURCES constant
**File:** `src/sadie/reference/models.py`
**Action:** Add module-level constant at top of file

```python
VALID_SOURCES = ["imgt", "ogrdb", "vdjbase", "custom"]
```

### Task 2: Update GeneEntry.check_source
**File:** `src/sadie/reference/models.py`
**Lines:** 40-42

**Before:**
```python
if v not in ["imgt", "custom"]:
    raise ValueError(f"{v} is not a valid source, chocies are 'imgt' or 'custom'")
```

**After:**
```python
if v not in VALID_SOURCES:
    raise ValueError(f"{v} is not a valid source, choices are {VALID_SOURCES}")
```

### Task 3: Update GeneEntries.check_source
**File:** `src/sadie/reference/models.py`
**Lines:** 68-70

**Before:**
```python
if v not in ["imgt", "custom"]:
    raise ValueError(f"{v} is not a valid source, chocies are 'imgt' or 'custom'")
```

**After:**
```python
if v not in VALID_SOURCES:
    raise ValueError(f"{v} is not a valid source, choices are {VALID_SOURCES}")
```

### Task 4: Add unit tests
**File:** `tests/unit/reference/test_models.py` (create if needed)
**Action:** Test all four sources pass validation, invalid sources fail

### Task 5: Run existing tests
**Command:** `pytest tests/unit/reference/ -v`

## Success Criteria
1. ✓ `models.py` accepts `source: ogrdb` and `source: vdjbase`
2. ✓ Existing `imgt` and `custom` sources continue to work
3. ✓ Error message shows all valid choices
4. ✓ Unit tests pass

## Decisions Made (from discussion)
- SRC-02 deferred to Phase 20 (validate at add_genes() time, not model time)
- Module-level VALID_SOURCES constant for easy future expansion
- Fix typo "chocies" → "choices"

---
*Created: 2026-01-23*
