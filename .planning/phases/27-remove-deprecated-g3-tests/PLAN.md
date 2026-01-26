# Phase 27: Remove Deprecated G3 Tests — Plan

---
wave: 1
depends_on: [26]
files_modified:
  - tests/unit/germlines/test_germline_data_legacy.py
  - .planning/codebase/CONCERNS.md
  - docs/G3-Deprecation.md
autonomous: true
tdd: false
---

## Context

Two tests are skipped with "G3 API deprecated, will be removed after 2026-06-01":
- `test_v_gene_dir_attribute_exists`
- `test_aux_path_attribute_exists`

Research confirms these should be **removed** (not migrated). Equivalent coverage exists in:
- `tests/unit/germlines/test_airr_integration.py::TestGermlineDataPaths::test_germlines_path_enabled`
- `tests/unit/airr/test_airr.py::test_germline_init`

The file contains two classes:
- `TestGermlineDataLegacyAPI` — **KEEP** (tests feature flag behavior, not deprecated)
- `TestGermlineDataPaths` — **REMOVE** (tests deprecated G3 API paths)

## Tasks

<task id="27-1">
<title>Remove TestGermlineDataPaths class</title>
<description>
Delete the entire `TestGermlineDataPaths` class from `tests/unit/germlines/test_germline_data_legacy.py`.

This class is marked `@pytest.mark.skip(reason="G3 API deprecated, will be removed after 2026-06-01")` and contains:
- `test_v_gene_dir_attribute_exists`
- `test_aux_path_attribute_exists`

Both tests use `patch("sadie.airr.igblast.germline._use_germlines_module", return_value=False)` to force the legacy G3 path — a code path being removed on 2026-06-01.

**Keep `TestGermlineDataLegacyAPI` class** — it tests feature flag behavior and is still valid.

**Action:**
1. Open `tests/unit/germlines/test_germline_data_legacy.py`
2. Remove the entire `TestGermlineDataPaths` class (lines 72-119)
3. Remove the `@pytest.mark.skip` decorator and class definition
4. Ensure the file still has proper structure (imports, `TestGermlineDataLegacyAPI` class)
</description>
<files>
  - tests/unit/germlines/test_germline_data_legacy.py
</files>
<verification>
```bash
# File should have no TestGermlineDataPaths class
grep -c "class TestGermlineDataPaths" tests/unit/germlines/test_germline_data_legacy.py
# Expected: 0

# File should still have TestGermlineDataLegacyAPI class
grep -c "class TestGermlineDataLegacyAPI" tests/unit/germlines/test_germline_data_legacy.py
# Expected: 1

# Tests should pass
pytest tests/unit/germlines/test_germline_data_legacy.py -v
```
</verification>
</task>

<task id="27-2">
<title>Create G3 deprecation documentation</title>
<description>
Create a dedicated deprecation documentation file at `docs/G3-Deprecation.md` to clearly document the G3 API deprecation timeline and migration guidance.

**Document must include:**
1. G3 API deprecation notice
2. Timeline: removal scheduled after 2026-06-01
3. Migration guidance (use germlines module)
4. Environment variable `SADIE_USE_GERMLINES_MODULE` behavior
5. What users need to do before the deadline

**Template structure:**
```markdown
# G3 API Deprecation Notice

## Overview
Brief description of the deprecation.

## Timeline
- **Deprecated:** [date]
- **Removal Date:** 2026-06-01

## Migration Guide
How to migrate from G3 to germlines module.

## Environment Variable
Explanation of SADIE_USE_GERMLINES_MODULE behavior.

## FAQ
Common questions about the transition.
```
</description>
<files>
  - docs/G3-Deprecation.md
</files>
<verification>
```bash
# File exists
test -f docs/G3-Deprecation.md && echo "PASS: File exists" || echo "FAIL: File missing"

# Contains required sections
grep -q "2026-06-01" docs/G3-Deprecation.md && echo "PASS: Contains deadline" || echo "FAIL: Missing deadline"
grep -q "SADIE_USE_GERMLINES_MODULE" docs/G3-Deprecation.md && echo "PASS: Contains env var" || echo "FAIL: Missing env var"
grep -q "Migration" docs/G3-Deprecation.md && echo "PASS: Contains migration section" || echo "FAIL: Missing migration"
```
</verification>
</task>

<task id="27-3">
<title>Update CONCERNS.md to remove reference to deleted tests</title>
<description>
Update `.planning/codebase/CONCERNS.md` to remove the reference to the deleted tests under the "Skipped Tests" section.

**Current line (in table):**
```markdown
| `tests/unit/germlines/test_germline_data_legacy.py:74` | G3 API deprecated | Low (expected) |
```

**Action:**
1. Open `.planning/codebase/CONCERNS.md`
2. Find the "Skipped Tests" table under "## Test Coverage Gaps"
3. Remove the row referencing `test_germline_data_legacy.py:74`
</description>
<files>
  - .planning/codebase/CONCERNS.md
</files>
<verification>
```bash
# Should not reference the deleted test
grep "test_germline_data_legacy.py:74" .planning/codebase/CONCERNS.md
# Expected: no output (line removed)
```
</verification>
</task>

## Verification Criteria

After all tasks complete:

```bash
# 1. No TestGermlineDataPaths class exists
grep -r "class TestGermlineDataPaths" tests/unit/germlines/test_germline_data_legacy.py && echo "FAIL" || echo "PASS"

# 2. TestGermlineDataLegacyAPI still exists
grep -c "class TestGermlineDataLegacyAPI" tests/unit/germlines/test_germline_data_legacy.py | grep -q "1" && echo "PASS" || echo "FAIL"

# 3. Deprecation documentation exists
test -f docs/G3-Deprecation.md && echo "PASS" || echo "FAIL"

# 4. Equivalent tests still pass (coverage verification)
pytest tests/unit/germlines/test_airr_integration.py::TestGermlineDataPaths::test_germlines_path_enabled -v
pytest tests/unit/airr/test_airr.py::test_germline_init -v

# 5. Full germlines test suite passes
pytest tests/unit/germlines/ -v --tb=short

# 6. No references to deleted tests in CONCERNS.md
grep -q "test_germline_data_legacy.py:74" .planning/codebase/CONCERNS.md && echo "FAIL" || echo "PASS"
```

## must_haves

Derived from phase goal "Remove or migrate tests that depend on deprecated G3 API":

1. **No tests reference deprecated G3 API** — `TestGermlineDataPaths` class removed from test file
2. **All test functionality covered** — Equivalent tests in `test_airr_integration.py` and `test_airr.py` verified to exist
3. **Clear documentation of G3 deprecation timeline** — `docs/G3-Deprecation.md` created with deadline, migration guidance, and env var documentation
4. **Documentation updated** — CONCERNS.md no longer references the deleted tests
5. **No regressions** — Remaining tests in `TestGermlineDataLegacyAPI` still pass

## Execution Notes

- This is primarily a deletion phase — minimal new code required
- Keep `TestGermlineDataLegacyAPI` class intact — it tests valid feature flag behavior
- The G3 deprecation timeline (2026-06-01) needs dedicated documentation per success criteria
- All three tasks can be executed in parallel (wave 1)
