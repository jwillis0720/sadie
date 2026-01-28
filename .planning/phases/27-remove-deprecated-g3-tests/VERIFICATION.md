---
status: gaps_found
verified_at: "2026-01-25T19:28:00Z"
verifier: gsd-verifier
gaps:
  - id: gap-27-1
    must_have: "Documentation updated — CONCERNS.md no longer references the deleted tests"
    evidence: "Line 46 still contains: '- tests/unit/germlines/test_germline_data_legacy.py:74 — Test skipped due to deprecation'"
    severity: low
    fix_action: "Remove line 46 reference from .planning/codebase/CONCERNS.md"
human_verification_needed: []
---

# Phase 27 Verification: Remove Deprecated G3 Tests

## Goal
Remove or migrate tests that depend on deprecated G3 API

## Verification Date
2026-01-25

## Must-Have Verification

### 1. No tests reference deprecated G3 API ✅ PASSED
**Requirement:** `TestGermlineDataPaths` class removed from test file

**Evidence:**
```bash
$ grep "class TestGermlineDataPaths" tests/unit/germlines/test_germline_data_legacy.py
# No matches found
```

**Verification:** The file `tests/unit/germlines/test_germline_data_legacy.py` no longer contains the deprecated `TestGermlineDataPaths` class. Only `TestGermlineDataLegacyAPI` remains, which tests valid feature flag behavior.

---

### 2. All test functionality covered ✅ PASSED
**Requirement:** Equivalent tests in `test_airr_integration.py` and `test_airr.py` verified to exist

**Evidence:**

**test_airr_integration.py::TestGermlineDataPaths (line 148):**
```python
class TestGermlineDataPaths:
    """Test GermlineData path switching with feature flag."""

    def test_germlines_path_enabled(self, monkeypatch):
        """Verify GermlineData uses germlines module paths when enabled."""
        # Tests v_gene_dir, j_gene_dir, aux_path existence

    def test_legacy_path_disabled(self, monkeypatch):
        """Verify GermlineData uses legacy paths when disabled."""
```

**test_airr.py::test_germline_init (line 92):**
```python
def test_germline_init() -> None:
    # Tests GermlineData initialization
```

**Test Execution Results:**
```
tests/unit/germlines/test_airr_integration.py::TestGermlineDataPaths::test_germlines_path_enabled PASSED
tests/unit/germlines/test_airr_integration.py::TestGermlineDataPaths::test_legacy_path_disabled PASSED
tests/unit/airr/test_airr.py::test_germline_init PASSED
```

---

### 3. Clear documentation of G3 deprecation timeline ✅ PASSED
**Requirement:** `docs/G3-Deprecation.md` created with deadline, migration guidance, and env var documentation

**Evidence:**
- File exists: `docs/G3-Deprecation.md` ✓
- Contains deadline: "2026-06-01" ✓
- Contains migration guidance: "Migration Path" section ✓
- Contains env var documentation: "SADIE_USE_GERMLINES_MODULE" ✓

**File Contents Verified:**
- Timeline table with deprecation/removal dates
- User migration instructions
- Developer migration code examples
- Environment variable behavior table
- Deprecation warning explanation

---

### 4. Documentation updated ❌ GAP FOUND
**Requirement:** CONCERNS.md no longer references the deleted tests

**Evidence:**
```bash
$ grep "test_germline_data_legacy.py:74" .planning/codebase/CONCERNS.md
- `tests/unit/germlines/test_germline_data_legacy.py:74` — Test skipped due to deprecation
```

**Gap:** Line 46 of `.planning/codebase/CONCERNS.md` still references the deleted test path.

**Expected:** This reference should have been removed per task 27-3.

**SUMMARY.md Claim (Commit e666e80d):**
> "Updated .planning/codebase/CONCERNS.md to remove the reference"

**Actual:** Reference still exists in the file.

---

### 5. No regressions ✅ PASSED
**Requirement:** Remaining tests in `TestGermlineDataLegacyAPI` still pass

**Evidence:**
```
$ poetry run pytest tests/unit/germlines/test_germline_data_legacy.py -v
tests/unit/germlines/test_germline_data_legacy.py::TestGermlineDataLegacyAPI::test_germline_data_import PASSED
tests/unit/germlines/test_germline_data_legacy.py::TestGermlineDataLegacyAPI::test_germline_data_has_required_attributes PASSED
tests/unit/germlines/test_germline_data_legacy.py::TestGermlineDataLegacyAPI::test_get_available_datasets_returns_set PASSED
tests/unit/germlines/test_germline_data_legacy.py::TestGermlineDataLegacyAPI::test_feature_flag_deprecation_warning PASSED
tests/unit/germlines/test_germline_data_legacy.py::TestGermlineDataLegacyAPI::test_feature_flag_default_true PASSED

============================== 5 passed in 0.09s ===============================
```

---

## ROADMAP Success Criteria Verification

| Criterion | Status |
|-----------|--------|
| No tests reference deprecated G3 API | ✅ PASSED |
| All test functionality covered by germlines module tests | ✅ PASSED |
| Clear documentation of G3 deprecation timeline | ✅ PASSED |

---

## Anti-Pattern Scan

| Pattern | Found | Notes |
|---------|-------|-------|
| Stubs without implementation | No | N/A (deletion phase) |
| Dead code references | Yes | CONCERNS.md line 46 references deleted test |
| Incomplete wiring | No | All replacement tests functional |

---

## Structured Gaps

```yaml
gaps:
  - id: gap-27-1
    must_have: "Documentation updated — CONCERNS.md no longer references the deleted tests"
    evidence: "Line 46 still contains reference to test_germline_data_legacy.py:74"
    severity: low
    fix_action: "Remove line 46 from .planning/codebase/CONCERNS.md Deprecation section"
```

---

## Overall Status: GAPS_FOUND

**Summary:**
- 4 of 5 must-haves verified ✅
- 1 documentation gap found (low severity)
- All code changes verified correct
- All tests pass
- Gap requires simple line deletion from CONCERNS.md

**Fix Required:**
Remove line 46 from `.planning/codebase/CONCERNS.md`:
```diff
- - `tests/unit/germlines/test_germline_data_legacy.py:74` — Test skipped due to deprecation
```

---

*Verification performed by gsd-verifier subagent*
