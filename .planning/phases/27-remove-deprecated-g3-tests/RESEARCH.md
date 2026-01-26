# Phase 27: Remove Deprecated G3 Tests — Research

**Mode:** Internal codebase analysis
**Confidence:** HIGH — All findings verified against source code
**Date:** 2026-01-25

---

## Summary

The deprecated G3 tests should be **removed** (not migrated). They test legacy G3 API behavior that is being deprecated on 2026-06-01. Equivalent test coverage already exists in the germlines module tests.

**Primary Recommendation:** Delete the `TestGermlineDataPaths` class from `test_germline_data_legacy.py`. Retain the `TestGermlineDataLegacyAPI` class as those tests are NOT deprecated and verify feature flag behavior.

---

## Standard Stack

| Component | Tool/Library | Version |
|-----------|-------------|---------|
| Testing | pytest | 8.x |
| Markers | `pytest.mark.skip` | Built-in |
| Mocking | `unittest.mock.patch` | Built-in |
| Environment | `monkeypatch` | pytest fixture |

---

## Architecture Patterns

### Test Removal Pattern

```python
# REMOVE: Entire class is skipped and tests deprecated functionality
@pytest.mark.skip(reason="G3 API deprecated, will be removed after 2026-06-01")
class TestGermlineDataPaths:
    """Tests for legacy G3 API path attributes."""
    def test_v_gene_dir_attribute_exists(self, tmp_path):
        ...
    def test_aux_path_attribute_exists(self, tmp_path):
        ...
```

### Verification Pattern

```python
# Verify equivalent tests exist BEFORE removal
# test_airr_integration.py::TestGermlineDataPaths::test_germlines_path_enabled
def test_germlines_path_enabled(self, monkeypatch):
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
    gd = GermlineData("human")
    assert gd.v_gene_dir.with_suffix(".nhr").exists()  # ← Covers removed test
    assert gd.j_gene_dir.with_suffix(".nhr").exists()  # ← Covers removed test
    assert gd.aux_path.exists()                         # ← Covers removed test
```

---

## Don't Hand-Roll

| Problem | Use Instead |
|---------|-------------|
| Test coverage verification | Run `pytest --collect-only` to list all tests |
| Deprecation timeline tracking | Update `.planning/codebase/CONCERNS.md` |
| Documentation updates | Follow existing patterns in `docs/germlines/migration-guide.md` |

---

## Common Pitfalls

### 1. Removing Tests Without Checking Coverage

**Pitfall:** Removing tests that have no equivalent coverage.

**Mitigation:** The deprecated tests verify:
- `v_gene_dir` attribute exists
- `d_gene_dir` attribute exists
- `j_gene_dir` attribute exists
- `aux_path` attribute exists

**Equivalent coverage verified in:**
- `tests/unit/germlines/test_airr_integration.py::TestGermlineDataPaths::test_germlines_path_enabled` — Tests `v_gene_dir`, `j_gene_dir`, `aux_path` with germlines backend
- `tests/unit/germlines/test_airr_integration.py::TestOfflineOperation::test_germline_data_paths_network_disabled` — Tests same attributes
- `tests/unit/airr/test_airr.py::test_germline_init` — Tests setter validation for all path attributes

### 2. Removing Wrong Tests

**Pitfall:** Removing non-deprecated tests in same file.

**The file contains TWO classes:**
```
TestGermlineDataLegacyAPI    ← NOT deprecated, KEEP
TestGermlineDataPaths        ← Deprecated, REMOVE
```

**Keep these tests in `TestGermlineDataLegacyAPI`:**
- `test_germline_data_import`
- `test_germline_data_has_required_attributes`
- `test_get_available_datasets_returns_set`
- `test_feature_flag_deprecation_warning`
- `test_feature_flag_default_true`

### 3. Missing Documentation Update

**Pitfall:** Not updating deprecation timeline documentation.

**Files to check:**
- `.planning/codebase/CONCERNS.md` — Remove mention of skipped test once removed
- No user-facing docs need changes — deprecation is internal

### 4. Breaking Test Collection

**Pitfall:** Leaving orphaned imports or fixtures after test removal.

**Check:** No fixtures are exclusive to the removed class. The `tmp_path` fixture is pytest built-in.

---

## Code Examples

### File Modification: Remove Deprecated Class

**File:** `tests/unit/germlines/test_germline_data_legacy.py`

**Before (lines 56-104):**
```python
@pytest.mark.skip(reason="G3 API deprecated, will be removed after 2026-06-01")
class TestGermlineDataPaths:
    """Tests for legacy G3 API path attributes.
    
    These tests are skipped because they test deprecated G3 API behavior.
    The G3 API will be removed after 2026-06-01.
    """

    def test_v_gene_dir_attribute_exists(self, tmp_path):
        ...  # ~25 lines

    def test_aux_path_attribute_exists(self, tmp_path):
        ...  # ~25 lines
```

**After:**
```python
# Remove entire TestGermlineDataPaths class (lines 56-104)
# File ends after TestGermlineDataLegacyAPI class
```

### Documentation Update

**File:** `.planning/codebase/CONCERNS.md`

**Current (line 50):**
```markdown
| `tests/unit/germlines/test_germline_data_legacy.py:74` | G3 API deprecated | Low (expected) |
```

**Update to:**
```markdown
# Remove this line — tests no longer exist
```

### Verification Command

```bash
# Run tests to ensure nothing breaks after removal
pytest tests/unit/germlines/test_germline_data_legacy.py -v

# Verify equivalent tests still pass
pytest tests/unit/germlines/test_airr_integration.py::TestGermlineDataPaths -v
pytest tests/unit/airr/test_airr.py::test_germline_init -v

# Full test suite smoke test
pytest tests/unit/germlines/ -v --tb=short
```

---

## Analysis: What These Tests Were Validating

### Test 1: `test_v_gene_dir_attribute_exists`

**Purpose:** Verify `GermlineData` has `v_gene_dir`, `d_gene_dir`, `j_gene_dir` attributes when using **legacy G3 paths**.

**Why it's deprecated:** 
- Tests use `patch("sadie.airr.igblast.germline._use_germlines_module", return_value=False)` to force G3 path
- G3 API is being removed on 2026-06-01
- After removal, `_use_germlines_module=False` code path won't exist

**Equivalent coverage:**
- `test_airr_integration.py::test_germlines_path_enabled` tests same attributes with germlines backend
- `test_airr.py::test_germline_init` tests setter validation for these attributes

### Test 2: `test_aux_path_attribute_exists`

**Purpose:** Verify `GermlineData` has `aux_path` attribute when using **legacy G3 paths**.

**Why it's deprecated:** Same reason as above — tests G3 code path that's being removed.

**Equivalent coverage:**
- `test_airr_integration.py::test_germlines_path_enabled` tests `gd.aux_path.exists()`
- `test_airr.py::test_germline_init` tests `gd.aux_path = ...` setter validation

---

## Sources

| Source | Confidence | Notes |
|--------|------------|-------|
| `tests/unit/germlines/test_germline_data_legacy.py` | HIGH | Direct analysis of deprecated tests |
| `tests/unit/germlines/test_airr_integration.py` | HIGH | Verified equivalent coverage exists |
| `tests/unit/airr/test_airr.py` | HIGH | Verified setter validation tests |
| `src/sadie/airr/igblast/germline.py` | HIGH | Understood G3 vs germlines code paths |
| `.planning/ROADMAP.md` | HIGH | Phase 27 requirements |
| `.planning/codebase/CONCERNS.md` | HIGH | Deprecation tracking |
| `docs/germlines/migration-guide.md` | HIGH | Deprecation timeline 2026-06-01 |

---

## Verification Checklist

- [x] All domains investigated (deprecated tests, equivalent coverage, documentation)
- [x] Claims verified with source code inspection
- [x] Multiple sources confirm equivalent coverage exists
- [x] Confidence levels assigned honestly (all HIGH — internal codebase analysis)
- [x] Section names match plan-phase expectations

---

## Decision Summary

| Action | Target | Reason |
|--------|--------|--------|
| **REMOVE** | `TestGermlineDataPaths` class | Tests deprecated G3 API; equivalent coverage exists |
| **KEEP** | `TestGermlineDataLegacyAPI` class | Tests feature flag behavior; still valid |
| **UPDATE** | `.planning/codebase/CONCERNS.md` | Remove reference to deleted tests |
| **VERIFY** | Run pytest after changes | Ensure no regressions |
