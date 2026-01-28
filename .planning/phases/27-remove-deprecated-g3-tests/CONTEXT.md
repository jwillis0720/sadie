# Phase 27 Context: Remove Deprecated G3 Tests

## Decisions

### 1. Test Removal Scope

**Decision:** Remove only the skipped `TestGermlineDataPaths` class

- **Remove:** `TestGermlineDataPaths` class (2 skipped tests: `test_v_gene_dir_attribute_exists`, `test_aux_path_attribute_exists`)
- **Keep:** `TestGermlineDataLegacyAPI` class (feature flag tests for migration mechanism)
- **Keep:** File name as `test_germline_data_legacy.py`
- **Defer:** Orphaned import cleanup to future pass

**Rationale:** Feature flag tests (`_use_germlines_module`) are still needed because users may still use `SADIE_USE_GERMLINES_MODULE=false` to fall back to G3. Keep these until G3 is fully removed.

### 2. Deprecation Documentation Location

**Decision:** Create dedicated deprecation doc at `docs/G3-Deprecation.md`

Document should include:
- G3 API deprecation notice
- Timeline: removal after 2026-06-01
- Migration guidance (use germlines module)
- Environment variable `SADIE_USE_GERMLINES_MODULE` behavior

### 3. G3 Code Removal Timing

**Decision:** This phase removes tests only, not G3 code

- G3 code paths remain until official removal (6 months from now)
- This phase focuses on cleaning up skipped tests
- G3 code removal will be a separate future phase

## Constraints

- Do not remove feature flag tests or migration mechanism tests
- Do not remove G3 code paths
- Do not rename the legacy test file

## Files to Modify

- `tests/unit/germlines/test_germline_data_legacy.py` — Remove `TestGermlineDataPaths` class
- `docs/G3-Deprecation.md` — Create deprecation documentation (new file)
