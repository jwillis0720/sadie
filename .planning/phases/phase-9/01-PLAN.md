# Phase 9 Plan: Compliance & Coverage Gaps

## Objective

Close remaining requirement gaps (16%) by implementing validation, error handling, and compliance tests for the germlines integration.

## Context

Phase 9 tasks from ROADMAP.md:
- T053: Single-provider enforcement (FR-014)
- T054: Tests rejecting per-segment providers (FR-014)
- T055: Clear error when provider lacks species (FR-006, NFR-002)
- T056: Custom germline ingestion validation (FR-012)
- T057: Species/chain/segment parity verification (FR-010)
- T058: Default priority order test (FR-004)
- T059: No G3 fallback negative test (NFR-002)
- T060: Fail-fast for missing gapped sequences (FR-013)

## Tasks

### Group A: Priority Order Fix (T058)

**Files**: `src/sadie/germlines/manager.py`

**Issue**: `DEFAULT_PROVIDERS = ["custom", "imgt", "ogrdb", "vdjbase"]`
**Spec**: custom > ogrdb > vdjbase > imgt

**Changes**:
1. Fix `DEFAULT_PROVIDERS` to `["custom", "ogrdb", "vdjbase", "imgt"]`
2. Add test verifying priority order

**Test**: `tests/unit/germlines/test_compliance.py::test_default_priority_order`

---

### Group B: Single-Provider Enforcement (T053, T054)

**Files**: 
- `src/sadie/germlines/manager.py`
- `tests/unit/germlines/test_compliance.py`

**Changes**:
1. `GermlineManager` validates that all segment queries use same provider config
2. Document that per-segment provider mixing is not supported
3. Add tests that verify rejection of per-segment provider params

**Implementation**:
- Manager already enforces this by design (one provider list for all queries)
- Need tests confirming this behavior and documenting in docstrings

---

### Group C: Species Error Handling (T055)

**Files**:
- `src/sadie/germlines/manager.py`
- `src/sadie/airr/igblast/germline.py`

**Changes**:
1. `GermlineManager.get_genes()` - raise clear error when no provider has species data (not just return empty)
2. `GermlineManager._fetch_from_provider()` - distinguish "no data" from "error"
3. `GermlineData.__init__()` - raise `ValueError` instead of warning + fallback when species missing

**Test**: `tests/unit/germlines/test_compliance.py::test_clear_error_missing_species`

---

### Group D: No-Fallback Enforcement (T059)

**Files**:
- `src/sadie/airr/igblast/germline.py`

**Changes**:
1. Remove silent fallback to legacy G3 paths in `GermlineData.__init__()`
2. When `SADIE_USE_GERMLINES_MODULE=true` and species not found, raise error with instructions
3. Add test verifying no fallback occurs

**Test**: `tests/unit/germlines/test_compliance.py::test_no_g3_fallback`

---

### Group E: Custom Ingestion Validation (T056)

**Files**:
- `src/sadie/germlines/providers/custom.py`
- `tests/unit/germlines/test_compliance.py`

**Current**: `_validate_sequence()` checks nucleotides only

**Changes**:
1. Add FASTA header validation (detect malformed headers)
2. Add AA sequence validation if AA sequences provided
3. Return detailed error messages identifying specific failures
4. Add test with malformed FASTA files

**Test**: `tests/unit/germlines/test_compliance.py::test_custom_ingestion_validation`

---

### Group F: Species/Chain/Segment Parity (T057)

**Files**:
- `tests/unit/germlines/test_compliance.py`

**Changes**:
1. Create test that queries all species/chain/segment combinations
2. Compare coverage against expected support (human, mouse, rat, rabbit, macaque, dog)
3. Verify H, K, L chains supported
4. Verify V, D, J segments supported

**Test**: `tests/unit/germlines/test_compliance.py::test_species_chain_segment_parity`

---

### Group G: Gapped Sequence Fail-Fast (T060)

**Files**:
- `src/sadie/germlines/renumbering_integration.py`

**Current**: `_get_vj_alignment_pairs()` silently skips genes without gapped sequences

**Changes**:
1. Track genes that have neither `sequence_aa_gapped` nor `sequence_gapped`
2. Raise `ValueError` listing specific genes lacking gapped data
3. Add test verifying fail-fast behavior

**Test**: `tests/unit/germlines/test_compliance.py::test_gapped_sequence_fail_fast`

---

## Execution Order

```
1. Group A (T058) - Priority fix (quick, foundational)
2. Group B (T053-T054) - Single-provider (tests existing behavior)
3. Group C (T055) - Species error handling
4. Group D (T059) - No-fallback enforcement
5. Group E (T056) - Custom validation
6. Group F (T057) - Parity verification
7. Group G (T060) - Gapped fail-fast
```

## New Test File Structure

```python
# tests/unit/germlines/test_compliance.py

class TestPriorityOrder:
    def test_default_priority_order(self): ...  # T058

class TestSingleProvider:
    def test_single_provider_enforcement(self): ...  # T053
    def test_rejects_per_segment_providers(self): ...  # T054

class TestErrorHandling:
    def test_clear_error_missing_species(self): ...  # T055
    def test_no_g3_fallback(self): ...  # T059

class TestCustomValidation:
    def test_rejects_malformed_fasta(self): ...  # T056
    def test_rejects_invalid_amino_acids(self): ...  # T056
    def test_detailed_error_messages(self): ...  # T056

class TestParity:
    def test_species_chain_segment_parity(self): ...  # T057

class TestGappedSequences:
    def test_fail_fast_missing_gapped(self): ...  # T060
```

## Success Criteria

- [ ] Default priority order is custom > ogrdb > vdjbase > imgt
- [ ] Single-provider enforcement documented and tested
- [ ] Clear error when provider lacks species (not empty result)
- [ ] No silent fallback to G3 when germlines selected
- [ ] Custom ingestion rejects malformed FASTA with detailed errors
- [ ] Germlines supports same species/chains/segments as G3
- [ ] Fail-fast when gapped sequences missing for HMM building
- [ ] All new tests pass
- [ ] Existing tests still pass

## Risk Assessment

| Risk | Mitigation |
|------|------------|
| Breaking existing workflows | Feature flag already protects legacy behavior |
| Missing species data | Clear error messages with setup instructions |
| Test flakiness | Use deterministic test data |

---
*Created: 2026-01-21*
