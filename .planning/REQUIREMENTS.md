# Requirements: Germline Database Integration

**Defined:** 2026-01-22
**Core Value:** Enable researchers to select germline database for AIRR annotation and renumbering

## v1.2 Requirements

Requirements for Reference Module Unification milestone.

### Source Validation

- [x] **SRC-01**: Expand VALID_SOURCES to include `ogrdb`, `vdjbase`
- [x] **SRC-02**: Validate source exists in germlines before processing

### Integration

- [x] **INT-01**: Add `use_germlines=True` parameter to `References.from_yaml()`
- [x] **INT-02**: Route source selection through GermlineManager (explicit source, no priority)
- [x] **INT-03**: Generate synthetic `_id` field in adapter

### Build CLI

- [x] **CLI-01**: Add `sadie reference build <yaml> --output <path>` command
- [x] **CLI-02**: Build generates complete IgBLAST database structure
- [x] **CLI-03**: Progress output during build

### Runtime Usage

- [x] **RUN-01**: Add `Airr(database=<path>)` parameter to use prebuilt database
- [x] **RUN-02**: Skip germlines/G3 lookup when using prebuilt database
- [x] **RUN-03**: Validate database structure on load

### Documentation

- [x] **DOC-01**: Create reference-sample.yml (mouse=imgt, human=ogrdb, macaque=vdjbase)
- [x] **DOC-02**: Document build → use workflow

---

## v1.1 Requirements (Complete)

Requirements for audit validation milestone.

### Phase 13: Audit

- [x] **AUDIT-01**: Run AIRR annotation with germlines backend on test sequences
- [x] **AUDIT-02**: Run AIRR annotation with G3 backend on same sequences
- [x] **AUDIT-03**: Compare results for column-level identity (excluding source column)
- [x] **AUDIT-04**: Document any discrepancies with root cause analysis

### Phase 14: C Region Integration

- [x] **CREG-01**: Update germlines sources to pull C region data from IMGT/OGRDB/VDJbase
- [x] **CREG-02**: Generate IgBLAST C gene databases in germlines module
- [x] **CREG-03**: Verify C gene columns present in AIRR output
- [x] **CREG-04**: Re-run audit to validate parity improvement

### Phase 15: J Gene Matching

- [x] **JFIX-01**: Investigate IgBLAST J gene database configuration
- [x] **JFIX-02**: Verify aux file format and content (fixed: 5-column format)
- [x] **JFIX-03**: Check internal_data directory structure
- [x] **JFIX-04**: Debug IgBLAST execution and parameters
- [x] **JFIX-05**: Re-run audit to validate CDR3 annotation

### Phase 16: NDM.IMGT FWR3 Fix

- [x] **NDM-01**: Fix build_internal_data.py to calculate correct FWR3 end position
- [x] **NDM-02**: Regenerate ndm.imgt files for human
- [x] **NDM-03**: Re-run audit to validate parity improvement

### Phase 17: complete_vdj Fix

- [x] **VDJ-01**: Investigate post-processing solution (AIRR-standard recalculation)
- [x] **VDJ-03**: Verify complete_vdj is accurate per AIRR standard
- [x] **VDJ-04**: Document the IgBLAST quirk in audit/igblast-quirk.md

## Out of Scope

| Feature | Reason |
|---------|--------|
| Renumbering parity audit | Focus on AIRR first, renumbering already tested in v1.0 |
| Performance benchmarking | Parity focus, not speed |
| Multi-species audit | Human only for initial validation |
| Multi-provider blending per analysis | Out of scope for v1.2 |
| Automatic G3 fallback | Explicit source only |
| Real-time provider sync | Offline databases |
| Per-gene source overrides | Reference-level source only |
| TR (T-cell receptor) support | IG focus for v1.2 |

## Traceability

### v1.2 Requirements Mapping

| Requirement | Phase | Status |
|-------------|-------|--------|
| SRC-01 | Phase 19 | Complete |
| SRC-02 | Phase 24 | Complete |
| INT-01 | Phase 20 | Complete |
| INT-02 | Phase 20 | Complete |
| INT-03 | Phase 20 | Complete |
| CLI-01 | Phase 21 | Complete |
| CLI-02 | Phase 21 | Complete |
| CLI-03 | Phase 21 | Complete |
| RUN-01 | Phase 22 | Complete |
| RUN-02 | Phase 22 | Complete |
| RUN-03 | Phase 22 | Complete |
| DOC-01 | Phase 23 | Complete |
| DOC-02 | Phase 23 | Complete |

**v1.2 Coverage:**
- Total requirements: 12
- Mapped to phases: 12
- Complete: 12 ✓
- Unmapped: 0 ✓

### v1.1 Requirements Mapping

| Requirement | Phase | Status |
|-------------|-------|--------|
| AUDIT-01 | Phase 13 | Complete |
| AUDIT-02 | Phase 13 | Complete |
| AUDIT-03 | Phase 13 | Complete |
| AUDIT-04 | Phase 13 | Complete |
| CREG-01 | Phase 14 | Complete |
| CREG-02 | Phase 14 | Complete |
| CREG-03 | Phase 14 | Complete |
| CREG-04 | Phase 14 | Complete |
| JFIX-01 | Phase 15 | Complete |
| JFIX-02 | Phase 15 | Complete |
| JFIX-03 | Phase 15 | Complete |
| JFIX-04 | Phase 15 | Complete |
| JFIX-05 | Phase 15 | Complete |
| NDM-01 | Phase 16 | Complete |
| NDM-02 | Phase 16 | Complete |
| NDM-03 | Phase 16 | Complete |
| VDJ-01 | Phase 17 | Complete |
| VDJ-03 | Phase 17 | Complete |
| VDJ-04 | Phase 17 | Complete |

**v1.1 Coverage:**
- Total requirements: 19
- Mapped to phases: 19
- Complete: 19 ✓
- Unmapped: 0 ✓

---

## v1.3 Requirements

Requirements for Test Infrastructure & Species Expansion milestone.

### Phase 25: Macaque Germlines Integration

- [x] **MAC-01**: Build macaque IgBLAST databases in germlines module
- [x] **MAC-02**: Generate macaque internal_data and aux files
- [x] **MAC-03**: Verify macaque AIRR annotation works
- [x] **MAC-04**: Remove skip markers from macaque tests

### Phase 26: AIRR Package Dependency

- [x] **AIRR-01**: Add `airr` package to pyproject.toml dependencies
- [x] **AIRR-02**: Verify airr package installs correctly
- [x] **AIRR-03**: Remove importorskip from test

### Phase 27: Remove Deprecated G3 Tests

- [ ] **G3-01**: Review what these tests are validating
- [ ] **G3-02**: Determine if equivalent germlines module tests exist
- [ ] **G3-03**: Either migrate tests to germlines or remove if redundant
- [ ] **G3-04**: Update deprecation timeline documentation

### Phase 28: Fix Germline Priority Order

- [ ] **PRIO-01**: Update default provider priority in GermlineManager
- [ ] **PRIO-02**: Document priority rationale in code comments
- [ ] **PRIO-03**: Verify priority order used in fallback resolution
- [ ] **PRIO-04**: Test priority order with multi-source queries

### v1.3 Requirements Mapping

| Requirement | Phase | Status |
|-------------|-------|--------|
| MAC-01 | Phase 25 | Complete |
| MAC-02 | Phase 25 | Complete |
| MAC-03 | Phase 25 | Complete |
| MAC-04 | Phase 25 | Complete |
| AIRR-01 | Phase 26 | Complete |
| AIRR-02 | Phase 26 | Complete |
| AIRR-03 | Phase 26 | Complete |
| G3-01 | Phase 27 | Pending |
| G3-02 | Phase 27 | Pending |
| G3-03 | Phase 27 | Pending |
| G3-04 | Phase 27 | Pending |
| PRIO-01 | Phase 28 | Pending |
| PRIO-02 | Phase 28 | Pending |
| PRIO-03 | Phase 28 | Pending |
| PRIO-04 | Phase 28 | Pending |

**v1.3 Coverage:**
- Total requirements: 15
- Mapped to phases: 15
- Complete: 7 (Phases 25-26)
- Pending: 8 (Phases 27-28)

---
*Requirements defined: 2026-01-22*
*Last updated: 2026-01-25 — v1.3 Phase 26 requirements complete*
