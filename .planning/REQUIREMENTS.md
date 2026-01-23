# Requirements: Germline Database Integration

**Defined:** 2026-01-22
**Core Value:** Enable researchers to select germline database for AIRR annotation and renumbering

## v1.1 Requirements

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

## Traceability

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

**Coverage:**
- v1.1 requirements: 19 total
- Mapped to phases: 19
- Complete: 19 ✓
- Unmapped: 0 ✓

---
*Requirements defined: 2026-01-22*
*Last updated: 2026-01-22 — v1.1 milestone complete*
