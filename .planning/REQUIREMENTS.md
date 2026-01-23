# Requirements: Germline Database Integration

**Defined:** 2026-01-22
**Core Value:** Enable researchers to select germline database for AIRR annotation and renumbering

## v1.1 Requirements

Requirements for audit validation milestone.

### Audit

- [x] **AUDIT-01**: Run AIRR annotation with germlines backend on test sequences
- [x] **AUDIT-02**: Run AIRR annotation with G3 backend on same sequences
- [x] **AUDIT-03**: Compare results for column-level identity (excluding source column)
- [x] **AUDIT-04**: Document any discrepancies with root cause analysis

### C Region Integration

- [ ] **CREG-01**: Update germlines sources to pull C region data from IMGT/OGRDB/VDJbase
- [ ] **CREG-02**: Generate IgBLAST C gene databases in germlines module
- [ ] **CREG-03**: Verify C gene columns present in AIRR output
- [ ] **CREG-04**: Re-run audit to validate parity improvement

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
| CREG-01 | Phase 14 | Pending |
| CREG-02 | Phase 14 | Pending |
| CREG-03 | Phase 14 | Pending |
| CREG-04 | Phase 14 | Pending |

**Coverage:**
- v1.1 requirements: 8 total
- Mapped to phases: 8
- Unmapped: 0 ✓

---
*Requirements defined: 2026-01-22*
*Last updated: 2026-01-22 after v1.1 milestone start*
