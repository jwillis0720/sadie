# Requirements: Germline Database Integration

## v1 Requirements

### Provider Selection

- [x] **PROV-01**: System MUST allow users to specify a germline provider (imgt, ogrdb, vdjbase, custom) when initializing AIRR analysis (FR-001)
- [x] **PROV-02**: System MUST allow users to specify a germline provider when initializing Renumbering operations (FR-002)
- [x] **PROV-03**: System MUST expose a `germline_backend` parameter accepting "g3" (default) or "germlines" (FR-003)
- [ ] **PROV-04**: System MUST use default provider priority (custom > ogrdb > vdjbase > imgt) when no provider specified (FR-004)
- [x] **PROV-05**: System MUST use local germlines module data when germlines backend selected (FR-005)
- [ ] **PROV-06**: System MUST apply single provider selection to all V/D/J segments within a run (FR-014)

### Error Handling

- [ ] **ERR-01**: System MUST provide clear error messages when provider has no data for specified species (FR-006)
- [ ] **ERR-02**: System MUST validate custom germlines at ingestion, rejecting invalid sequences with detailed errors (FR-012)
- [ ] **ERR-03**: System MUST ensure gapped AA sequences available for all V/J genes in HMM building; fail with clear error if missing (FR-013)
- [ ] **ERR-04**: System MUST NOT silently fall back to G3 when germlines backend fails (NFR-002)

### Backwards Compatibility

- [x] **COMPAT-01**: G3 remains default backend; germlines is opt-in via explicit parameter (FR-007)
- [x] **COMPAT-02**: Existing tests in tests/unit/airr/ continue to pass
- [x] **COMPAT-03**: Existing tests in tests/unit/renumbering/ continue to pass
- [x] **COMPAT-04**: Output format/schema identical regardless of provider (FR-011)

### Testing

- [x] **TEST-01**: New test directory tests/unit/germlines/ with mirrored AIRR tests (FR-008)
- [x] **TEST-02**: Mirrored renumbering tests using germlines backend (FR-009)
- [ ] **TEST-03**: Tests verify same species/chains/segments supported as existing modules (FR-010)
- [ ] **TEST-04**: Gapped AA fallback translation test when only gapped NT available

### Performance

- [x] **PERF-01**: Germline lookup performance equivalent to G3 backend (NFR-001)

## v2 Requirements (Deferred)

- Multi-provider blending per analysis (currently single provider per run)
- GUI for provider selection
- Provider-specific analytics/reporting

## Out of Scope

- T-cell receptor (TR) germlines — focus on immunoglobulin (IG) only
- Real-time provider synchronization — manual update via update_databases()
- Provider switching mid-linked-analysis consistency — user responsibility

## Traceability

| Requirement | Phase | Tasks |
|-------------|-------|-------|
| PROV-01 | Phase 3 (US1: AIRR) | T010-T015 |
| PROV-02 | Phase 4 (US2: Renumbering) | T016-T023 |
| PROV-03 | Phase 2 (Foundational) | T005-T009 |
| PROV-04 | Phase 9 (Compliance) | T058 |
| PROV-05 | Phase 3, 4 | T010-T023 |
| PROV-06 | Phase 9 (Compliance) | T053-T054 |
| ERR-01 | Phase 9 (Compliance) | T055 |
| ERR-02 | Phase 9 (Compliance) | T056 |
| ERR-03 | Phase 9 (Compliance) | T004a, T035a, T060 |
| ERR-04 | Phase 9 (Compliance) | T055, T059 |
| COMPAT-01 | Phase 2 (Foundational) | T009, T012 |
| COMPAT-02 | Phase 8 (Polish) | T052 |
| COMPAT-03 | Phase 8 (Polish) | T052 |
| COMPAT-04 | Phase 5 (Reference) | T030 |
| TEST-01 | Phase 6 (US3: Tests) | T031-T035 |
| TEST-02 | Phase 6 (US3: Tests) | T036-T038 |
| TEST-03 | Phase 9 (Compliance) | T057 |
| TEST-04 | Phase 6 (US3: Tests) | T035a |
| PERF-01 | Phase 8 (Polish) | T049-T050 |

---
*Last updated: 2026-01-21 — converted from spec-kit FR/NFR requirements*
