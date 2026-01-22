# Tasks: Germlines Module Completion

**Feature Branch**: `001-germline-completion`
**Generated**: 2026-01-08
**Spec**: [spec.md](./spec.md) | **Plan**: [plan.md](./plan.md)

## Overview

This document breaks down the implementation into actionable tasks organized by user story. Each user story phase is independently testable and delivers incremental value.

**Total Tasks**: 60 (plus 2 backlog)
**Parallelizable Tasks**: 29
**User Stories**: 6 (4 × P1, 2 × P2)
**Estimated Duration**: 24-28 hours

## Task Format

```
- [ ] T### [P] [US#] Description with file path
```

- `T###`: Sequential task ID
- `[P]`: Parallelizable (optional marker)
- `[US#]`: User story label (US1-US6)
- Description: Clear action with exact file path

## Implementation Strategy

**MVP Scope**: User Story 1 (Add Custom Germline Sequences) + User Story 6 (Populate Reference Data)
- Delivers core value: custom germline injection
- Establishes data population workflow
- Estimated: 8-10 hours

**Incremental Delivery**:
1. MVP (US1 + US6): Custom sequences + data population
2. Phase 2 (US2 + US4): Offline operation + backward compatibility
3. Phase 3 (US3 + US5): Priority system + VDJbase provider

## Dependencies Graph

```
Phase 1 (Setup)
    ↓
Phase 2 (Foundational)
    ↓
├── US1 (P1) ← Independent
├── US2 (P1) ← Requires US1 (needs populated data)
├── US3 (P2) ← Independent
├── US4 (P1) ← Requires US1 (needs working module)
├── US5 (P2) ← Independent
└── US6 (P1) ← Independent
    ↓
Phase 8 (Polish)
```

**Parallel Execution Opportunities**:
- US1, US3, US5, US6 can be developed in parallel
- US2, US4 require US1 completion
- Within each story: Data model + contracts + services can be parallel if tests not blocking

---

## Phase 1: Project Setup

**Goal**: Initialize project infrastructure and tooling

**Duration**: 1-2 hours

### Setup Tasks

- [X] T001 Create VDJbase provider directory structure in src/sadie/germlines/sources/vdjbase/
- [X] T002 Create VDJbase provider human subdirectory in src/sadie/germlines/sources/vdjbase/human/
- [X] T003 [P] Create test data directory structure in src/sadie/germlines/tests/data/{provider}/{species}/ (custom/, imgt/, ogrdb/, vdjbase/ with human/)
- [X] T004 [P] Create curated test dataset with 5-10 genes per segment in src/sadie/germlines/tests/data/{provider}/human/
- [X] T004a [P] Create FR-025b regression test sequences (IGHV1-69*01, IGHV3-23*01, IGHD3-3*01, IGHJ4*01) in src/sadie/germlines/tests/data/regression/
- [X] T004b [P] Create G3 regression test comparing germlines module output vs expected G3 format for FR-025b sequences in src/sadie/germlines/tests/test_g3_regression.py
- [X] T005 Create validation script template in src/sadie/germlines/scripts/validate.py
- [X] T006 Update pyproject.toml to ensure all dependencies listed (BioPython, pytest)
- [X] T007 Create feature flag utility module in src/sadie/germlines/utils/feature_flags.py
- [X] T096 Create research.md (Phase 0 deliverable) in specs/001-germline-completion/research.md per plan
- [X] T097 Create data-model.md (Phase 1 deliverable) in specs/001-germline-completion/data-model.md per plan
- [X] T098 Create quickstart.md (Phase 1 deliverable) in specs/001-germline-completion/quickstart.md per plan

---

## Phase 2: Foundational Components

**Goal**: Implement shared infrastructure required by all user stories

**Duration**: 3-4 hours

**Prerequisites**: Phase 1 complete

### Foundational Tasks

- [X] T008 Implement feature flag function `use_germlines_module()` in src/sadie/germlines/utils/feature_flags.py
- [X] T009 [P] Implement auto-gapping service using BioPython alignment against IMGT-gapped templates (per-gene fallback to per-segment consensus) in src/sadie/germlines/builders/gapper.py
- [X] T010 [P] Add logging configuration for germlines module in src/sadie/germlines/__init__.py
- [X] T011 [P] Update GermlineManager to support vdjbase provider in src/sadie/germlines/manager.py
- [X] T012 Create VDJbase provider stub with base interface in src/sadie/germlines/providers/vdjbase.py
- [X] T013 Implement VDJbase FASTA parsing logic in src/sadie/germlines/providers/vdjbase.py
- [X] T014 Implement VDJbase provider metadata methods in src/sadie/germlines/providers/vdjbase.py
- [X] T015 Add timing metrics logging to pipeline.py in src/sadie/germlines/pipeline.py
- [X] T091 [P] Ensure normalized outputs follow FR-022/022a/023 (normalized/{species}/gapped|ungapped; D segments ungapped-only) and add tests validating paths/content
- [X] T094 [P] Gap amino-acid then back-map to nucleotide per FR-021b in src/sadie/germlines/builders/gapper.py; add acceptance test for codon-aware gaps
- [X] T095 [P] Reuse existing pipeline utilities (no duplicated gapping/build logic) and add tests validating normalized layout per FR-024

---

## Phase 3: User Story 1 - Add Custom Germline Sequences (P1)

**Goal**: Enable researchers to add novel IGHV alleles to local database for immediate use

**Independent Test**: User adds FASTA to `src/sadie/germlines/sources/custom/human/IGHV.fasta`, runs Sadie, sees custom allele in results

**Duration**: 2-3 hours

**Prerequisites**: Phase 2 complete

### US1 Tasks

- [X] T016 [US1] Verify custom provider handles new sequences in src/sadie/germlines/providers/custom.py (integrated GapperService for auto-gapping)
- [X] T017 [US1] Implement change detection for custom sequences in src/sadie/germlines/pipeline.py
- [X] T018 [US1] Add validation for custom FASTA files (nucleotides, format) in src/sadie/germlines/providers/custom.py
- [X] T019 [US1] Implement auto-rebuild trigger on custom file change in src/sadie/germlines/pipeline.py
- [X] T020 [P] [US1] Write unit test for custom sequence priority in src/sadie/germlines/tests/test_custom_provider.py
- [X] T021 [P] [US1] Write integration test for custom sequence end-to-end in src/sadie/germlines/tests/test_integration.py
- [X] T022 [US1] Add logging for custom sequence load events in src/sadie/germlines/providers/custom.py
- [X] T023 [US1] Document custom sequence addition process in src/sadie/germlines/sources/custom/README.md

**Acceptance Criteria**:
- [X] User can add novel IGHV to src/sadie/germlines/sources/custom/human/IGHV.fasta
- [X] Pipeline auto-detects change and rebuilds (<5 minutes)
- [X] Custom version takes priority over IMGT per Constitution Principle II
- [X] Invalid sequences log warning but continue with valid ones

---

## Phase 4: User Story 6 - Populate Reference Data Sources (P1)

**Goal**: Enable new users to set up germlines module with standard IMGT and OGRDB data

**Independent Test**: Fresh install, user follows README, downloads IMGT/OGRDB, runs validation, receives confirmation

**Duration**: 3-4 hours

**Prerequisites**: Phase 2 complete

### US6 Tasks

- [X] T024 [US6] Complete IMGT download script implementation in src/sadie/germlines/scripts/download_imgt.py
- [X] T025 [US6] Add species parameter support to IMGT download script in src/sadie/germlines/scripts/download_imgt.py
- [X] T026 [US6] Implement resume capability for IMGT downloads in src/sadie/germlines/scripts/download_imgt.py
- [X] T027 [US6] Implement OGRDB download script in src/sadie/germlines/scripts/download_ogrdb.py
- [X] T028 [US6] Add species parameter support to OGRDB download script in src/sadie/germlines/scripts/download_ogrdb.py
- [X] T029 [US6] Implement validation for downloaded FASTA files in src/sadie/germlines/scripts/validate.py
- [X] T030 [P] [US6] Create VDJbase manual download instructions in src/sadie/germlines/sources/vdjbase/README.md
- [X] T031 [P] [US6] Update IMGT data documentation in src/sadie/germlines/sources/imgt/IMGT_DATA.md
- [X] T032 [P] [US6] Update OGRDB data documentation in src/sadie/germlines/sources/ogrdb/OGRDB_DATA.md
- [X] T033 [US6] Add progress indicators to download scripts (INFO logging) in src/sadie/germlines/scripts/download_imgt.py
- [X] T034 [US6] Add timing metrics to download scripts in src/sadie/germlines/scripts/download_ogrdb.py
- [X] T090 [US6] Implement resume capability for OGRDB downloads in src/sadie/germlines/scripts/download_ogrdb.py (checkpoint file; idempotent retries)
- [X] T093 [P] [US6] Standardize progress logging cadence to INFO every 10 files in IMGT/OGRDB download scripts ("Downloaded {completed}/{total} files ({percentage}%)")

**Acceptance Criteria**:
- [X] User runs `python src/sadie/germlines/scripts/download_imgt.py human` and gets validated FASTA files
- [X] Pipeline builds BLAST databases automatically (~1-2 min) with timing logs
- [X] Download script resumes from checkpoint if interrupted
- [X] Validation script confirms data ready with clear success message

---

## Phase 5: User Story 2 - Use Local Germline Databases Offline (P1)

**Goal**: Enable bioinformaticians in air-gapped environments to annotate sequences without internet

**Independent Test**: Populate data once, disconnect network, run full Sadie pipeline successfully

**Duration**: 2-3 hours

**Prerequisites**: Phase 3 (US1), Phase 4 (US6) complete

### US2 Tasks

- [X] T035 [US2] Verify pipeline operates without network calls in src/sadie/germlines/pipeline.py
- [X] T036 [US2] Add offline mode detection and logging in src/sadie/germlines/pipeline.py
- [X] T037 [US2] Implement clear error messages for missing data in src/sadie/germlines/providers/base.py
- [X] T038 [US2] Add README references in error messages for data population in src/sadie/germlines/manager.py
- [X] T039 [P] [US2] Write offline integration test (network disabled) in src/sadie/germlines/tests/test_offline_operation.py
- [X] T040 [US2] Verify cached data usage (6 months old) works correctly in src/sadie/germlines/pipeline.py

**Acceptance Criteria**:
- [X] IMGT and OGRDB data populated in src/sadie/germlines/sources/
- [X] Sadie annotation completes without internet access
- [X] Empty src/sadie/germlines/sources/ gives clear error with remediation steps
- [X] Pipeline uses cached data without update checks

---

## Phase 6: User Story 4 - Integrate with Existing Sadie Workflows (P1)

**Goal**: Ensure existing Sadie users' workflows work transparently with new germlines module

**Independent Test**: Run existing Sadie IgBLAST test suite, all tests pass without modification

**Duration**: 4-5 hours

**Prerequisites**: Phase 3 (US1) complete

### US4 Tasks

- [X] T041 [US4] Update IgBLAST germline paths to new module structure in src/sadie/airr/igblast/germline.py
- [X] T042 [US4] Add feature flag check to IgBLAST integration in src/sadie/airr/igblast/germline.py
- [X] T043 [US4] Update Reference system to query germlines module in src/sadie/reference/reference.py
- [X] T044 [US4] Implement G3 API response format adapter (regions fields) in src/sadie/reference/reference.py
- [X] T045 [US4] Add feature flag check to Reference system in src/sadie/reference/reference.py
- [X] T046 [US4] Update HMM builder to use germlines module in src/sadie/renumbering/aligners/hmmer.py
- [X] T047 [US4] Add feature flag check to HMM builder in src/sadie/renumbering/aligners/hmmer.py
- [X] T048 [US4] Verify gapped sequences available for Stockholm alignment in src/sadie/germlines/builders/gapper.py
- [X] T049 [P] [US4] Run existing Sadie IgBLAST test suite (regression test)
- [X] T050 [P] [US4] Test feature flag SADIE_USE_GERMLINES_MODULE=false (G3 mode)
- [X] T051 [US4] Document backward compatibility approach in src/sadie/germlines/INTEGRATION_GUIDE.md

**Acceptance Criteria**:
- [X] `GermlineData("human")` resolves to new database locations
- [X] Reference system returns G3-compatible format
- [X] HMM builder gets gapped sequences successfully
- [X] Feature flag=false falls back to G3 without errors

---

## Phase 7A: User Story 3 - Priority-Based Database Selection (P2)

**Goal**: Allow researchers to use OGRDB preferentially over IMGT with fallback

**Independent Test**: Configure `["ogrdb", "imgt"]`, run annotation, verify OGRDB used when present, IMGT for missing genes

**Duration**: 1-2 hours

**Prerequisites**: Phase 2 complete (independent from US1)

### US3 Tasks

- [X] T052 [P] [US3] Verify priority ordering logic in GermlineManager in src/sadie/germlines/manager.py
- [X] T053 [P] [US3] Write unit test for priority ordering scenarios in src/sadie/germlines/tests/test_priority_ordering.py
- [X] T054 [P] [US3] Test deduplication rules (same name, same sequence, novel) in src/sadie/germlines/tests/test_priority_ordering.py
- [X] T055 [US3] Document priority configuration in src/sadie/germlines/README.md
- [X] T056 [US3] Add logging for priority-based gene selection in src/sadie/germlines/manager.py

**Acceptance Criteria**:
- [X] Both OGRDB and IMGT have IGHV1-69*01, OGRDB version used with priority `["ogrdb", "imgt"]`
- [X] Gene in IMGT but not OGRDB is included in merged database
- [X] Two providers with same sequence but different names keeps both

---

## Phase 7B: User Story 5 - Add VDJbase Provider (P2)

**Goal**: Enable researchers to use VDJbase genotype data for population-specific analysis

**Independent Test**: Populate VDJbase sample data, configure `["vdjbase", "imgt"]`, verify VDJbase sequences used

**Duration**: 2-3 hours

**Prerequisites**: Phase 2 complete (independent from US1)

### US5 Tasks

- [X] T057 [P] [US5] Complete VDJbase provider implementation in src/sadie/germlines/providers/vdjbase.py
- [X] T058 [P] [US5] Add VDJbase to default provider list in src/sadie/germlines/manager.py (verified: manager.py DEFAULT_PROVIDERS includes vdjbase, __init__.py exports VDJbaseProvider, data exists for human/rhesus_macaque)
- [X] T059 [P] [US5] Write unit tests for VDJbase provider in src/sadie/germlines/tests/test_vdjbase_provider.py
- [X] T060 [P] [US5] Create VDJbase test data in src/sadie/germlines/tests/data/vdjbase/human/
- [X] T061 [US5] Test VDJbase in priority ordering in src/sadie/germlines/tests/test_priority_ordering.py
- [X] T062 [US5] Add VDJbase error handling for format changes in src/sadie/germlines/providers/vdjbase.py

*Note: T063 removed - duplicates T030 (VDJbase README created in Phase 4)*

**Acceptance Criteria**:
- [X] VDJbase FASTA files in src/sadie/germlines/sources/vdjbase/human/ are parsed successfully
- [X] VDJbase allele used instead of IMGT when priority includes vdjbase first
- [X] Format errors show clear message with documentation link

---

## Phase 8: Polish & Cross-Cutting Concerns

**Goal**: Finalize integration, testing, and documentation

**Duration**: 3-4 hours

**Prerequisites**: All user story phases complete

### Polish Tasks

- [ ] T064 Run full Sadie test suite and ensure 100% pass rate
- [ ] T065 Verify SC-001: Custom germline injection in <5 minutes (stopwatch test)
- [ ] T066 Verify SC-002: Offline operation works (network disabled test)
- [ ] T067 Verify SC-003: Download scripts <10 minutes for human (timing test)
- [ ] T068 Verify SC-005: Disk usage <500MB for human data (du -sh test)
- [ ] T069 Verify SC-009: Change detection triggers rebuild correctly (file modification test)
- [ ] T070 Verify SC-010: Clear error messages for all failure modes (error handling review)
- [ ] T071 Verify SC-011: Timing metrics logged for major operations (log inspection)
- [ ] T072 Verify SC-012: CI tests complete in <5 minutes (GitHub Actions check)
- [X] T073 Update main germlines README with completion status in src/sadie/germlines/README.md
- [X] T074 Update INTEGRATION_GUIDE with actual integration points in src/sadie/germlines/INTEGRATION_GUIDE.md
- [X] T075 Add performance profiling for critical paths in src/sadie/germlines/pipeline.py
- [ ] T076 Run pre-commit hooks and fix any linting issues
- [ ] T077 Update CHANGELOG or release notes with germlines module completion
- [ ] T099 Benchmark rebuild time (<2 minutes for human dataset) and record timing logs to satisfy SC-001/SC-011 in specs/001-germline-completion/validation-tracking.md

### Supplemental Tasks (Added 2026-01-15)

The following tasks were added after initial task generation to address IgBLAST auxiliary requirements (FR-037-039) and logging standardization (FR-032, FR-035, FR-036).

- [X] T078 [P] Complete AuxFileBuilder CDR/FWR boundary detection and IgBLAST format output in src/sadie/germlines/builders/aux.py
- [X] T079 [P] Update BlastDBBuilder to use -hash_index and validate output naming/paths per FR-038/FR-038a in src/sadie/germlines/builders/blast.py
- [X] T080 [P] Generate igblast/internal_data/organism.yaml per FR-039/FR-039a in src/sadie/germlines/pipeline.py
- [X] T081 [P] Enforce structured logging format and key-value fields per FR-032a/FR-032b in src/sadie/germlines/__init__.py
- [X] T082 [P] Add change-detection log details (path/hash/change type) per FR-035a in src/sadie/germlines/pipeline.py
- [X] T083 [P] Standardize error message templates per FR-036a/FR-036b in src/sadie/germlines/providers/*.py and src/sadie/germlines/builders/blast.py
- [X] T084 [P] Add legacy GermlineData API compatibility tests in src/sadie/germlines/tests/test_germline_data_legacy.py
- [X] T085 [P] Add regression tests comparing germlines vs G3 output in src/sadie/germlines/tests/test_g3_regression.py
- [X] T086 [P] Replace .specify/memory/constitution.md placeholder with actual constitution text
- [X] T087 [US4] Create validation period tracking document at specs/001-germline-completion/validation-tracking.md with: start date, release count, bug tracker, performance baseline comparison template
- [X] T088 [US4] Add deprecation warning log when SADIE_USE_GERMLINES_MODULE=false: "G3 API is deprecated. Set SADIE_USE_GERMLINES_MODULE=true. G3 will be removed after {date}."
- [X] T089 [P] [US6] Verify OGRDB provider correctly loads _gapped.fasta files alongside ungapped files in src/sadie/germlines/providers/ogrdb.py
- [X] T092 [US4] Enforce validation period per FR-017a/b: track start/end, release counts, success criteria, and deprecation notice schedule; log metrics and update validation-tracking.md
- [X] T100 [P] Consolidate G3 regression parity into a single test at src/sadie/germlines/tests/test_g3_regression.py (supersedes T004b and T085)
- [X] T101 [P] Implement HMM builder for Stockholm alignment generation per FR-013a-c in src/sadie/germlines/builders/hmm.py

### Backlog (Post-Validation Period)

These tasks execute after validation period success criteria (FR-017b) are met:

- [ ] T-BACKLOG-001 Issue deprecation notice per FR-019a (CHANGELOG entry, GitHub discussion, runtime warning active)
- [ ] T-BACKLOG-002 Remove G3 dependencies per FR-019b (src/sadie/renumbering/clients/g3.py, all G3 imports, feature flag code, G3-related tests)

---

## Parallel Execution Plan

### Within Phases (Independent Work)

**Phase 1 Setup**:
- T003, T004 (test data) parallel with T001, T002 (directories)

**Phase 2 Foundational**:
- T009 (gapper), T010 (logging), T011 (manager update), T015 (timing) - all parallel
- New: T091 (normalized outputs), T094 (AA→codon gapping), T095 (utility reuse) can run in parallel with T009/T015

**Phase 3 US1**:
- T020, T021 (tests) parallel with T023 (docs)

**Phase 4 US6**:
- T030, T031, T032 (documentation) - all parallel
- New: T090 (OGRDB resume) and T093 (logging cadence) parallel with T026 (IMGT resume)

**Phase 7 (Independent Stories)**:
- US3 (T052-T056) and US5 (T057-T063) can be done in parallel

### Across Phases (Team Distribution)

```
Developer 1: US1 → US2 → US4
Developer 2: US6 → US3 → US5
 Developer 3: Phase 1 → Phase 2 → Phase 8 (support + polish; owns T096-T098 deliverables and T099 benchmark)
```

---

## Testing Strategy

### Unit Tests

- Test each provider in isolation with curated test dataset
- Test priority ordering with various configurations
- Test gapping service with sample sequences
- Test feature flag toggling

### Integration Tests

- End-to-end custom sequence addition
- Full pipeline with multiple providers
- Offline operation (network disabled)
- Feature flag mode uses G3 when false (no automatic fallback)

### Regression Tests

- Existing Sadie IgBLAST test suite
- Reference system compatibility
- HMM builder with new paths

**Test Data**: Curated dataset in `src/sadie/germlines/tests/data/` with 5-10 genes per segment covering edge cases

---

## Success Metrics Checklist

From spec.md success criteria (SC-001 to SC-015):

- [ ] SC-001: Custom germline in <5 minutes (T065)
- [ ] SC-002: Offline operation works (T066)
- [ ] SC-003: Download scripts <10 minutes (T067)
- [ ] SC-004: 100% test pass rate (T064)
- [ ] SC-005: <500MB disk usage (T068)
- [ ] SC-006: Priority ordering verified (T053, T054)
- [ ] SC-007: Feature flag works (T050)
- [ ] SC-008: Setup in <30 minutes (user testing post-release)
- [ ] SC-009: Change detection works (T069)
- [ ] SC-010: Clear error messages (T070)
- [ ] SC-011: Timing metrics logged (T071)
- [ ] SC-012: CI tests <5 minutes (T072)
- [ ] SC-013-015: Qualitative (code review + user testing)

---

## Risk Mitigation Checklist

- [ ] Biopython gapping accuracy: T009 aligns to IMGT templates; fallback to ungapped on failure
- [ ] VDJbase format changes: T062 adds error handling with clear remediation
- [ ] G3 parity: T044 adapter maintains format, T050 tests G3 mode
- [ ] Test data completeness: T004 curated dataset, documented limitations
- [ ] Performance: T075 profiling, T065 <5min verification

---

## Notes for Implementation

**MVP First**: Implement US1 + US6 first for fastest time-to-value

**Test Coverage**: Unit tests cover ~80% with curated test dataset; full-scale testing with complete reference data done manually

**Feature Flag Timeline**:
1. Initial: SADIE_USE_GERMLINES_MODULE=true (default), G3 available via feature flag (no automatic fallback)
2. Validation period: Monitor usage, collect feedback
3. Deprecation: Announce G3 removal timeline
4. Removal: Delete G3 client code after validation complete

**Documentation**: Update germlines README and INTEGRATION_GUIDE as implementation progresses, not just at end

---

**Generated**: 2026-01-08
**Total Duration Estimate**: 24-28 hours
**MVP Estimate**: 8-10 hours (US1 + US6)
