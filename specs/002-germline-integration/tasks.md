# Tasks: Germline Database Integration

**Input**: Design documents from `/specs/002-germline-integration/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, contracts/

**Tests**: Tests are included as specified in FR-008, FR-009 (mirrored test suite requirements).

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story?] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (US1, US2, US3, US4)
- Include exact file paths in descriptions

---

## Summary

**Feature**: 002-germline-integration
**Branch**: `002-germline-integration`
**Tech Stack**: Python 3.10+, pyhmmer, Biopython, pydantic, pytest
**Integration Points**: 3 (IgBLAST paths, HMM generation, Reference system)

**User Stories** (from spec.md):
- **US1 (P1)**: Select Germline Provider for AIRR Analysis
- **US2 (P1)**: Select Germline Provider for Renumbering
- **US3 (P2)**: Consistent Test Suite Using Germlines Backend
- **US4 (P3)**: Offline Germline Operation

**Clarifications Applied** (from spec.md Session 2026-01-21):
- Custom germlines with invalid sequences: Validate at ingestion; reject with detailed error (FR-012)
- Performance expectation: Equivalent to G3 (NFR-001)
- Fallback behavior: No fallback to G3; fail with clear error (NFR-002)
- Backend parameter: `germline_backend` with values "g3" (default) or "germlines" (FR-003)
- Single-provider rule: One provider per run for all V/D/J segments; no per-segment mixing (FR-014)

## Constitution Alignment

- Principle I (Provider-Based Architecture): Providers remain independent; no cross-provider dependencies.
- Principle II (Priority-Based Merging, NON-NEGOTIABLE): Default priority custom > ogrdb > vdjbase > imgt; no per-segment mixing; duplicates drop lower priority with warning.
- Principle III (Local-First Operation): Runtime uses local data only; offline flows validated in US4.
- Principle IV (Staged Pipeline): sources → normalized → igblast respected; no bypass.
- Principle V (Integration Compatibility): G3 remains default via feature flag; public formats unchanged across backends.

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Verify germlines module is ready and databases are built

**⚠️ External Dependency**: BLAST+ tools (makeblastdb) required for database building

- [x] T001 Verify germlines module installation at src/sadie/germlines/__init__.py
- [x] T002 Verify normalized germline data exists in src/sadie/germlines/normalized/human/
- [x] T003 Build IgBLAST databases using update_databases("human") in src/sadie/germlines/pipeline.py
- [x] T004 Verify BLAST database files exist at src/sadie/germlines/igblast/database/human/
- [ ] T004a [US1/US2] Verify gapped AA or gapped NT sequences exist for all V/J genes via: `python -c "from sadie.germlines import get_manager; m=get_manager(); genes=[g for g in m.get_genes('human','V','H') if not g.sequence_aa_gapped and not g.sequence_gapped]; print(f'Missing gapped data: {[g.name for g in genes]}' if genes else 'All V/J genes have gapped sequences')"`

**Checkpoint**: Setup verified - databases exist and are built

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Create shared adapters and integration utilities used by ALL user stories

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [x] T005 [P] Create G3 adapter class in src/sadie/germlines/g3_adapter.py
- [x] T006 [P] Implement GermlineToG3Adapter.to_g3_format() method with region extraction in src/sadie/germlines/g3_adapter.py
- [x] T007 [P] Implement GermlineToG3Adapter.to_g3_format_batch() for bulk conversion in src/sadie/germlines/g3_adapter.py
- [x] T008 [P] Create LocalHMMBuilder class in src/sadie/germlines/renumbering_integration.py
- [x] T009 [P] Implement feature flag function use_germlines_module() in src/sadie/germlines/utils/feature_flags.py

**Checkpoint**: Foundation ready - user story implementation can now begin in parallel

---

## Phase 3: User Story 1 - Select Germline Provider for AIRR Analysis (Priority: P1) 🎯 MVP

**Goal**: Enable researchers to select germline database providers (IMGT, OGRDB, VDJbase, custom) for AIRR annotation

**Independent Test**: Run AIRR annotation with different provider parameters and verify gene calls come from the specified database source

**Acceptance Criteria**:
- User runs AIRR with provider="imgt" → annotations reference IMGT germlines
- User runs AIRR with provider="custom" → custom genes included with priority
- User runs AIRR without provider → default priority order (custom > ogrdb > vdjbase > imgt)
- Existing tests in tests/unit/airr/ continue to pass

### Implementation for User Story 1

- [x] T010 [P] [US1] Update GermlineData.__init__() with feature flag check in src/sadie/airr/igblast/germline.py
- [x] T011 [P] [US1] Add germlines module path logic for v_gene_dir, d_gene_dir, j_gene_dir in src/sadie/airr/igblast/germline.py
- [x] T012 [P] [US1] Add legacy G3 path fallback logic when feature flag disabled in src/sadie/airr/igblast/germline.py
- [x] T013 [US1] Test IgBLAST paths switch correctly when SADIE_USE_GERMLINES_MODULE=true
- [x] T014 [US1] Test IgBLAST paths use legacy G3 when SADIE_USE_GERMLINES_MODULE=false (backwards compatibility)
- [x] T015 [US1] Run AIRR annotation with germlines backend and validate v_call results

**Checkpoint**: User Story 1 complete - AIRR annotation works with germlines backend

---

## Phase 4: User Story 2 - Select Germline Provider for Renumbering (Priority: P1)

**Goal**: Enable researchers to select germline database providers for renumbering/HMM alignment

**Independent Test**: Run renumbering with different provider parameters and verify HMM alignments use the specified germline source

**Acceptance Criteria**:
- User runs renumbering with germlines backend → HMMs built from germline sequences
- User switches providers → results reflect new provider
- Existing tests in tests/unit/renumbering/ continue to pass

### Implementation for User Story 2

- [x] T016 [P] [US2] Implement LocalHMMBuilder.get_hmm() with cache check in src/sadie/germlines/renumbering_integration.py
- [x] T017 [P] [US2] Implement LocalHMMBuilder._build_hmm() with Stockholm generation in src/sadie/germlines/renumbering_integration.py
- [x] T018 [P] [US2] Implement LocalHMMBuilder._get_vj_alignment_pairs() querying GermlineManager in src/sadie/germlines/renumbering_integration.py
- [x] T019 [P] [US2] Add _use_local_hmm_builder() feature flag check in src/sadie/renumbering/aligners/hmmer.py
- [x] T020 [US2] Update HMMER.get_hmm_models() with LocalHMMBuilder integration in src/sadie/renumbering/aligners/hmmer.py
- [x] T021 [US2] Test renumbering HMM generation with feature flag enabled
- [x] T022 [US2] Test renumbering fallback with feature flag disabled (backwards compatibility)
- [x] T023 [US2] Validate renumbering results match between G3 and germlines backends

**Checkpoint**: User Story 2 complete - Renumbering works with germlines backend

---

## Phase 5: Reference System Integration (Supports US1/US2)

**Goal**: Enable Reference system to query germlines module with G3 format compatibility

**Note**: This supports both US1 and US2 by providing consistent germline data access

- [x] T024 [P] Add use_germlines parameter to Reference.__init__() in src/sadie/reference/reference.py
- [x] T025 [P] Initialize GermlineManager and G3Adapter when use_germlines=True in src/sadie/reference/reference.py
- [x] T026 Update Reference._get_gene() to query germlines module when enabled in src/sadie/reference/reference.py
- [x] T027 Update Reference._get_genes() to query germlines module when enabled in src/sadie/reference/reference.py
- [x] T028 Test Reference with use_germlines=True returns G3-compatible format
- [x] T029 Test Reference with use_germlines=False uses G3 API (default behavior)
- [x] T030 Validate Reference output format consistency between backends

**Checkpoint**: Reference system supports both backends with format consistency

---

## Phase 6: User Story 3 - Consistent Test Suite Using Germlines Backend (Priority: P2)

**Goal**: Create test suite in tests/unit/germlines/ that mirrors critical path tests from airr and renumbering

**Independent Test**: Run the new test suite and compare results with existing G3-based tests

**Acceptance Criteria**:
- All mirrored tests in tests/unit/germlines/ pass with germlines backend
- Test results equivalent to existing G3-based tests
- Tests run in CI alongside existing test suites

### Tests for User Story 3

- [x] T031 [P] [US3] Create tests/unit/germlines/ directory structure
- [x] T032 [P] [US3] Create test_airr_integration.py in tests/unit/germlines/test_airr_integration.py
- [x] T033 [P] [US3] Implement test_airr_annotation_with_germlines() in tests/unit/germlines/test_airr_integration.py
- [x] T034 [P] [US3] Implement test_provider_selection() in tests/unit/germlines/test_airr_integration.py
- [x] T035 [P] [US3] Implement test_backwards_compatibility() in tests/unit/germlines/test_airr_integration.py
- [ ] T035a [P] [US3] Implement test_gapped_aa_fallback_translation() in tests/unit/germlines/test_renumbering_integration.py - verify HMM builds correctly when only gapped NT available (no gapped AA)
- [x] T036 [P] [US3] Create test_renumbering_integration.py in tests/unit/germlines/test_renumbering_integration.py
- [x] T037 [P] [US3] Implement test_hmmer_with_local_builder() in tests/unit/germlines/test_renumbering_integration.py
- [x] T038 [P] [US3] Implement test_hmm_caching() in tests/unit/germlines/test_renumbering_integration.py
- [x] T039 [P] [US3] Create test_reference_integration.py in tests/unit/germlines/test_reference_integration.py
- [x] T040 [US3] Run pytest tests/unit/germlines/ and ensure all tests pass

**Checkpoint**: User Story 3 complete - Mirrored test suite validates integration

---

## Phase 7: User Story 4 - Offline Germline Operation (Priority: P3)

**Goal**: Ensure AIRR and renumbering work completely offline using local germlines database

**Independent Test**: Disable network access and verify all analysis workflows complete successfully

**Acceptance Criteria**:
- With germlines populated and network disabled → AIRR annotation succeeds
- With cached HMMs and network disabled → renumbering succeeds
- Clear error messages if germlines not populated (first-time setup)

### Implementation for User Story 4

- [x] T041 [P] [US4] Add network isolation fixture in tests/unit/germlines/conftest.py
- [x] T042 [US4] Test AIRR annotation with network disabled (must succeed with germlines)
- [x] T043 [US4] Test renumbering with network disabled (must succeed with cached HMMs)
- [x] T044 [US4] Test clear error message when germlines not populated
- [x] T045 [US4] Test clear error message when BLAST databases not built

**Checkpoint**: User Story 4 complete - Offline operation verified

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, performance validation, and finalization

- [x] T046 [P] Update INTEGRATION_GUIDE.md with Reference system details in src/sadie/germlines/INTEGRATION_GUIDE.md
- [x] T047 [P] Document SADIE_USE_GERMLINES_MODULE environment variable in README.md
- [x] T048 [P] Document germline_backend parameter in API docs
- [x] T049 [P] Add performance benchmark comparing germlines vs G3 backends
- [x] T050 Run performance comparison and document results in specs/002-germline-integration/
- [x] T051 Update quickstart.md with any implementation changes in specs/002-germline-integration/quickstart.md
- [x] T052 Run full existing test suite (tests/unit/airr/, tests/unit/renumbering/) to verify no regressions

---

## Phase 9: Compliance & Coverage Gaps

**Purpose**: Close remaining requirement/test gaps and enforce constitution-aligned behaviors

- [x] T053 [P] [US1/US2] Enforce single-provider selection validation (reject per-segment provider parameters) in germlines backend and callers (FR-014)
- [x] T054 [P] [US1/US2] Add AIRR and renumbering tests that reject per-segment provider parameters (FR-014)
- [x] T055 [P] [US1/US2] Implement and test clear error when requested provider lacks species data (FR-006, NFR-002)
- [x] T056 [P] [US1/US2] Validate custom germline ingestion rejects malformed FASTA/invalid amino acids with detailed errors; add tests (FR-012)
- [x] T057 [P] [US1/US2] Verify germlines backend supports same species/chains/segments as existing modules; add assertion tests (FR-010)
- [x] T058 [P] [US1] Test default priority order custom > ogrdb > vdjbase > imgt (FR-004)
- [x] T059 [P] [Cross] Add negative test ensuring germlines backend errors do not fall back to G3 when selected (NFR-002)
- [x] T060 [P] [US3] Add fail-fast logic and tests when both gapped AA and gapped NT sequences are missing for a gene (FR-013)

---

## Phase 10: Species Expansion

**Purpose**: Populate germlines IgBLAST databases for all IMGT-supported species to enable multi-species analysis

**⚠️ External Dependency**: BLAST+ tools (makeblastdb) required for database building

**Species Coverage** (33 species from download_imgt.py SPECIES_MAP):
- Primates: human, rhesus_macaque, cynomolgus, gorilla, chimpanzee, orangutan_sumatran, orangutan_bornean, lemur, owl_monkey
- Rodents: mouse, mouse_c57bl6j, rat, naked_mole_rat
- Carnivores: dog, cat, ferret, mink
- Ungulates: rabbit, pig, cow, sheep, goat, horse, alpaca, camel
- Birds: chicken
- Fish: zebrafish, atlantic_salmon, rainbow_trout, atlantic_cod, channel_catfish
- Marine mammals: dolphin
- Monotremes: platypus

### Data Acquisition

- [ ] T061 [P] Download IMGT data for all SPECIES_MAP species using: `python -m sadie.germlines.scripts.download_imgt --species human mouse rat rabbit rhesus_macaque cynomolgus dog cat pig cow sheep goat horse alpaca camel chicken zebrafish gorilla chimpanzee orangutan_sumatran orangutan_bornean ferret mink dolphin platypus atlantic_salmon rainbow_trout atlantic_cod channel_catfish lemur owl_monkey naked_mole_rat mouse_c57bl6j`

### Database Building

- [ ] T062 [P] Create auxiliary file generator script at src/sadie/germlines/scripts/build_aux_files.py for J gene CDR3 start positions
- [ ] T063 Build IgBLAST BLAST databases for all downloaded species using BlastDBBuilder in src/sadie/germlines/builders/blast.py
- [ ] T064 Generate auxiliary files (*.aux) for each species in src/sadie/germlines/igblast/aux_db/
- [ ] T065 Update organism.yaml with all species configurations in src/sadie/germlines/igblast/internal_data/organism.yaml

### Validation

- [ ] T066 Verify BLAST database integrity for all species (check .nhr, .nin, .nsq files exist)
- [ ] T067 Test AIRR annotation with mouse species in tests/unit/germlines/test_multi_species.py
- [ ] T068 Test AIRR annotation with non-human primate (rhesus_macaque) in tests/unit/germlines/test_multi_species.py
- [ ] T069 Test AIRR annotation with non-mammalian species (chicken or zebrafish) in tests/unit/germlines/test_multi_species.py

### Renumbering Integration

- [ ] T070 Test renumbering HMM generation for mouse in tests/unit/germlines/test_multi_species.py
- [ ] T071 Test renumbering HMM generation for rabbit in tests/unit/germlines/test_multi_species.py
- [ ] T072 Add multi-species integration test suite in tests/unit/germlines/test_multi_species.py

**Checkpoint**: Phase 10 complete - Multi-species support fully operational

---

## Dependencies & Execution Order

### Phase Dependencies

```
Phase 1 (Setup) → Phase 2 (Foundational)
                    ↓
                    ├→ Phase 3 (US1: AIRR) ─────────────┐
                    │                                    │
                    ├→ Phase 4 (US2: Renumbering) ──────┼→ Phase 5 (Reference)
                    │                                    │
                    └────────────────────────────────────┼→ Phase 6 (US3: Tests)
                                                         │
                                                         └→ Phase 7 (US4: Offline)
                                                              ↓
                                                         Phase 8 (Polish)
                                                               ↓
                                                         Phase 9 (Compliance & Coverage)
```

### User Story Dependencies

- **US1 (P1)**: Can start after Phase 2 (Foundational) - No dependencies on other stories
- **US2 (P1)**: Can start after Phase 2 (Foundational) - Parallel with US1
- **US3 (P2)**: Depends on US1 and US2 completion (mirrors their functionality)
- **US4 (P3)**: Depends on US3 (uses test infrastructure)

### Within Each User Story

- Implementation tasks before validation tasks
- Models/utilities before services
- Services before integration
- Core implementation before error handling
- Story complete before moving to next priority

### Parallel Opportunities

**Phase 2 (Foundational)**: T005, T006, T007, T008, T009 can all run in parallel

**Phase 3 (US1)**: T010, T011, T012 can run in parallel; then T013 → T014 → T015

**Phase 4 (US2)**: T016, T017, T018, T019 can run in parallel; then T020 → T021 → T022 → T023

**Phase 5 (Reference)**: T024, T025 can run in parallel; then T026 → T027 → T028 → T029 → T030

**Phase 6 (US3)**: T031-T039 can ALL run in parallel (different files); then T040

**Phase 8 (Polish)**: T046, T047, T048, T049 can run in parallel; then T050 → T051 → T052

---

## Parallel Execution Examples

### Phase 2 - All Parallel
```bash
# All adapter/builder tasks can run simultaneously
parallel_tasks=(T005 T006 T007 T008 T009)
```

### US1 & US2 - Parallel Development
```bash
# After foundational phase, US1 and US2 can run in parallel
# Developer A: US1 tasks T010-T015
# Developer B: US2 tasks T016-T023
```

### Phase 6 - All Test Files Parallel
```bash
# All test file creation can happen in parallel
parallel_tasks=(T031 T032 T033 T034 T035 T036 T037 T038 T039)
# Then sequential: T040 (run all tests)
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup
2. Complete Phase 2: Foundational (CRITICAL - blocks all stories)
3. Complete Phase 3: User Story 1
4. **STOP and VALIDATE**: Test User Story 1 independently
5. Deploy/demo if ready

### Incremental Delivery

1. Complete Setup + Foundational → Foundation ready
2. Add User Story 1 → Test independently → Deploy/Demo (MVP!)
3. Add User Story 2 → Test independently → Deploy/Demo
4. Add User Story 3 → Test independently → Deploy/Demo
5. Add User Story 4 → Test independently → Final release
6. Polish phase → Documentation complete

### Parallel Team Strategy

With multiple developers:

1. Team completes Setup + Foundational together
2. Once Foundational is done:
   - Developer A: User Story 1 (AIRR)
   - Developer B: User Story 2 (Renumbering)
3. Stories complete and integrate independently
4. Single developer: User Story 3 (Tests) then User Story 4 (Offline)

---

## Task Progress

**Total Tasks**: 74
**Completed**: 60 (81%)
**Remaining**: 14 (19%)

| Phase | Total | Complete | Status |
|-------|-------|----------|--------|
| Phase 1: Setup | 5 | 4 | ⏳ 1 remaining (T004a) |
| Phase 2: Foundational | 5 | 5 | ✅ Complete |
| Phase 3: US1 (AIRR) | 6 | 6 | ✅ Complete |
| Phase 4: US2 (Renumbering) | 8 | 8 | ✅ Complete |
| Phase 5: Reference | 7 | 7 | ✅ Complete |
| Phase 6: US3 (Tests) | 11 | 10 | ⏳ 1 remaining (T035a) |
| Phase 7: US4 (Offline) | 5 | 5 | ✅ Complete |
| Phase 8: Polish | 7 | 7 | ✅ Complete |
| Phase 9: Compliance & Coverage Gaps | 8 | 8 | ✅ Complete |
| Phase 10: Species Expansion | 12 | 0 | 🚧 Not started |

---

## Success Criteria Mapping

| Criterion | Tasks | Status |
|-----------|-------|--------|
| SC-001: Mirrored AIRR tests | T032-T035 | ✅ Complete |
| SC-002: Mirrored renumbering tests | T036-T038 | ✅ Complete |
| SC-003: AIRR works with any provider | T013-T015 | ✅ Complete |
| SC-004: Renumbering works with any provider | T021-T023 | ✅ Complete |
| SC-005: Results match G3 IMGT | T015, T023, T030 | ✅ Complete |
| SC-006: Offline operation | T041-T045 | ✅ Complete |
| SC-007: No breaking changes | T014, T022, T029, T052 | ✅ Complete |
| FR-004: Default priority order | T058 | 🚧 Not started |
| FR-006: Missing-species error | T055 | 🚧 Not started |
| FR-010: Species/chain/segment parity | T057 | 🚧 Not started |
| FR-012: Custom ingestion validation | T056 | 🚧 Not started |
| FR-013: Gapped AA availability | T004a, T035a, T060 | ⏳ Partial |
| FR-014: Single provider per run | T053, T054 | 🚧 Not started |
| NFR-002: No silent G3 fallback on germlines failure | T055, T059 | 🚧 Not started |

---

## Notes

**Feature Flag**: `SADIE_USE_GERMLINES_MODULE` environment variable
- Default: unset/false -> uses G3 backend (backwards compatible)
- Set to "true" to opt-in to germlines backend
- If `germline_backend` is provided, it overrides the env flag

**Backend Parameter**: `germline_backend` (per FR-003)
- Values: "g3" (default) or "germlines"
- Available on AIRR and Renumbering classes

**Performance Targets** (NFR-001):
- Gene lookup: <200ms (10x faster than G3)
- HMM build (uncached): <15s (equivalent to G3)
- HMM load (cached): <100ms (100x faster than G3)

**Error Handling** (NFR-002):
- No silent fallback to G3 when germlines backend selected
- Clear error messages with setup instructions
- Validation at ingestion for custom germlines (FR-012)
