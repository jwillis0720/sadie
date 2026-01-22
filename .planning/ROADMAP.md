# Roadmap: Germline Database Integration

## Overview

| Phase | Name | Goal | Status |
|-------|------|------|--------|
| 1 | Setup | Verify germlines module ready and databases built | ✅ Complete |
| 2 | Foundational | Create shared adapters and integration utilities | ✅ Complete |
| 3 | US1: AIRR | Enable germline provider selection for AIRR annotation | ✅ Complete |
| 4 | US2: Renumbering | Enable germline provider selection for renumbering | ✅ Complete |
| 5 | Reference | Enable Reference system to query germlines module | ✅ Complete |
| 6 | US3: Tests | Create mirrored test suite using germlines backend | ✅ Complete |
| 7 | US4: Offline | Verify offline operation capability | ✅ Complete |
| 8 | Polish | Documentation, performance validation, finalization | ✅ Complete |
| 9 | Compliance | Close requirement gaps and enforce constitution | ✅ Complete (8/8) |
| 10 | Species Expansion | Populate IgBLAST databases for all IMGT-supported species | ✅ Complete |
| 11 | IMGT Gapped Fix | Fix IMGT provider to load gapped sequences from *_gapped.fasta files | ✅ Complete |
| 12 | Provider Auto-Population | Make all providers auto-populate via CLI command | ✅ Complete |

**Progress**: 92/92 tasks (100%)

---

## Phase 1: Setup

**Goal**: Verify germlines module installation and database availability

**Requirements**: Foundation for all integration work

**Success Criteria**:
1. Germlines module installed at src/sadie/germlines/
2. Normalized data exists for human species
3. BLAST databases built and verified
4. Gapped sequences available for V/J genes

**Tasks**:
- [x] T001: Verify germlines module installation
- [x] T002: Verify normalized germline data exists
- [x] T003: Build IgBLAST databases
- [x] T004: Verify BLAST database files exist
- [x] T004a: Verify gapped AA/NT sequences for all V/J genes

---

## Phase 2: Foundational

**Goal**: Create shared adapters and utilities used by all user stories

**Requirements**: PROV-03, COMPAT-01

**Success Criteria**:
1. G3 adapter converts GermlineGene to G3 format
2. LocalHMMBuilder generates HMMs from germlines
3. Feature flag controls backend selection

**Tasks**:
- [x] T005: Create G3 adapter class
- [x] T006: Implement to_g3_format() method
- [x] T007: Implement to_g3_format_batch()
- [x] T008: Create LocalHMMBuilder class
- [x] T009: Implement feature flag function

---

## Phase 3: US1 - AIRR Provider Selection

**Goal**: Enable researchers to select germline database for AIRR annotation

**Requirements**: PROV-01, PROV-05

**Success Criteria**:
1. User can run AIRR with provider="imgt" and get IMGT-sourced results
2. User can run AIRR with provider="custom" and get custom gene priority
3. Default provider priority works when unspecified
4. Existing AIRR tests pass

**Tasks**:
- [x] T010: Update GermlineData.__init__() with feature flag
- [x] T011: Add germlines module path logic
- [x] T012: Add legacy G3 fallback logic
- [x] T013: Test IgBLAST paths with feature flag enabled
- [x] T014: Test IgBLAST paths with feature flag disabled
- [x] T015: Run AIRR annotation with germlines backend

---

## Phase 4: US2 - Renumbering Provider Selection

**Goal**: Enable researchers to select germline database for renumbering

**Requirements**: PROV-02, PROV-05

**Success Criteria**:
1. User can run renumbering with germlines backend
2. HMMs built from selected provider's sequences
3. Provider switching reflects in results
4. Existing renumbering tests pass

**Tasks**:
- [x] T016: Implement LocalHMMBuilder.get_hmm()
- [x] T017: Implement LocalHMMBuilder._build_hmm()
- [x] T018: Implement LocalHMMBuilder._get_vj_alignment_pairs()
- [x] T019: Add feature flag check in hmmer.py
- [x] T020: Update HMMER.get_hmm_models()
- [x] T021: Test HMM generation with flag enabled
- [x] T022: Test fallback with flag disabled
- [x] T023: Validate results match between backends

---

## Phase 5: Reference System

**Goal**: Enable Reference system to query germlines with G3 format compatibility

**Requirements**: COMPAT-04

**Success Criteria**:
1. Reference with use_germlines=True returns G3-compatible format
2. Reference with use_germlines=False uses G3 API
3. Output format consistent between backends

**Tasks**:
- [x] T024: Add use_germlines parameter to Reference
- [x] T025: Initialize GermlineManager when enabled
- [x] T026: Update Reference._get_gene()
- [x] T027: Update Reference._get_genes()
- [x] T028: Test Reference with use_germlines=True
- [x] T029: Test Reference with use_germlines=False
- [x] T030: Validate output format consistency

---

## Phase 6: US3 - Test Suite

**Goal**: Create mirrored test suite validating germlines integration

**Requirements**: TEST-01, TEST-02, TEST-04

**Success Criteria**:
1. tests/unit/germlines/ directory exists
2. AIRR integration tests pass
3. Renumbering integration tests pass
4. Reference integration tests pass

**Tasks**:
- [x] T031: Create tests/unit/germlines/ directory
- [x] T032: Create test_airr_integration.py
- [x] T033: Implement test_airr_annotation_with_germlines()
- [x] T034: Implement test_provider_selection()
- [x] T035: Implement test_backwards_compatibility()
- [x] T035a: Implement test_gapped_aa_fallback_translation()
- [x] T036: Create test_renumbering_integration.py
- [x] T037: Implement test_hmmer_with_local_builder()
- [x] T038: Implement test_hmm_caching()
- [x] T039: Create test_reference_integration.py
- [x] T040: Run pytest tests/unit/germlines/

---

## Phase 7: US4 - Offline Operation

**Goal**: Verify AIRR and renumbering work offline

**Requirements**: Offline capability

**Success Criteria**:
1. AIRR annotation succeeds with network disabled
2. Renumbering succeeds with cached HMMs offline
3. Clear errors when germlines not populated
4. Clear errors when BLAST databases not built

**Tasks**:
- [x] T041: Add network isolation fixture
- [x] T042: Test AIRR annotation offline
- [x] T043: Test renumbering offline
- [x] T044: Test error message for unpopulated germlines
- [x] T045: Test error message for missing BLAST databases

---

## Phase 8: Polish

**Goal**: Documentation, performance validation, finalization

**Requirements**: PERF-01, COMPAT-02, COMPAT-03

**Success Criteria**:
1. INTEGRATION_GUIDE.md updated
2. Environment variable documented
3. Performance benchmark completed
4. Full test suite passes

**Tasks**:
- [x] T046: Update INTEGRATION_GUIDE.md
- [x] T047: Document environment variable
- [x] T048: Document germline_backend parameter
- [x] T049: Add performance benchmark
- [x] T050: Run and document performance comparison
- [x] T051: Update quickstart.md
- [x] T052: Run full existing test suite

---

## Phase 9: Compliance & Coverage Gaps

**Goal**: Close remaining requirement gaps and enforce constitution-aligned behaviors

**Requirements**: PROV-04, PROV-06, ERR-01, ERR-02, ERR-03, ERR-04, TEST-03

**Success Criteria**:
1. Single-provider enforcement implemented and tested
2. Missing-species errors clear and informative
3. Custom ingestion validation rejects invalid data
4. No silent G3 fallback when germlines selected
5. Species/chain/segment parity verified
6. Default priority order tested
7. Gapped AA fail-fast implemented

**Tasks**:
- [x] T053: Enforce single-provider selection validation (FR-014)
- [x] T054: Add tests rejecting per-segment provider parameters
- [x] T055: Implement clear error when provider lacks species data (FR-006, NFR-002)
- [x] T056: Validate custom germline ingestion (FR-012)
- [x] T057: Verify species/chain/segment parity (FR-010)
- [x] T058: Test default priority order (FR-004)
- [x] T059: Add negative test for no G3 fallback (NFR-002)
- [x] T060: Add fail-fast when both gapped AA and gapped NT missing (FR-013)

---

## Phase 10: Species Expansion

**Goal**: Populate germlines IgBLAST databases for all IMGT-supported species

**Requirements**: Enable AIRR/renumbering analysis for all species supported by IMGT V-QUEST reference directory

**Success Criteria**:
1. IMGT data downloaded for all 33 mapped species
2. IgBLAST BLAST databases built for each species with available data
3. Auxiliary files (*.aux) created for J gene CDR3 positions
4. organism.yaml updated with all species configurations
5. Multi-species AIRR annotation tests pass
6. Multi-species renumbering tests pass

**Species Coverage** (from download_imgt.py SPECIES_MAP):
- Primates: human, rhesus_macaque, cynomolgus, gorilla, chimpanzee, orangutan_sumatran, orangutan_bornean, lemur, owl_monkey
- Rodents: mouse, mouse_c57bl6j, rat, naked_mole_rat
- Carnivores: dog, cat, ferret, mink
- Ungulates: rabbit, pig, cow, sheep, goat, horse, alpaca, camel
- Birds: chicken
- Fish: zebrafish, atlantic_salmon, rainbow_trout, atlantic_cod, channel_catfish
- Marine mammals: dolphin
- Monotremes: platypus

**Tasks**:
- [x] T061: Download IMGT data for all SPECIES_MAP species using download_imgt.py
- [x] T062: Create auxiliary file generator for J gene CDR3 start positions
- [x] T063: Build IgBLAST databases for all downloaded species
- [x] T064: Generate auxiliary files (*.aux) for each species
- [x] T065: Update organism.yaml with all species configurations
- [x] T066: Verify BLAST database integrity for all species
- [x] T067: Test AIRR annotation with mouse species
- [x] T068: Test AIRR annotation with non-human primate (rhesus_macaque)
- [x] T069: Test AIRR annotation with non-mammalian species (chicken or zebrafish)
- [x] T070: Test renumbering HMM generation for mouse
- [x] T071: Test renumbering HMM generation for rabbit
- [x] T072: Add multi-species integration test suite

---

## Phase 11: IMGT Gapped Fix

**Goal**: Fix IMGT provider to load gapped sequences from `*_gapped.fasta` files, enabling HMM generation for all species with gapped data

**Requirements**: Enable rabbit, chicken, and other species to build HMMs for renumbering

**Background**:
The IMGT provider currently only reads `IGHV.fasta` (ungapped) and ignores `IGHV_gapped.fasta` files.
- Human works because human's main file happens to contain dots (gapped)
- Rabbit/chicken fail because their main files are ungapped
- All 29 species have gapped files that should be loaded

**Success Criteria**:
1. IMGT provider loads `*_gapped.fasta` files (like OGRDB provider does)
2. GermlineGene.sequence_gapped populated for all species with gapped data
3. HMM generation works for rabbit, chicken, and other previously failing species
4. All existing tests continue to pass

**Tasks**:
- [x] T073: Add `_get_gapped_fasta_path()` method to IMGT provider
- [x] T074: Add `_load_gapped_sequences()` method to IMGT provider
- [x] T075: Update `fetch_genes()` to merge gapped sequences
- [x] T076: Test rabbit HMM generation now works
- [x] T077: Test chicken HMM generation now works
- [x] T078: Verify all 29 species have gapped sequences loaded (26/29 have IGHV gapped)
- [x] T079: Run full test suite to ensure no regressions (64/64 passed)

**Depends on**: Phase 10

---

## Phase 12: Provider Auto-Population

**Goal**: Make all germline providers (IMGT, OGRDB, VDJbase) auto-populate via CLI command

**Requirements**: Enable programmatic data population for all providers

**Background**:
Currently, provider coverage is uneven and IMGT requires manual script runs:
- IMGT: 33 species but `download()` raises `NotImplementedError`
- OGRDB: 2 species (human, mouse) - `download()` works
- VDJbase: 2 species (human, rhesus_macaque) - `download()` works

**Decisions** (from 12-CONTEXT.md):
- Trigger: Explicit CLI command (`sadie germlines populate`)
- Scope: All species available from each provider's API
- Errors: Fail fast with checkpoint for resume
- Updates: Version-check and update if newer

**Success Criteria**:
1. `sadie germlines populate` CLI command works
2. IMGT provider `download()` method implemented (not NotImplementedError)
3. All providers download all available species
4. Version checking prevents redundant downloads
5. Fail-fast with checkpoint enables resume after errors

**Tasks**:
- [x] T080: Implement `sadie germlines populate` CLI command
- [x] T081: Implement `IMGTProvider.download()` from existing script logic
- [x] T082: Add version tracking for IMGT releases
- [x] T083: Audit and download all OGRDB available species
- [x] T084: Audit and download all VDJbase available species
- [x] T085: Add `--force` flag for re-download
- [x] T086: Add checkpoint/resume for fail-fast recovery
- [x] T087: Add rich progress bars for download tracking
- [x] T088: Integrate post-download build pipeline (BLAST DBs, aux, internal_data)
- [x] T089: Test CLI command with all providers
- [x] T090: Verify downloaded data integrity

**Depends on**: Phase 11

---

## Dependencies

```
Phase 1 → Phase 2
              ↓
    ┌─────────┼─────────┐
    ↓         ↓         ↓
Phase 3   Phase 4   Phase 5
    ↓         ↓         ↓
    └─────────┼─────────┘
              ↓
          Phase 6
              ↓
          Phase 7
              ↓
          Phase 8
              ↓
          Phase 9
              ↓
          Phase 10
              ↓
          Phase 11
              ↓
          Phase 12
```

---
*Last updated: 2026-01-22*
