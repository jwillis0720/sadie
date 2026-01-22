# Feature Specification: Germline Database Integration

**Feature Branch**: `002-germline-integration`  
**Created**: 2026-01-19  
**Status**: Planning Complete  
**Input**: Connect SADIE's new germline database to existing sadie.airr.Airr and sadie.numbering.Renumbering to allow picking new germlines database, plus create tests in tests/unit/germlines mirroring airr and renumbering tests using backend germline instead of default IMGT from G3

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Select Germline Provider for AIRR Analysis (Priority: P1)

As a researcher running immune repertoire analysis, I want to select which germline database (IMGT, OGRDB, VDJbase, or custom) my AIRR annotation uses, so that I can leverage alternative gene databases or my own custom germlines instead of being limited to the default IMGT from G3.

**Why this priority**: This is the core integration point - enabling germline provider selection in the primary AIRR analysis workflow unlocks all downstream functionality.

**Independent Test**: Can be tested by running AIRR annotation with different provider parameters and verifying gene calls come from the specified database source.

**Acceptance Scenarios**:

1. **Given** a user has sequence data and the new germlines module is populated, **When** they run AIRR analysis with `provider="imgt"`, **Then** gene annotations reference IMGT-sourced germlines from the local database.

2. **Given** a user has custom germline sequences in `sources/custom/`, **When** they run AIRR analysis with `provider="custom"`, **Then** gene annotations include their custom genes with priority over default providers.

3. **Given** a user does not specify a provider, **When** they run AIRR analysis, **Then** the system uses the default priority order (custom > ogrdb > vdjbase > imgt) from the germlines module.

---

### User Story 2 - Select Germline Provider for Renumbering (Priority: P1)

As a researcher performing antibody numbering, I want to select which germline database my renumbering/HMM alignment uses, so that my numbering schemes align against the same germlines used in my analysis pipeline.

**Why this priority**: Renumbering is a parallel critical path alongside AIRR - both need germline integration for a complete solution.

**Independent Test**: Can be tested by running renumbering with different provider parameters and verifying HMM alignments use the specified germline source.

**Acceptance Scenarios**:

1. **Given** a user has antibody sequences, **When** they run renumbering with the germlines backend, **Then** HMM models are built from the selected germline provider's sequences.

2. **Given** a user switches germline providers between analyses, **When** they run renumbering, **Then** results reflect the currently selected provider's gene database.

---

### User Story 3 - Consistent Test Suite Using Germlines Backend (Priority: P2)

As a developer maintaining SADIE, I want a test suite in `tests/unit/germlines/` that mirrors the existing airr and renumbering tests but uses the local germlines backend, so that I can validate the germlines integration produces equivalent results to the G3-based tests.

**Why this priority**: Test coverage ensures the integration is correct and provides regression protection for future changes.

**Independent Test**: Can be tested by running the new test suite and comparing results with existing G3-based tests.

**Acceptance Scenarios**:

1. **Given** the existing `tests/unit/airr/` test suite, **When** I run the mirrored tests in `tests/unit/germlines/`, **Then** all tests pass using the local germlines backend instead of G3.

2. **Given** the existing `tests/unit/renumbering/` test suite, **When** I run the mirrored renumbering tests with germlines backend, **Then** all tests pass with equivalent results.

---

### User Story 4 - Offline Germline Operation (Priority: P3)

As a researcher working in an environment without reliable internet, I want AIRR and renumbering to work completely offline using the local germlines database, so that I am not dependent on G3 API availability.

**Why this priority**: Offline capability is a key benefit of the germlines module but depends on P1/P2 being complete first.

**Independent Test**: Can be tested by disabling network access and verifying all analysis workflows complete successfully.

**Acceptance Scenarios**:

1. **Given** the germlines module is populated and network is unavailable, **When** I run AIRR analysis, **Then** annotation completes successfully using local data.

2. **Given** cached HMM models exist, **When** I run renumbering offline, **Then** numbering completes without network errors.

---

### Edge Cases

- Requested provider has no data for specified species: Fail with clear error message (FR-005); no silent fallback to G3.
- Conflicting gene names across providers: Highest-priority provider wins; lower-priority duplicates are dropped with a WARNING log (per Principle II).
- Provider switching mid-analysis: Allowed; user responsible for maintaining consistency across linked datasets.
- Custom germlines with invalid sequences (malformed FASTA, invalid amino acids): Validate at ingestion; reject with detailed error message identifying the specific validation failure.
- First-time setup with unpopulated germlines: Fail with clear error message and setup instructions (no silent fallback).
- Germlines backend query failure: No automatic fallback to G3; fail with clear error since user explicitly selected germlines backend.
- Gene missing both gapped AA and gapped NT sequences: Fail with clear error message listing the gene name and suggesting data population via `update_databases()`.
- Provider missing a required segment (e.g., D only in another provider): Fail with clear error; do not mix providers per segment.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST allow users to specify a germline provider (e.g., "imgt", "ogrdb", "vdjbase", "custom") when initializing AIRR analysis.

- **FR-002**: System MUST allow users to specify a germline provider when initializing Renumbering operations.

- **FR-003**: System MUST expose a `germline_backend` parameter accepting values "g3" (default) or "germlines" to select the backend system for germline data.

- **FR-004**: System MUST use a default provider selection list (custom > ogrdb > vdjbase > imgt) when no provider is explicitly specified; duplicates are resolved by priority with a WARNING log.

- **FR-005**: System MUST use local germlines module data instead of G3 API when germlines backend is selected.

- **FR-006**: System MUST provide clear error messages when a requested provider has no data for the specified species.

- **FR-007**: System MUST support backwards compatibility - G3 remains the default backend; germlines backend is opt-in via an explicit parameter. Existing code using default settings continues to work without changes.

- **FR-008**: System MUST include a new test directory `tests/unit/germlines/` containing critical path tests that mirror core AIRR annotation functionality using germlines backend.

- **FR-009**: System MUST include critical path tests that mirror core renumbering functionality using germlines backend.

- **FR-010**: System MUST support the same species, chains, and segments currently supported by the existing AIRR and Renumbering modules.

- **FR-011**: System MUST maintain identical output format/schema regardless of which germline provider is used.

- **FR-012**: System MUST validate custom germline sequences at ingestion time, rejecting invalid sequences (malformed FASTA, invalid amino acids) with detailed error messages identifying the specific validation failure.

- **FR-013**: System MUST ensure gapped amino acid sequences (`sequence_aa_gapped`) are available for all V and J genes used in HMM building. If pre-computed gapped AA is unavailable, the system MUST translate from gapped nucleotide sequences (`sequence_gapped`) at runtime. If neither is available, the system MUST fail with a clear error message identifying the specific gene lacking sequence data.

- **FR-014**: System MUST apply a single provider selection configuration to all V/D/J segments within a run; users cannot specify different providers per segment. The selection may be a single provider or an ordered provider list.

### Non-Functional Requirements

- **NFR-001**: Germline lookup performance using the local germlines backend MUST be equivalent to G3 backend (no regression from current behavior).

- **NFR-002**: When the germlines backend is explicitly selected, the system MUST NOT silently fall back to G3 on failure; it MUST fail with a clear error message.

### Key Entities

- **GermlineProvider**: Represents a source of germline data (IMGT, OGRDB, VDJbase, custom) with its priority level and available species.

- **GermlineSelection**: Configuration specifying which provider(s) to use for an analysis run, including priority ordering.

- **GermlineData**: The germline gene information used by IgBLAST including V, D, J segment databases and auxiliary files - now sourced from local germlines module.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Critical path AIRR annotation tests are mirrored in `tests/unit/germlines/` with passing results using germlines backend.

- **SC-002**: Critical path renumbering tests are mirrored with passing results using germlines backend.

- **SC-003**: Users can complete AIRR analysis with any available germline provider without errors.

- **SC-004**: Users can complete renumbering with any available germline provider without errors.

- **SC-005**: Gene annotation results using IMGT provider from germlines module match results from G3-based IMGT for the same input sequences.

- **SC-006**: System operates fully offline after initial germlines database population.

- **SC-007**: No breaking changes to existing user workflows - all existing tests in `tests/unit/airr/` and `tests/unit/renumbering/` continue to pass.

## Clarifications

### Session 2026-01-21

- Q: How should the system handle custom germlines with invalid sequences? → A: Validate at ingestion; reject invalid sequences with detailed error
- Q: Should there be a performance expectation for germline lookups compared to G3? → A: Equivalent to G3 (no regression from current behavior)
- Q: If the germlines backend fails for a specific query, should the system fall back to G3? → A: No fallback; fail with clear error (user chose germlines explicitly)
- Q: What parameter name should be used to specify the germline backend? → A: `germline_backend` (values: "g3" or "germlines")
- Q: Can V/D/J segments use different providers within one run? → A: No; provider selection is made once per run and applies uniformly across segments (single provider or ordered provider list).

**Planning Phase Complete** - See [plan.md](./plan.md) for implementation artifacts.

### Session 2026-01-19

- Q: How does the system handle conflicting gene names across multiple providers? → A: Highest-priority provider wins; lower-priority duplicates dropped with WARNING log.
- Q: How does backwards compatibility work for existing code? → A: G3 remains default; germlines is opt-in via parameter
- Q: What scope of tests should be mirrored? → A: Mirror critical path tests only (AIRR annotation, renumbering core)
- Q: What happens when germlines module is not populated (first-time setup)? → A: Fail with clear error message and setup instructions
- Q: What happens when switching providers mid-analysis for linked datasets? → A: Allow switching; user responsibility for consistency

## Constitution Alignment

- **Principle I — Provider-Based Architecture**: Providers remain independent and self-contained; integration uses GermlineProvider interfaces without cross-provider dependencies.
- **Principle II — Priority-Based Merging (NON-NEGOTIABLE)**: Default priority (custom > ogrdb > vdjbase > imgt) enforced; no per-segment mixing; highest-priority wins for duplicates with WARNING log.
- **Principle III — Local-First Operation**: All runtime operations use local data only; offline scenarios covered in User Story 4 with no network reliance.
- **Principle IV — Staged Pipeline Architecture**: Processing follows sources → normalized → igblast; integration does not bypass pipeline stages.
- **Principle V — Integration Compatibility**: Backward compatibility preserved (G3 remains default, feature-flag opt-in); API formats remain consistent across backends.
- **Status**: No known constitution violations.

## Assumptions

- The germlines module (`src/sadie/germlines/`) is already implemented and populated with IMGT data as per the `001-germline-completion` feature.
- The INTEGRATION_GUIDE.md in the germlines module provides accurate integration patterns to follow.
- Existing test fixtures are compatible with the germlines backend or can be adapted.
- The priority ordering (custom > ogrdb > vdjbase > imgt) is the correct default based on project requirements.
