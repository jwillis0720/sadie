# Project: Germline Database Integration

## What This Is

Connect SADIE's germline database module to existing `sadie.airr.Airr` and `sadie.renumbering.Renumbering` modules, enabling selection of germline providers (IMGT, OGRDB, VDJbase, custom) via a `germline_backend` parameter. The integration uses a feature flag for backwards compatibility, with G3 API remaining default while local germlines is opt-in.

## Core Value

Enable researchers to select which germline database their AIRR annotation and antibody renumbering uses, supporting offline operation and custom germline databases instead of being limited to the default G3/IMGT source.

## Context

**Repository**: sadie (Python bioinformatics library)
**Branch**: `002-germline-integration`
**Prerequisite**: 001-germline-completion (germlines module exists and is populated)
**Tech Stack**: Python 3.10+, pyhmmer, Biopython, pydantic, pytest

## User Stories

### US1 (P1): Select Germline Provider for AIRR Analysis
As a researcher running immune repertoire analysis, I want to select which germline database (IMGT, OGRDB, VDJbase, or custom) my AIRR annotation uses.

### US2 (P1): Select Germline Provider for Renumbering
As a researcher performing antibody numbering, I want to select which germline database my renumbering/HMM alignment uses.

### US3 (P2): Consistent Test Suite Using Germlines Backend
As a developer, I want a test suite in `tests/unit/germlines/` that mirrors critical path tests using the local germlines backend.

### US4 (P3): Offline Germline Operation
As a researcher in an environment without reliable internet, I want AIRR and renumbering to work completely offline.

## Constraints

- G3 remains default backend (backwards compatibility)
- No silent fallback to G3 on germlines failure (NFR-002)
- Single provider selection per run (no per-segment mixing)
- Default priority: custom > ogrdb > vdjbase > imgt

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Feature flag approach | Enables gradual migration without breaking existing code | SADIE_USE_GERMLINES_MODULE env var |
| G3 format adapter | Maintains compatibility with Reference system | GermlineToG3Adapter class |
| Stockholm HMM building | Matches G3 workflow exactly for result parity | LocalHMMBuilder with pyhmmer |
| Mirrored test suite | Validates integration without modifying existing tests | tests/unit/germlines/ |

## Constitution Principles

1. **Provider-Based Architecture**: Providers remain independent; no cross-provider dependencies
2. **Priority-Based Merging (NON-NEGOTIABLE)**: custom > ogrdb > vdjbase > imgt; no per-segment mixing
3. **Local-First Operation**: Runtime uses local data only; offline-capable
4. **Staged Pipeline**: sources → normalized → igblast pipeline respected
5. **Integration Compatibility**: Backward compatibility preserved; API formats consistent

## Source Specification

Full specification from spec-kit: `specs/002-germline-integration/`
- spec.md - Feature specification
- plan.md - Implementation plan
- tasks.md - Task breakdown
- research.md - Research findings
- data-model.md - Data model documentation

---
*Last updated: 2026-01-21 after conversion from spec-kit*
