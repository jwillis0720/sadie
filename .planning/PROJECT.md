# Project: Germline Database Integration

## What This Is

Connect SADIE's germline database module to existing `sadie.airr.Airr` and `sadie.renumbering.Renumbering` modules, enabling selection of germline providers (IMGT, OGRDB, VDJbase, custom) via a `germline_backend` parameter. Supports 29 species with offline operation capability.

## Core Value

Enable researchers to select which germline database their AIRR annotation and antibody renumbering uses, supporting offline operation and custom germline databases instead of being limited to the default G3/IMGT source.

## Requirements

### Validated

- ✓ PROV-01: Specify germline provider for AIRR — v1.0
- ✓ PROV-02: Specify germline provider for Renumbering — v1.0
- ✓ PROV-03: Expose `germline_backend` parameter — v1.0
- ✓ PROV-04: Default priority (custom > ogrdb > vdjbase > imgt) — v1.0
- ✓ PROV-05: Use local germlines module data — v1.0
- ✓ PROV-06: Single provider selection per run — v1.0
- ✓ ERR-01: Clear error when provider lacks species — v1.0
- ✓ ERR-02: Validate custom germline ingestion — v1.0
- ✓ ERR-03: Gapped AA sequence availability — v1.0
- ✓ ERR-04: No silent G3 fallback — v1.0
- ✓ COMPAT-01: G3 remains default backend — v1.0
- ✓ COMPAT-02: Existing AIRR tests pass — v1.0
- ✓ COMPAT-03: Existing renumbering tests pass — v1.0
- ✓ COMPAT-04: Output format identical — v1.0
- ✓ TEST-01: tests/unit/germlines/ directory — v1.0
- ✓ TEST-02: Mirrored renumbering tests — v1.0
- ✓ TEST-03: Species/chains/segments parity — v1.0
- ✓ PERF-01: Performance equivalent to G3 — v1.0

### Active

(None — v1.0 complete)

### Out of Scope

- T-cell receptor (TR) germlines — focus on immunoglobulin (IG) only
- Real-time provider synchronization — manual update via CLI
- Multi-provider blending per analysis — single provider per run
- GUI for provider selection — CLI and API only

## Context

**Shipped:** v1.0 MVP (2026-01-22)
**Codebase:** ~12,440 LOC Python (germlines module + tests)
**Tech Stack:** Python 3.10+, pyhmmer, Biopython, pydantic, pytest, click
**Species Coverage:** 29 species with IgBLAST databases

**Key Integration Points:**
- `src/sadie/airr/igblast/germline.py` — IgBLAST path resolution
- `src/sadie/renumbering/aligners/hmmer.py` — HMM generation
- `src/sadie/reference/reference.py` — Reference system
- `src/sadie/germlines/cli.py` — CLI commands

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Feature flag approach | Enables gradual migration without breaking existing code | ✓ Good — SADIE_USE_GERMLINES_MODULE env var |
| G3 format adapter | Maintains compatibility with Reference system | ✓ Good — GermlineToG3Adapter class |
| Stockholm HMM building | Matches G3 workflow exactly for result parity | ✓ Good — LocalHMMBuilder with pyhmmer |
| Mirrored test suite | Validates integration without modifying existing tests | ✓ Good — tests/unit/germlines/ |
| Single provider per run | Simplifies logic, prevents confusing mixed results | ✓ Good — Constitution principle enforced |
| CLI population command | Programmatic data management | ✓ Good — `sadie germlines populate` |

## Constitution Principles

1. **Provider-Based Architecture**: Providers remain independent; no cross-provider dependencies
2. **Priority-Based Merging (NON-NEGOTIABLE)**: custom > ogrdb > vdjbase > imgt; no per-segment mixing
3. **Local-First Operation**: Runtime uses local data only; offline-capable
4. **Staged Pipeline**: sources → normalized → igblast pipeline respected
5. **Integration Compatibility**: Backward compatibility preserved; API formats consistent

## Constraints

- G3 remains default backend (backwards compatibility)
- No silent fallback to G3 on germlines failure
- Single provider selection per run (no per-segment mixing)
- Default priority: custom > ogrdb > vdjbase > imgt

---
*Last updated: 2026-01-22 after v1.0 milestone*
