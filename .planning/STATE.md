# State: Germline Database Integration

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-22)

**Core value:** Enable researchers to select germline database for AIRR annotation and renumbering
**Current focus:** v1.3 Test Infrastructure & Species Expansion

## Current Position

Phase: 25 — Macaque Germlines Integration
Plan: .planning/phases/25-macaque-germlines-integration/PLAN.md
Status: Complete
Last activity: 2026-01-25 — Completed Phase 25 (macaque germlines)

Progress: █████░░░░░░░░░░░░░░░ 25% (1/4 phases complete)

**Milestone v1.3: Test Infrastructure & Species Expansion** — In Progress (phases 25-28)

**Next Phase:** Phase 26 — Add AIRR Package Dependency

## Milestone v1.3 Overview

**Goal:** Fix skipped tests by adding macaque germlines, airr package dependency, removing deprecated G3 tests, and fix germline priority order

### Phase 25: Macaque Germlines Integration ✓
- ✓ Build macaque IgBLAST databases
- ✓ Generate internal_data and aux files
- ✓ Enable 6 previously skipped tests (4 pass, 2 fail due to pre-existing bug)

### Phase 26: Add AIRR Package Dependency
- [ ] Add airr to pyproject.toml
- [ ] Remove importorskip from test

### Phase 27: Remove Deprecated G3 Tests
- [ ] Review and migrate/remove 2 G3 legacy tests

### Phase 28: Fix Germline Priority Order
- [ ] Update default priority to ['vdjbase', 'ogrdb', 'imgt', 'custom']
- [ ] Document priority rationale

---

## Milestone v1.2 Overview (Complete)

**Goal:** Enable reference.yml to select alleles from all germline sources (imgt, ogrdb, vdjbase, custom), using germlines module as data provider instead of G3 API.

### Phase 19: Source Validation ✓
- ✓ Expand VALID_SOURCES in models.py (imgt, ogrdb, vdjbase, custom)
- SRC-02 moved to Phase 24 (validate source/species in germlines path)

### Phase 20: Integration Foundation ✓
- ✓ Add `use_germlines=True` to `References.from_yaml()`
- ✓ Route explicit source to GermlineManager (providers=[source])
- ✓ Generate synthetic `_id` in adapter (SHA-256 hash)

### Phase 21: Build CLI ✓
- ✓ Add `sadie reference build` command
- ✓ Generate complete IgBLAST database structure
- ✓ Progress output during build

### Phase 22: Runtime Usage ✓
- ✓ Add `Airr(database=<path>)` parameter
- ✓ Skip germlines/G3 lookup with prebuilt
- ✓ Validate database structure on load

### Phase 23: Documentation ✓
- ✓ Create reference-sample.yml (multi-source examples)
- ✓ Document build → use workflow

### Phase 24: Gap Closure ✓
- ✓ Implement SRC-02 source/species validation for germlines
- ✓ Add IMGT region positions so `--use-germlines` build succeeds
- ✓ Create SUMMARY/VERIFICATION artifacts for phases 19-23

## Milestone v1.1 Summary (Complete)

**Final Result:** 98.29% structural parity between germlines and G3 backends

| Phase | Description | Status |
|-------|-------------|--------|
| 13 | Backend Parity Audit | ✓ Complete |
| 14 | C Region Data Integration | ✓ Complete |
| 15 | J Gene Matching & CDR3 Fix | ✓ Complete |
| 16 | Fix NDM.IMGT FWR3 End | ✓ Complete |
| 17 | Fix complete_vdj Quirk | ✓ Complete |
| 18 | Document D-region Variance | ✓ Complete |

## Key Files

### v1.2 Target Files
- `src/sadie/reference/models.py` — Source validation expansion
- `src/sadie/reference/reference.py` — from_yaml() germlines integration
- `src/sadie/germlines/g3_adapter.py` — Add `_id` field generation
- `src/sadie/reference/cli.py` — Build CLI command (new)
- `src/sadie/airr/airr.py` — Database path parameter

### v1.1 Artifacts
- `audit/audit.md` — Detailed audit report
- `audit/parity-notes.md` — Parity explanation
- `audit/igblast-quirk.md` — IgBLAST quirk documentation

## Session Continuity

Last session: 2026-01-25T19:34:16Z
Stopped at: Completed phase-24/PLAN.md
Resume file: None

---
*Last updated: 2026-01-25*
