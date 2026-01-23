# State: Germline Database Integration

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-22)

**Core value:** Enable researchers to select germline database for AIRR annotation and renumbering
**Current focus:** v1.2 Reference Module Unification

## Current Position

Phase: 22 — Runtime Usage
Plan: PLAN.md
Status: Complete
Last activity: 2026-01-23 — Phase 22 implemented

Progress: ████████████████░░░░ 80%

**Milestone v1.2: Reference Module Unification** — In progress (phases 19-23)

## Milestone v1.2 Overview

**Goal:** Enable reference.yml to select alleles from all germline sources (imgt, ogrdb, vdjbase, custom), using germlines module as data provider instead of G3 API.

### Phase 19: Source Validation ✓
- ✓ Expand VALID_SOURCES in models.py (imgt, ogrdb, vdjbase, custom)
- SRC-02 deferred to Phase 20 (validate at add_genes time)

### Phase 20: Integration Foundation ✓
- ✓ Add `use_germlines=True` to `References.from_yaml()`
- ✓ Route explicit source to GermlineManager (providers=[source])
- ✓ Generate synthetic `_id` in adapter (SHA-256 hash)

### Phase 21: Build CLI ✓
- ✓ Add `sadie reference build` command
- ✓ Generate complete IgBLAST database structure
- ✓ Progress output during build
- Note: --use-germlines has gap (missing IMGT region fields)

### Phase 22: Runtime Usage ✓
- ✓ Add `Airr(database=<path>)` parameter
- ✓ Skip germlines/G3 lookup with prebuilt
- ✓ Validate database structure on load

### Phase 23: Documentation
- Create reference-sample.yml
- Document build → use workflow

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

---
*Last updated: 2026-01-23*
