# State: Germline Database Integration

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-22)

**Core value:** Enable researchers to select germline database for AIRR annotation and renumbering
**Current focus:** Phase 31 — Add Database Parameter Support to Renumbering

## Current Position

Phase: 31 — Add Database Parameter Support to Renumbering
Plan: .planning/phases/31-add-database-parameter-support-to-renumbering/31-01-PLAN.md
Status: Complete
Last activity: 2026-01-27 — Completed Plan 31-01 (HMM Building)

Progress: ████████████████████ 100% (1/? plans complete)

**Phase 31-01: Add HMM Building to Reference Database Build** — COMPLETE

### Phase 31-01: Add HMM Building to Reference Database Build ✓
- ✓ Add _make_hmm_files() method to References class
- ✓ Add _write_stockholm_file() helper for Stockholm format
- ✓ Add _translate_gapped_nt_to_aa() for NT-to-AA fallback
- ✓ Integrate _make_hmm_files() into make_airr_database()
- ✓ Add 3 unit tests for HMM building functionality

### Phase 30: Add _gapped.fasta Support to CustomProvider ✓
- ✓ Add _get_gapped_fasta_path() and _load_gapped_sequences() methods
- ✓ Modify fetch_genes() to load gapped sequences from _gapped.fasta
- ✓ Update _create_gene_from_record() to use pre-loaded gapped sequences
- ✓ Add 3 unit tests for _gapped.fasta support

### Phase 29: Germline Source Tracking ✓
- ✓ Add get_source_lookup() to GermlineData
- ✓ Add _lookup_source() and _add_source_columns() to Airr
- ✓ Integrate source columns into run_fasta() and _run_scfv()
- ✓ Add 4 unit tests for source tracking

**Next Phase:** Phase 31 — Add Database Parameter Support to Renumbering

## Milestone v1.3 Overview

**Goal:** Fix skipped tests by adding macaque germlines, airr package dependency, removing deprecated G3 tests, and fix germline priority order

### Phase 25: Macaque Germlines Integration ✓
- ✓ Build macaque IgBLAST databases
- ✓ Generate internal_data and aux files
- ✓ Enable 6 previously skipped tests (4 pass, 2 fail due to pre-existing bug)

### Phase 26: Add AIRR Package Dependency ✓
- ✓ Add airr to pyproject.toml
- ✓ Remove importorskip from test

### Phase 27: Remove Deprecated G3 Tests ✓
- ✓ Remove TestGermlineDataPaths class (deprecated G3 API tests)
- ✓ Create G3 deprecation documentation
- ✓ Update CONCERNS.md to remove deleted test reference

### Phase 28: Fix Germline Priority Order ✓
- ✓ Update default priority to ['vdjbase', 'ogrdb', 'imgt', 'custom']
- ✓ Document priority rationale

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

### v1.3 Target Files (Phase 31)
- `src/sadie/reference/reference.py` — HMM building methods added
- `src/sadie/germlines/renumbering_integration.py` — LocalHMMBuilder reference implementation

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

## Accumulated Context

### Roadmap Evolution
- Phase 30 complete: Add _gapped.fasta support to CustomProvider
- Phase 31 added: Add Database Parameter Support to Renumbering

### Pending Todos
- 1 pending todo(s) in `.planning/todos/pending/`

## Session Continuity

Last session: 2026-01-27
Stopped at: Completed Phase 31-01 — HMM building added to reference database build
Resume file: None

---
*Last updated: 2026-01-27*
