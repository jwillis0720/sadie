# Phase 25: Macaque Germlines Integration — Context

## Overview

**Phase Goal:** Build macaque germline databases from VDJbase to enable 6 skipped tests.

**Decision Date:** 2026-01-25

---

## Decisions

### 1. Data Source Selection

**Decision:** Use VDJbase "Rhesus Macaque" as the authoritative source.

**Rationale:** VDJbase is the canonical, maintained source for macaque germlines. The custom sources in `sources/custom/macaque/` are legacy data.

**Implications:**
- Fetch germlines from VDJbase API (`/genomic/data_sets/Rhesus%20Macaque`)
- VDJbase provides IGH, IGK, IGL datasets for Rhesus Macaque
- Must include C genes (constant regions) — required, fail if missing

---

### 2. Species Naming Convention

**Decision:** Rename `rhesus_macaque/` → `macaque/` in VDJbase provider mapping.

**Rationale:** 
- `macaque` is the shorthand name used throughout the codebase
- Tests check for `internal_data/macaque/` directory
- Other databases use `macaque` consistently
- Creating an alias might break other database mappings

**Implementation:**
- Update VDJbase provider species mapping: `"Rhesus Macaque": "macaque"` (not `rhesus_macaque`)
- Build databases under `macaque/` directory names

---

### 3. Custom Sources Handling

**Decision:** Keep existing custom sources as legacy data.

**Rationale:** The data in `sources/custom/macaque/` (~4000 V genes) is valuable for future curation. It may contain internal lab data that should be preserved.

**Implementation:**
- Do not delete `sources/custom/macaque/`
- VDJbase is primary; custom remains for potential future use
- No integration of custom sources in this phase

---

### 4. Test Fixture Strategy

**Decision:** Regenerate fixtures from VDJbase-annotated data.

**Rationale:** 
- Tests should verify annotation quality (specific gene calls), not just functionality
- Fixtures were created for an older germline set; VDJbase may have different allele names
- Existing regeneration scripts already handle this workflow

**Implementation:**
- After building VDJbase macaque databases, run:
  - `python scripts/regenerate_igl_reference.py`
  - `python scripts/regenerate_linked_igl_reference.py`
- Update expected gene calls in fixtures to match VDJbase annotations

---

### 5. C Gene Requirements

**Decision:** VDJbase must provide C genes; fail if missing.

**Rationale:**
- `test_airr_constant_region_macaque` requires C gene annotation
- Each germline source is required to have complete V/D/J/C coverage
- The test uses human sequences (`HD_w_constant.fasta`) to verify C gene annotation capability

**Implementation:**
- Verify VDJbase macaque includes C genes before proceeding
- If C genes missing from VDJbase, phase cannot complete (escalate)

---

## Constraints

- **Tests expect `macaque` species name** — `Airr("macaque")` must resolve
- **Tests check specific path** — `internal_data/macaque/` must exist
- **6 tests must pass** — Not just run, but pass with correct assertions

---

## Out of Scope (Deferred)

### Provider Priority Reordering

**Idea:** Reorder provider priority to `['vdjbase', 'ogrdb', 'imgt', 'custom']`

**Rationale:** VDJbase is best for human and macaque; OGRDB is good for mouse; IMGT provides species diversity; custom fills gaps.

**Status:** Captured for future phase. Not part of Phase 25.

---

## Dependencies

- VDJbase API availability
- VDJbase macaque data includes V, D, J, and C genes
- `makeblastdb` (BLAST+) installed for database generation

---

## Success Criteria

1. `GermlineData("macaque")` resolves without error
2. All 6 macaque tests pass (not skipped)
3. Macaque annotation produces valid AIRR output with correct gene calls
4. C gene annotation works for macaque species

---

## Files to Modify

**VDJbase Provider:**
- `src/sadie/germlines/providers/vdjbase.py` — Change species mapping to `macaque`

**Database Generation:**
- `src/sadie/germlines/igblast/database/macaque/` — BLAST databases
- `src/sadie/germlines/igblast/Ig/internal_data/macaque/` — NDM files + symlinks
- `src/sadie/germlines/igblast/aux_db/macaque_gl.aux` — Auxiliary file

**Tests:**
- `tests/unit/airr/test_airr.py` — Remove `@skip_no_macaque` decorators
- `tests/unit/airr/test_methods.py` — Remove `@skip_no_macaque` decorators

**Fixtures (regenerate):**
- `tests/data/fixtures/airr_tables/igl_out.feather`
- `tests/data/fixtures/airr_tables/bum_link_solution.feather`
