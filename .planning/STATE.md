# State: Germline Database Integration

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-22)

**Core value:** Enable researchers to select germline database for AIRR annotation and renumbering
**Current focus:** v1.1 Audit — CDR3 annotation investigation

## Current Position

Phase: 14 Complete
Plan: Completed
Status: C gene integration done; CDR3 issue identified
Last activity: 2026-01-22 — Phase 14 completed

Progress: ██████████░░░░░░░░░░ 55% (5/8 requirements)

## Milestone v1.1 Progress

### Phase 13: Backend Parity Audit ✓
- Audit completed: 72.19% parity
- Root cause identified: Missing C gene data
- Documented in `audit/audit.md`

### Phase 14: C Region Data Integration ✓
- Added C gene data from IMGT GENE-DB (684 IGHC, 4 IGKC, 16 IGLC)
- Generated IgBLAST C gene databases (human_C.*)
- No more "C gene directory not found" warnings
- All 129 columns now present
- CDR3 annotation issue identified as pre-existing (separate from C genes)

### Phase 15: CDR3 Annotation Fix ○ (Next)
- Investigate J gene matching failure
- Debug IgBLAST configuration
- Target: CDR3/junction annotation working

## Key Discovery

**CDR3 annotation failure is PRE-EXISTING** (not caused by C gene absence):
- J genes not being matched (j_call = NaN)
- CDR3, junction, fwr4 all return NaN
- Affects 99% of sequences
- Requires IgBLAST configuration investigation

## Key Files

### Phase 14 Completed
- `src/sadie/germlines/scripts/download_imgt.py` — GENE-DB C gene download
- `src/sadie/germlines/igblast/database/human/human_C.*` — C gene BLAST DBs
- `.planning/phases/phase-14/SUMMARY.md` — Phase summary

### Phase 15 Targets
- `src/sadie/airr/igblast/igblast.py` — IgBLAST execution
- `src/sadie/germlines/igblast/aux_db/` — Aux files
- `src/sadie/germlines/igblast/Ig/internal_data/` — Internal data

### Audit Artifacts
- `audit/audit.ipynb` — Comparison notebook
- `audit/audit.md` — Detailed report
- `audit/audit.py` — Quick audit script
- `audit/20260112_HCV_DB_example.csv` — Test data

---
*Last updated: 2026-01-22*
