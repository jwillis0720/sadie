# State: Germline Database Integration

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-22)

**Core value:** Enable researchers to select germline database for AIRR annotation and renumbering
**Current focus:** v1.1 Audit — Backend parity improvement

## Current Position

Phase: 18 Complete
Plan: Completed
Status: D-region IMGT version variance documented; 98.29% structural parity finalized
Last activity: 2026-01-23 — Phase 18 documentation completed

Progress: ████████████████████ 100% (12/12 requirements)

**Milestone v1.1 Complete** — All 6 phases (13-18) completed successfully

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

### Phase 15: J Gene Matching & CDR3 Annotation Fix ✓
- Root cause: Aux file J gene format was 3 columns instead of 5
- Fix: Created `j_gene_data.py` module with CDR3 end positions
- Fix: Modified `aux.py` to generate correct 5-column format
- Results: j_call 100%, cdr3 98.7%, junction 98.7%, fwr4 98.7%, complete_vdj 97.4%
- Backend parity improved: 72.19% → 77.60%

### Phase 16: Fix NDM.IMGT FWR3 End Position ✓
- Root cause: Column 11 used seq_len instead of FR3 end position
- Fix: Modified `build_internal_data.py` to use `regions['FR3'][1]`
- Fix: Regenerated `human.ndm.imgt` with correct values
- Results: FWR3/CDR3 boundaries now match G3
- Backend parity improved: 77.60% → 86.71%

### Phase 17: Fix complete_vdj IgBLAST Quirk ✓
- Root cause: IgBLAST complete_vdj varies based on allele selection
- Fix: Added `_recalculate_complete_vdj()` method using AIRR standard definition
- Fix: Added `J_GENE_LENGTHS` dictionary for J gene expected lengths
- Results: complete_vdj differences reduced from 22 to 4
- Remaining 4 differences: germlines=True (correct), G3=False (incorrect)

### Phase 18: Document D-region IMGT Version Variance ✓
- Root cause: IMGT database version differences (germlines: 40 D alleles, G3: 34 D alleles)
- Decision: Accept 98.29% structural parity as final
- Documentation: Created `audit/parity-notes.md` explaining IMGT version variance
- Result: Germlines produces MORE accurate annotations due to current IMGT data

## Key Achievements (Phase 15)

| Metric | Before | After | Target |
|--------|--------|-------|--------|
| j_call populated | 0% | 100% | >95% ✓ |
| cdr3 populated | 0% | 98.7% | >95% ✓ |
| junction populated | 0% | 98.7% | >95% ✓ |
| fwr4 populated | 0% | 98.7% | >95% ✓ |
| complete_vdj = True | 0% | 97.4% | >95% ✓ |

## Key Files

### Phase 15 Completed
- `src/sadie/germlines/builders/j_gene_data.py` — J gene reference data module
- `src/sadie/germlines/builders/aux.py` — Fixed aux file builder
- `src/sadie/germlines/igblast/aux_db/human_gl.aux` — Corrected aux file
- `audit/validate_j_gene_fix.py` — Validation script
- `.planning/phases/phase-15/SUMMARY.md` — Phase summary

### Audit Artifacts
- `audit/audit.ipynb` — Comparison notebook
- `audit/audit.md` — Detailed report
- `audit/audit.py` — Quick audit script
- `audit/20260112_HCV_DB_example.csv` — Test data
- `audit/igblast-quirk.md` — IgBLAST complete_vdj quirk documentation

## Roadmap Evolution

- Phase 17 completed: Fix complete_vdj IgBLAST Quirk (22→4 differences, now AIRR-compliant)
- Phase 18 added: Fix D-region Boundary Alignment (target: 98.29% → 99%+ structural parity)

---
*Last updated: 2026-01-22*
