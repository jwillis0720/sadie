# State: Germline Database Integration

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-22)

**Core value:** Enable researchers to select germline database for AIRR annotation and renumbering
**Current focus:** v1.1 Audit — validate germlines backend parity with G3

## Current Position

Phase: 13
Plan: Not started
Status: Defining requirements
Last activity: 2026-01-22 — Milestone v1.1 started

Progress: ░░░░░░░░░░░░░░░░░░░░ 0% (0/4 requirements)

## Milestone v1.1 Goals

- Validate germlines backend produces identical AIRR results to G3
- Create audit notebook with side-by-side comparison
- Document any discrepancies with root cause analysis

## Test Data

- `audit/20260112_HCV_DB_example.csv` — 95 HCV antibody sequences

## Key Files

### v1.1 Audit
- `audit/audit.ipynb` — Comparison notebook (to create)

### Integration Points (from v1.0)
- `src/sadie/airr/igblast/germline.py` — IgBLAST paths
- `src/sadie/renumbering/aligners/hmmer.py` — HMM generation
- `src/sadie/reference/reference.py` — Reference system
- `src/sadie/germlines/cli.py` — CLI commands

---
*Last updated: 2026-01-22*
