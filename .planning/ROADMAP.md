# Roadmap: Germline Database Integration

**Milestone:** v1.1 Audit
**Phases:** 13 (continuing from v1.0)

## Phase 13: Backend Parity Audit

**Goal:** Validate germlines backend produces identical AIRR results to G3 backend

**Requirements:**
- AUDIT-01: Run AIRR annotation with germlines backend
- AUDIT-02: Run AIRR annotation with G3 backend
- AUDIT-03: Compare results for column-level identity
- AUDIT-04: Document discrepancies with root cause analysis

**Success Criteria:**
1. Audit notebook executes without errors
2. Both backends process all 95 test sequences
3. Column comparison completes and reports parity percentage
4. Any differences are documented with explanation

**Deliverables:**
- `audit/audit.ipynb` — Comparison notebook
- `audit/20260112_HCV_DB_example.csv` — Test data
- Documented results in notebook output

---
*Created: 2026-01-22*
