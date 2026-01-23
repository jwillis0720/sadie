# Phase 18 Summary: Document D-region IMGT Version Variance

## Decision

**Accepted 98.29% structural parity as final.**

No code changes required. The 1.71% structural difference is caused by IMGT database version differences, not a bug.

## Rationale

The germlines module uses current IMGT GENE-DB data (40 D alleles) while G3 uses a legacy snapshot (34 D alleles). This difference is **expected behavior** — the germlines module produces more accurate results.

## Key Findings

| Aspect | Germlines | G3 |
|--------|-----------|-----|
| D gene alleles | 40 | 34 |
| OR15 orphon genes | ✓ (5 alleles) | ✗ |
| Newer *03 alleles | ✓ (3 alleles) | ✗ |
| Sequence deduplication | ✓ | ✗ |

### Root Cause Analysis

1. **Different database content** — Germlines has 6 more unique allele sequences
2. **Cascade effect** — D allele differences affect np1/np2 and boundary calculations
3. **Legitimate variance** — Not a bug, reflects database version differences

## Files Created/Modified

| Action | File | Description |
|--------|------|-------------|
| CREATE | `audit/parity-notes.md` | Parity explanation documentation |
| MODIFY | `.planning/ROADMAP.md` | Phase 18 marked complete |
| MODIFY | `.planning/STATE.md` | Updated to Phase 18 complete |
| CREATE | `.planning/phases/phase-18/SUMMARY.md` | This summary |

## Detailed Research

See [RESEARCH.md](./RESEARCH.md) for complete analysis including:
- D gene database comparison
- Allele-by-allele differences
- Cascade effect on downstream fields
- Decision rationale

## Milestone v1.1 Complete

Phase 18 concludes the v1.1 Audit milestone:

| Phase | Goal | Status |
|-------|------|--------|
| 13 | Backend Parity Audit | ✓ Complete |
| 14 | C Region Data Integration | ✓ Complete |
| 15 | J Gene Matching & CDR3 Annotation Fix | ✓ Complete |
| 16 | Fix NDM.IMGT FWR3 End Position | ✓ Complete |
| 17 | Fix complete_vdj IgBLAST Quirk | ✓ Complete |
| 18 | Document D-region IMGT Version Variance | ✓ Complete |

**Final Result:** 98.29% structural parity with scientifically valid remaining differences.

---

*Completed: 2026-01-23*
