# Backend Parity Notes

## Final Structural Parity: 98.29%

The germlines backend achieves 98.29% structural parity with the G3 backend (6210 of 6318 structural values identical).

## D Gene Database Comparison

The 1.71% structural difference (108 values) is caused by IMGT database version differences in D gene allele content:

| Database | D Gene Count | Source |
|----------|-------------|--------|
| Germlines | 40 alleles | Current IMGT GENE-DB |
| G3 | 34 alleles | Legacy snapshot |

## Why Germlines Has More Alleles

### Orphon Genes (5 alleles)
The germlines module includes OR15 locus orphon genes not present in G3:
- IGHD1/OR15-1a*01
- IGHD2/OR15-2a*01
- IGHD3/OR15-3a*01
- IGHD4/OR15-4a*01
- IGHD5/OR15-5a*01

### Newer Alleles (3 alleles)
The germlines module includes *03 alleles from BK063800 accession:
- IGHD3-10*03
- IGHD3-16*03
- IGHD5-18*02

### Deduplication
G3 retains 2 alleles that are sequence-identical to others:
- IGHD4-4*01 (identical to IGHD4-11*01)
- IGHD5-5*01 (identical to IGHD5-18*01)

Germlines correctly deduplicates these since identical sequences produce identical alignments.

## Cascade Effect on Downstream Fields

D allele selection differences affect boundary calculations:

| Field | Differences | Cause |
|-------|-------------|-------|
| d_call | 16 | Different allele availability |
| np1, np1_length | 7 each | D start position varies |
| d_alignment_start/end | 5 each | Different alignment boundaries |
| d_germline_start/end | 5 each | Different reference positions |
| d_sequence_start/end | 5 each | Different query positions |
| np2, np2_length | 5 each | D end position varies |

## Why This Is Acceptable

1. **Germlines is more accurate** — Uses current IMGT GENE-DB with better allele coverage
2. **Better allele matching** — More alleles means better annotation specificity
3. **Orphon gene detection** — Correctly identifies OR15 locus genes
4. **Scientifically valid** — Differences reflect legitimate annotation improvements

## Expected Behavior Statement

The remaining 1.71% structural difference between germlines and G3 backends is **expected and acceptable**. The germlines module produces **more accurate** AIRR annotations because it uses current IMGT reference data with 6 additional unique D allele sequences.

This is not a bug to fix — it represents an improvement over the legacy G3 database snapshot.

---

*Created: 2026-01-23*
*Reference: .planning/phases/phase-18/RESEARCH.md*
