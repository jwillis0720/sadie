# Codebase Concerns

## Overview

This document tracks technical debt, issues, and concerns in the SADIE codebase, with particular focus on segment position discovery differences between the germlines module and G3 backend.

**Current Backend Parity:** 77.60% (improved from 72.19% after Phase 15 J gene fix)

---

## Critical: NDM.IMGT FWR3 End Position Bug

### Problem Statement

The `build_internal_data.py` script generates ndm.imgt files with **incorrect column 11 values** (V gene endpoint). This causes 100% of sequences to have wrong FWR3/CDR3 boundaries.

### Root Cause

**Current (Wrong):** Uses full ungapped sequence length
```python
# src/sadie/germlines/scripts/build_internal_data.py line 97
entry = (
    f"{gene_name}\t"
    ...
    f"{regions['FR3'][0]}\t{seq_len}\t"  # seq_len = full ungapped length
    f"{chain_type}\t0"
)
```

**Expected:** Should use ungapped position of IMGT gapped position 312 (FR3 end)

### Example: IGHV1-69*01

| Backend | Column 11 Value | Meaning |
|---------|-----------------|---------|
| G3 | 288 | Ungapped position of IMGT position 312 (Cys-104 codon) |
| Germlines | 296 | Full ungapped sequence length |
| Difference | 8 nucleotides | Extra 3' sequence after Cys-104 |

### Impact

| Column | Diff Count | Cause |
|--------|-----------|-------|
| fwr3 | 78/78 (100%) | FWR3 is 8nt longer in germlines |
| fwr3_aa | 78/78 (100%) | ~3 extra amino acids |
| fwr3_end | 78/78 (100%) | 8 position offset |
| cdr3 | 77/78 (99%) | CDR3 is 8nt shorter in germlines |
| cdr3_aa | 77/78 (99%) | ~3 fewer amino acids |
| cdr3_start | 77/78 (99%) | 8 position offset |
| junction | 77/78 (99%) | Shorter junction |
| junction_aa | 77/78 (99%) | Shorter junction |

### Required Fix

```python
# Calculate ungapped position of IMGT gapped position 312 (FR3 end)
def get_fr3_end_ungapped(gapped_seq: str) -> int:
    """Calculate ungapped position of IMGT position 312."""
    ungapped_pos = 0
    for gapped_pos, char in enumerate(gapped_seq, 1):
        if char not in ".-":
            ungapped_pos += 1
        if gapped_pos == 312:
            return ungapped_pos
    return ungapped_pos  # Sequence shorter than 312

# In generate_ndm_entry():
fr3_end = get_fr3_end_ungapped(gapped_seq)
entry = f"...{regions['FR3'][0]}\t{fr3_end}\t{chain_type}\t0"
```

### Files to Modify

| File | Change |
|------|--------|
| `src/sadie/germlines/scripts/build_internal_data.py` | Fix column 11 calculation |
| `src/sadie/germlines/igblast/Ig/internal_data/human/human.ndm.imgt` | Regenerate |

### Expected Improvement

~15-20% parity increase (fwr3, cdr3, junction fields will match)

---

## High: J Gene Allele Database Differences

### Problem Statement

The germlines and G3 backends have **different J gene allele sets**, causing 92% of sequences to get different j_call assignments.

### Root Cause

**Aux file allele comparison:**

| Germlines Aux | G3 Aux |
|---------------|--------|
| IGHJ3*01 | IGHJ3*02 |
| IGHJ4*01, IGHJ4*03 | IGHJ4*02 |
| IGHJ5*01, IGHJ5*03, IGHJ5*04 | IGHJ5*02 |
| IGHJ6*01, IGHJ6*04 | IGHJ6*02 |
| IGLJ2A*01 | (missing) |
| IGLJ3*01 | (missing in germlines) |

**File locations:**
- Germlines: `src/sadie/germlines/igblast/aux_db/human_gl.aux`
- G3: `src/sadie/airr/data/germlines/aux_db/imgt/human_gl.aux`

### Impact

| Column | Diff Count | Impact |
|--------|-----------|--------|
| j_call | 72/78 (92%) | Different allele assignments |
| j_call_top | 64/78 | Different top allele |
| j_germline_alignment | 62/78 | Different reference sequence |
| j_identity | 62/78 | Different identity calculation |
| j_score | 60/78 | Different BLAST scores |
| j_cigar | 32/78 | Different alignment |

### Severity

**LOW** - Gene-level calls match; only allele-level differs. This is expected behavior when using a richer allele database.

### Options

1. **Accept differences** (recommended): Gene-level calls are correct
2. **Synchronize J alleles**: Use G3's J gene alleles for strict parity

---

## High: V Gene Allele Database Differences

### Problem Statement

Germlines has **MORE V gene alleles** (707 vs 399 in ndm.imgt), causing different tie-breaking and variant allele assignments.

### Root Cause

Germlines pulls from multiple sources (IMGT, OGRDB, VDJbase) while G3 uses only IMGT base alleles.

### Examples

| Germlines | G3 | Category |
|-----------|-----|----------|
| IGHV1-69*01,IGHV1-69*17_a244g | IGHV1-69*01 | Extra variant allele |
| IGHV4-34*01_a318g | IGHV4-34*01 | Variant as primary |
| IGHV3-48*01,IGHV3-48*01_a85c,IGHV3-48*02 | IGHV3-48*01,IGHV3-48*02 | Extra variant |
| IGHV3-20*01 | IGHV3-20*04 | Different primary allele |

### Impact

| Column | Diff Count |
|--------|-----------|
| v_call | 40/78 (51%) |
| v_identity | 19/78 |
| v_score | 18/78 |

### Severity

**LOW** - This is intentional behavior. The germlines module has a richer allele database, which provides more specific annotations.

---

## Medium: Score/Support Statistical Differences

### Problem Statement

BLAST scores and E-values depend on database size, causing 100% of sequences to have different support values.

### Analysis

| Metric | Mean Diff | Max Diff | Non-matching |
|--------|-----------|----------|--------------|
| v_score | +1.08 | 14.02 | 18/78 |
| j_score | -5.30 | 42.30 | 60/78 |
| d_score | +0.19 | 3.17 | 5/78 |
| v_identity | +0.0014 | 0.019 | 19/78 |
| j_identity | -0.016 | 0.085 | 62/78 |

**Support (E-values):** All 78 sequences differ because database size affects E-value calculation.

### Severity

**LOW** - Expected statistical variation, not indicative of incorrect results.

---

## Medium: IMGT Provider Missing Gapped Sequence Loading

### Problem Statement

**Rabbit and chicken HMM generation fails** with "no genes have gapped sequences" error, despite IMGT gapped files existing.

### Root Cause

The **IMGT provider** (`src/sadie/germlines/providers/imgt.py`) reads **only** `IGHV.fasta` (ungapped) and does NOT read `IGHV_gapped.fasta`.

**Verified via debugging:**
```python
provider.fetch_genes("rabbit", "V", "H")
# Found 49 genes
# Genes with sequence_gapped: 0/49  <-- ALL MISSING GAPPED DATA!
```

**Why human works but rabbit doesn't:**
- Human's `IGHV.fasta` contains dots (IMGT-gapped format)
- Rabbit's `IGHV.fasta` does NOT contain dots (ungapped)
- Both have separate `IGHV_gapped.fasta` files with proper gaps

### File Structure Evidence

```
sources/imgt/rabbit/
├── IGHV.fasta         # Ungapped (18KB)
├── IGHV_gapped.fasta  # Gapped (20KB)  <-- NOT BEING LOADED
├── IGHJ.fasta         # Ungapped
├── IGHJ_gapped.fasta  # Gapped  <-- NOT BEING LOADED
...
```

### Required Fix

Update `src/sadie/germlines/providers/imgt.py` to match OGRDB pattern:

1. Add `_get_gapped_fasta_path()` method
2. Add `_load_gapped_sequences()` method  
3. In `fetch_genes()`, load gapped sequences and merge with main data
4. Update `_create_imgt_gene()` to accept gapped sequence parameter

### Affected Species

Species with separate gapped files that would benefit:
- rabbit (7 gapped files)
- chicken (5 gapped files)

---

## Low: Biopython Deprecation Warnings

### Problem Statement

Deprecated parameter names in Biopython alignment functions.

### Warnings

- `target_end_gap_score` renamed to `end_insertion_score`
- `query_end_gap_score` renamed to `end_deletion_score`

### Severity

**LOW** - Functional but will break in future Biopython versions.

---

## Low: Pandas Deprecation

### Problem Statement

Passing BlockManager to AirrTable is deprecated.

### Severity

**LOW** - Functional but will break in future pandas versions.

---

## Summary: Parity Gap Analysis

### Current State (77.60% parity)

| Issue | Impact | Fix Priority |
|-------|--------|-------------|
| NDM.IMGT FWR3 end bug | 100% sequences | **P1 - CRITICAL** |
| J gene allele differences | 92% sequences | P3 - Accept |
| V gene allele differences | 51% sequences | P3 - Accept |
| Score/support differences | 100% sequences | P3 - Accept |

### Expected Parity After Fixes

| Fix Applied | Expected Parity |
|-------------|-----------------|
| Current | 77.60% |
| + NDM.IMGT fix | ~92-95% |
| + J allele sync (optional) | ~98-99% |

### Files Summary

**Critical to fix:**
- `src/sadie/germlines/scripts/build_internal_data.py` - Column 11 calculation

**Reference data to regenerate:**
- `src/sadie/germlines/igblast/Ig/internal_data/human/human.ndm.imgt`

**Optional sync:**
- `src/sadie/germlines/igblast/aux_db/human_gl.aux`
- `src/sadie/germlines/builders/j_gene_data.py`

---
*Last updated: 2026-01-22*
*Phase 15 complete: J gene aux format fixed (0% → 100% j_call)*
