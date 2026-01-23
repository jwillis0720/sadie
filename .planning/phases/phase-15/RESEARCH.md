# Phase 15: Parity Gap Analysis — Research

**Created:** 2026-01-22  
**Updated:** 2026-01-22  
**Status:** Research Complete (Extended Analysis)

---

## Summary

Phase 15 fixed J gene matching and improved parity from 72.19% → 77.60%. However, a ~22% gap remains. This document analyzes the root causes of all remaining differences.

**Primary Root Cause Identified:** The `ndm.imgt` file generator uses full V gene sequence length instead of the IMGT-standard FWR3 end position (gapped position 312).

**Secondary Issues:**
- Different J gene allele databases
- Additional V gene alleles in germlines database

---

## Detailed Findings

### 1. FWR3/CDR3 Boundary Issue (CRITICAL - 100% of sequences)

**Root Cause:** The `build_internal_data.py` script generates ndm.imgt files with incorrect column 11 values.

**Technical Details:**
- Column 11 of ndm.imgt specifies the V gene endpoint (where FWR3 ends, CDR3 begins)
- G3 uses ungapped position of IMGT gapped position 312 (standard V gene endpoint)
- Germlines uses full ungapped sequence length

**Example for IGHV1-69*01:**
```
G3:        288 (ungapped position of IMGT position 312, the conserved Cys)
Germlines: 296 (full ungapped sequence length)
Difference: 8 nucleotides
```

**Position Mapping:**
```
IMGT gapped position 312 → ungapped position 288
Nucleotides 286-288: TGT (Cysteine codon - end of V gene proper)
Nucleotides 289-296: GCGAGAGA (3' sequence after Cys-104)
```

**Columns Affected:**
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

**Code Location:** `src/sadie/germlines/scripts/build_internal_data.py`

**Current (Wrong):**
```python
entry = (
    f"{gene_name}\t"
    ...
    f"{regions['FR3'][0]}\t{seq_len}\t"  # seq_len = full ungapped length
    f"{chain_type}\t0"
)
```

**Fix Required:**
```python
# Calculate ungapped position of IMGT gapped position 312 (FR3 end)
fr3_end_ungapped = calculate_ungapped_position(gapped_seq, 312)
entry = (
    f"{gene_name}\t"
    ...
    f"{regions['FR3'][0]}\t{fr3_end_ungapped}\t"  # Use FR3 end, not seq_len
    f"{chain_type}\t0"
)
```

---

### 2. J Gene Allele Differences (72/78 sequences = 92%)

**Root Cause:** Different J gene allele sets between backends.

**Examples:**
| Germlines | G3 |
|-----------|-----|
| IGHJ3*01 | IGHJ3*02 |
| IGHJ4*01 | IGHJ4*02 |
| IGHJ5*01,IGHJ5*04 | IGHJ5*02 |
| IGHJ6*01,IGHJ6*04 | IGHJ6*02 |

**Analysis:** The germlines database has different J gene alleles available, leading to different best-match alleles and tied allele lists.

**Columns Affected:**
| Column | Diff Count | Impact |
|--------|-----------|--------|
| j_call | 72/78 | Different allele assignments |
| j_call_top | 64/78 | Different top allele |
| j_germline_alignment | 62/78 | Different reference sequence |
| j_identity | 62/78 | Different identity calculation |
| j_score | 60/78 | Different BLAST scores |
| j_cigar | 32/78 | Different alignment |

**Severity:** LOW - Gene-level calls match; only allele-level differs.

---

### 3. V Gene Allele Differences (40/78 sequences = 51%)

**Root Cause:** Germlines has MORE alleles (707 vs 399 in ndm.imgt).

**Allele Count Comparison:**
| Backend | V Gene Alleles in ndm.imgt |
|---------|---------------------------|
| G3 | 399 |
| Germlines | 707 |

**Types of Differences:**
- **Same gene, different allele/order:** 38/40 (95%)
- **Different gene:** 2/40 (5%)

**Examples:**
| Germlines | G3 | Category |
|-----------|-----|----------|
| IGHV1-69*01,IGHV1-69*17_a244g | IGHV1-69*01 | Extra variant allele |
| IGHV4-34*01_a318g | IGHV4-34*01 | Variant as primary |
| IGHV3-48*01,IGHV3-48*01_a85c,IGHV3-48*02 | IGHV3-48*01,IGHV3-48*02 | Extra variant allele |
| IGHV3-20*01 | IGHV3-20*04 | Different primary allele |

**Severity:** LOW - Most are expected differences due to larger allele database.

---

### 4. Score/Support Differences (All sequences affected)

**Root Cause:** BLAST scores and E-values depend on database size.

**Score Analysis:**
| Metric | Mean Diff | Max Diff | Non-matching |
|--------|-----------|----------|--------------|
| v_score | +1.08 | 14.02 | 18/78 |
| j_score | -5.30 | 42.30 | 60/78 |
| d_score | +0.19 | 3.17 | 5/78 |
| v_identity | +0.0014 | 0.019 | 19/78 |
| j_identity | -0.016 | 0.085 | 62/78 |

**Support (E-values):** All 78 sequences differ because database size affects E-value calculation.

**Severity:** LOW - Expected statistical variation, not indicative of incorrect results.

---

## Priority Ranking

### Priority 1: Fix ndm.imgt FWR3 End Position (CRITICAL)

**Impact:** 100% of sequences affected  
**Effort:** Low (single file change)  
**Fix:** Modify `build_internal_data.py` to calculate ungapped position of IMGT position 312

**Expected Improvement:** ~15-20% parity increase
- fwr3, fwr3_aa, fwr3_end will match
- cdr3, cdr3_aa, cdr3_start will match
- junction fields will match

### Priority 2: Align J Gene Allele Set (HIGH)

**Impact:** 92% of sequences affected  
**Effort:** Medium (database synchronization)  
**Options:**
1. Use G3's J gene alleles for germlines
2. Accept allele-level differences (gene-level matches)

**Expected Improvement:** If aligned, ~5-8% parity increase

### Priority 3: V Gene Allele Consistency (MEDIUM)

**Impact:** 51% of sequences affected  
**Effort:** Medium-High  
**Note:** Most differences are due to germlines having MORE alleles, which is intentional for the richer database.

**Recommendation:** Accept these differences as expected behavior.

### Priority 4: Score/Support Normalization (LOW)

**Impact:** 100% of sequences affected  
**Effort:** High (would require database size normalization)  
**Recommendation:** Do not fix - expected statistical variation.

---

## Validation Checklist (After Fixes)

After implementing Priority 1 fix:

- [ ] Rebuild human ndm.imgt file
- [ ] Verify IGHV1-69*01 column 11 = 288 (not 296)
- [ ] Run audit and verify:
  - [ ] fwr3_end differences: 0
  - [ ] cdr3_start differences: 0
  - [ ] CDR3 length matches G3
  - [ ] Parity >= 90%

---

## Technical Details

### NDM.IMGT File Format

```
Column 1:  Gene name (e.g., IGHV1-69*01)
Column 2:  FR1 start
Column 3:  FR1 end
Column 4:  CDR1 start
Column 5:  CDR1 end
Column 6:  FR2 start
Column 7:  FR2 end
Column 8:  CDR2 start
Column 9:  CDR2 end
Column 10: FR3 start
Column 11: V gene end (IMGT position 312 ungapped) <-- THE BUG
Column 12: Chain type (VH, VK, VL)
Column 13: Flag (0)
```

### IMGT Position Mapping

IMGT uses a standardized gapping system:
- Complete V gene: positions 1-312 (gapped)
- FR3 ends at position 312 (Cys-104 codon)
- For IGHV1-69*01: gapped 312 → ungapped 288

The gapped sequence has insertions at CDR1 and CDR2, but positions 1-312 always correspond to the same structural elements.

---

## Sources

| Source | Confidence | Notes |
|--------|------------|-------|
| G3 ndm.imgt file | HIGH | Reference implementation |
| Audit output (audit/audit.py) | HIGH | Direct comparison data |
| build_internal_data.py analysis | HIGH | Identified exact bug |
| IMGT numbering documentation | HIGH | Standard position definitions |

---

## Files to Modify

| File | Change |
|------|--------|
| `src/sadie/germlines/scripts/build_internal_data.py` | Fix column 11 calculation |
| `src/sadie/germlines/igblast/Ig/internal_data/human/human.ndm.imgt` | Regenerate |

---

*Research complete. Ready for Phase 16 planning.*
