# Phase 16: Fix NDM.IMGT FWR3 End Position

## Overview

**Goal:** Fix ndm.imgt file generation to use correct FWR3 end position (IMGT position 312) instead of full V gene sequence length.

**Impact:** Improve backend parity from 77.6% to ~95%+

## Problem Analysis

### Current Behavior
The `build_internal_data.py` script generates ndm.imgt files with column 11 set to `seq_len` (full ungapped sequence length):
```
IGHV1-69*01  1  75  76  99  ...  175  296  VH  0
                                      ^^^
                              seq_len (wrong)
```

### Expected Behavior (G3 Reference)
Column 11 should be the ungapped position corresponding to IMGT gapped position 312 (FR3 end):
```
IGHV1-69*01  1  75  76  99  ...  175  288  VH  0
                                      ^^^
                              FR3 end (correct)
```

### Root Cause
In `generate_ndm_entry()` at line ~97:
```python
entry = (
    f"{gene_name}\t"
    ...
    f"{regions['FR3'][0]}\t{seq_len}\t"  # <-- Using seq_len instead of FR3 end
    f"{chain_type}\t0"
)
```

The code already calculates `regions['FR3'][1]` (FR3 end position) but doesn't use it.

---

## Tasks

### Task 1: Fix Column 11 Calculation
**File:** `src/sadie/germlines/scripts/build_internal_data.py`

**Change:** In `generate_ndm_entry()`, replace `seq_len` with `regions['FR3'][1]` for column 11.

**Before (line ~97):**
```python
f"{regions['FR3'][0]}\t{seq_len}\t"
```

**After:**
```python
f"{regions['FR3'][0]}\t{regions['FR3'][1]}\t"
```

**Validation:**
- Run script to confirm it executes without errors
- Verify output format matches G3 reference

---

### Task 2: Regenerate human.ndm.imgt
**Execute:**
```bash
cd /Users/tmsincomb/sadie/src/sadie/germlines
python scripts/build_internal_data.py human
```

**Validation:**
- Compare IGHV1-69*01 column 11: should be 288 (was 296)
- Compare IGHV1-18*01 column 11: should be 288 (was 296)
- Spot-check additional alleles against G3 reference

---

### Task 3: Verify Backend Parity Improvement
**Execute:**
```bash
cd /Users/tmsincomb/sadie/audit
python audit.py
```

**Expected Results:**
- FWR3 boundaries should now match G3
- CDR3 start positions should align
- Overall backend parity should improve from 77.6% toward 95%+

---

## Success Criteria

| Metric | Before | After |
|--------|--------|-------|
| IGHV1-69*01 col 11 | 296 | 288 |
| FWR3 boundary match | ~0% | ~95%+ |
| CDR3 start alignment | ~1% | ~95%+ |
| Backend parity | 77.6% | ~95%+ |

## Dependencies

- None (self-contained fix)

## Risks

| Risk | Mitigation |
|------|------------|
| Incomplete V genes without FR3 | Existing code already handles this with `required` check |
| Edge cases with truncated sequences | FR3 end will reflect actual sequence end if truncated before position 312 |

## Execution Notes

1. The fix is a single-line change
2. Regeneration is idempotent - safe to run multiple times
3. Audit script provides immediate validation feedback
