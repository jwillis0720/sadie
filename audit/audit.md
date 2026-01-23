# AIRR Backend Parity Audit Report

**Date**: 2026-01-22  
**Dataset**: `20260112_HCV_DB_example.csv` (96 sequences)  
**Comparison**: Germlines module vs G3 backend

## Executive Summary

The germlines backend achieves **72.19% parity** with the G3 backend. The primary cause is **missing C gene (constant region) reference data** in the germlines module, which cascades into failures in CDR3/junction annotation.

## Critical Issue: Missing C Gene Data

### Symptoms

1. **Warning during execution**:
   ```
   UserWarning: C gene directory not found for human
   UserWarning: /Users/tmsincomb/sadie/src/sadie/germlines/igblast/Ig/internal_data/human/human_C is not found, No C gene segment
   ```

2. **15 columns missing** in germlines output (all C gene related):
   - `c_call`, `c_cigar`, `c_score`, `c_support`, `c_identity`
   - `c_alignment_start`, `c_alignment_end`
   - `c_sequence_start`, `c_sequence_end`
   - `c_germline_start`, `c_germline_end`
   - `c_sequence_alignment`, `c_sequence_alignment_aa`
   - `c_germline_alignment`, `c_germline_alignment_aa`

3. **CDR3/Junction annotation failure**: 95/96 sequences have `NaN` for CDR3 fields

### Root Cause

The germlines module directory structure is missing C gene blast databases:

**Germlines module path** (`src/sadie/germlines/igblast/Ig/internal_data/human/`):
```
human_V.fasta, human_V.n*   ✓ Present
human_D.fasta, human_D.n*   ✓ Present
human_J.fasta, human_J.n*   ✓ Present
human.ndm.imgt              ✓ Present
human_C.*                   ✗ MISSING
```

**C gene data exists** in the older location (`src/sadie/airr/data/germlines/Ig/blastdb/human/`):
```
human_C.fasta, human_C.nin, human_C.nhr, etc.  ✓ Present
```

The germlines module code looks for C gene data in `internal_data/human/human_C` but the blast databases were never copied/generated there.

## Detailed Differences

### Column Count
| Backend   | Columns |
|-----------|---------|
| Germlines | 114     |
| G3        | 129     |
| Difference| 15 (all C gene) |

### Value Differences in Common Columns

**78 columns** have value differences between backends:

| Column | Differences | Notes |
|--------|-------------|-------|
| `j_support` | 96/96 | All rows differ |
| `v_frameshift` | 96/96 | All rows differ |
| `v_support` | 96/96 | All rows differ |
| `cdr3` | 95/96 | Germlines returns NaN |
| `cdr3_aa` | 95/96 | Germlines returns NaN |
| `junction` | 95/96 | Germlines returns NaN |
| `junction_aa` | 95/96 | Germlines returns NaN |
| `fwr4` | 95/96 | Germlines returns NaN |
| `fwr4_aa` | 95/96 | Germlines returns NaN |
| `complete_vdj` | 56/96 | Germlines=False, G3=True |
| `j_call` | 88/96 | J gene assignment differs |
| `v_call` | 50/96 | V gene assignment differs |
| ... | ... | 66 more columns with differences |

### Sample Discrepancies

**CDR3 annotation**:
```
Sequence: 212-1-1
  Germlines: NaN
  G3:        GCGGGAGTAAGGGAGGGTATGGCAGCAATTAGTGGGAAGAATGCTTTTGATATC

Sequence: hcab1
  Germlines: NaN  
  G3:        GCGAGTGTTACGACGAGACAATGGTTCGGGAGGGGTGATGCTTTTGATCTC
```

**complete_vdj flag**:
```
Sequence: 212-1-1
  Germlines: False
  G3:        True
```

## Impact Assessment

| Severity | Impact |
|----------|--------|
| **Critical** | CDR3/junction regions not annotated (95% failure rate) |
| **Critical** | Constant region annotation completely absent |
| **High** | `complete_vdj` incorrectly marked False for productive sequences |
| **Medium** | V/J gene calls differ in ~50-88 sequences |
| **Medium** | Support/identity scores differ (format or calculation differences) |

## Recommended Fix

### Option 1: Copy C gene data (Quick fix)
Copy the existing C gene blast databases to the germlines module location:
```bash
cp src/sadie/airr/data/germlines/Ig/blastdb/human/human_C.* \
   src/sadie/germlines/igblast/Ig/internal_data/human/
```

### Option 2: Update path resolution (Proper fix)
Modify the germlines module code to find C gene data in its existing location, or update the data generation pipeline to include C genes.

### Files to investigate:
- `src/sadie/airr/igblast/germline.py` (line 221 - warns about missing C gene)
- `src/sadie/airr/igblast/igblast.py` (line 625 - path resolution)

## Test Verification

After fix, re-run audit notebook to verify:
1. No C gene warnings
2. All 129 columns present
3. CDR3/junction fields populated
4. `complete_vdj` matches G3
5. Parity approaches 100%

## Appendix: Notebook Issue

The original audit notebook (`audit.ipynb`) has a bug on cell 4:
```python
# Bug: sequence_id_heavy has duplicates (18 duplicated values)
sequences = df[["sequence_id_heavy", "sequence_heavy"]].dropna()

# Fix: Use 'id' column which is unique (96 unique values)
sequences = df[["id", "sequence_heavy"]].dropna()
sequences = sequences.rename(columns={"id": "sequence_id", "sequence_heavy": "sequence"})
```
