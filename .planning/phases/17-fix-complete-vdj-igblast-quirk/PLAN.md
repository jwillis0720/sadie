# Phase 17 Plan: Fix complete_vdj IgBLAST Quirk

## Goal
Ensure `complete_vdj=True` for sequences with valid VDJ alignments by post-processing IgBLAST results using the AIRR standard definition.

## Approach
Post-processing recalculation in `Airr.run_fasta()` after IgBLAST execution, using:
- `v_germline_start == 1` (V alignment complete)
- `j_germline_end == expected_j_length` (J alignment complete)

## Tasks

### Task 17-1: Add J Gene Sequence Lengths to j_gene_data.py

**File:** `src/sadie/germlines/builders/j_gene_data.py`

**Action:** Add `J_GENE_LENGTHS` dictionary with expected nucleotide lengths for each J allele.

```python
# Add after HUMAN_J_GENE_DATA

# J gene sequence lengths (nucleotides) for complete_vdj calculation
J_GENE_LENGTHS = {
    # Heavy chain
    "IGHJ1*01": 52, "IGHJ2*01": 55, "IGHJ3*01": 50, "IGHJ3*02": 50,
    "IGHJ4*01": 48, "IGHJ4*02": 48, "IGHJ4*03": 48,
    "IGHJ5*01": 51, "IGHJ5*02": 51, "IGHJ5*03": 51, "IGHJ5*04": 51,
    "IGHJ6*01": 62, "IGHJ6*02": 62, "IGHJ6*03": 62, "IGHJ6*04": 62,
    # Kappa chain
    "IGKJ1*01": 39, "IGKJ2*01": 36, "IGKJ2*02": 39, "IGKJ2*03": 36, "IGKJ2*04": 36,
    "IGKJ3*01": 36, "IGKJ4*01": 38, "IGKJ4*02": 38, "IGKJ4*03": 38, "IGKJ5*01": 39,
    # Lambda chain
    "IGLJ1*01": 37, "IGLJ2*01": 38, "IGLJ2A*01": 38, "IGLJ3*01": 36, "IGLJ3*02": 36,
    "IGLJ4*01": 35, "IGLJ5*01": 36, "IGLJ5*02": 36, "IGLJ6*01": 35,
    "IGLJ7*01": 38, "IGLJ7*02": 38,
}


def get_j_gene_length(allele_name: str) -> int | None:
    """Get expected J gene length for complete_vdj calculation."""
    if allele_name in J_GENE_LENGTHS:
        return J_GENE_LENGTHS[allele_name]
    # Try gene family (remove allele)
    gene = allele_name.split('*')[0] + '*01'
    return J_GENE_LENGTHS.get(gene)
```

**Verification:**
- Check all J alleles in `human_J.fasta` are covered
- Verify lengths match actual sequences

---

### Task 17-2: Add _recalculate_complete_vdj Method to Airr Class

**File:** `src/sadie/airr/airr.py`

**Action:** Add private method to recalculate complete_vdj from position data.

```python
def _recalculate_complete_vdj(self, result: AirrTable) -> AirrTable:
    """
    Recalculate complete_vdj based on AIRR standard definition.

    IgBLAST's complete_vdj calculation can vary based on allele selection.
    This post-processing ensures consistent results based purely on
    alignment positions per AIRR Standards v1.6.

    AIRR Definition: True if alignment spans entire V(D)J region:
    - v_germline_start == 1 (starts at V gene beginning)
    - j_germline_end == expected_j_length (extends to J gene end)
    """
    from sadie.germlines.builders.j_gene_data import get_j_gene_length

    def calc_complete(row):
        # Handle missing position data
        v_start = row.get('v_germline_start')
        j_end = row.get('j_germline_end')
        j_call = row.get('j_call')

        if pd.isna(v_start) or pd.isna(j_end) or pd.isna(j_call) or not j_call:
            return None

        # Get first allele if multiple
        j_allele = str(j_call).split(',')[0].strip()
        expected_j_len = get_j_gene_length(j_allele)

        if expected_j_len is None:
            return None

        # AIRR standard: spans entire V(D)J region
        v_complete = int(v_start) == 1
        j_complete = int(j_end) == expected_j_len

        return v_complete and j_complete

    result['complete_vdj'] = result.apply(calc_complete, axis=1)
    return result
```

**Location:** Add near other private methods (around line 645)

---

### Task 17-3: Integrate Recalculation in run_fasta

**File:** `src/sadie/airr/airr.py`

**Action:** Call `_recalculate_complete_vdj` before returning from `run_fasta`.

**Before (around line 642-644):**
```python
        if self.correct_indel:
            result.correct_indel()

        return result
```

**After:**
```python
        if self.correct_indel:
            result.correct_indel()

        # Fix IgBLAST complete_vdj quirk (see audit/igblast-quirk.md)
        result = self._recalculate_complete_vdj(result)

        return result
```

---

### Task 17-4: Integrate Recalculation in _run_scfv

**File:** `src/sadie/airr/airr.py`

**Action:** Also apply recalculation to scfv results (linked AirrTable).

**Location:** In `_run_scfv` method, before creating LinkedAirrTable.

Find the return statement for scfv (around line 530-545) and add recalculation.

---

### Task 17-5: Run Audit to Verify Fix

**Command:**
```bash
cd /Users/tmsincomb/sadie && conda run -n sadie-dev python audit/audit.py
```

**Expected Results:**
- complete_vdj differences: 0 (was 22)
- Pure structural parity: 99%+ (was 98.29%)

---

### Task 17-6: Update igblast-quirk.md with Resolution

**File:** `audit/igblast-quirk.md`

**Action:** Add "Resolution" section documenting the fix.

---

## Success Criteria

1. [ ] `complete_vdj` matches G3 for all 78 test sequences
2. [ ] No more False negatives for productive sequences
3. [ ] Pure structural parity >= 99%
4. [ ] Tests pass (if any exist for complete_vdj)

## Files Modified

| File | Change |
|------|--------|
| `src/sadie/germlines/builders/j_gene_data.py` | Add J_GENE_LENGTHS dict |
| `src/sadie/airr/airr.py` | Add _recalculate_complete_vdj, integrate in run_fasta |
| `audit/igblast-quirk.md` | Add resolution section |

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| Missing J alleles in lookup | Fallback to gene family (*01) |
| Performance impact | Single apply() operation, negligible |
| Breaking existing behavior | Only changes complete_vdj, all other fields unchanged |

---
*Plan created: 2026-01-22*
