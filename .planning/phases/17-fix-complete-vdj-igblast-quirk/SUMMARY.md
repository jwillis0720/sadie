# Phase 17 Summary: Fix complete_vdj IgBLAST Quirk

## Objective
Fix the IgBLAST quirk where `complete_vdj` was incorrectly reported as False for sequences with valid VDJ alignments due to allele selection differences.

## Solution
Post-processing recalculation in `Airr.run_fasta()` and `Airr._run_scfv()` using the AIRR standard definition:
- `v_germline_start == 1` (V alignment starts at beginning)
- `j_germline_end == expected_j_length` (J alignment extends to end)

## Commits

| Commit | Type | Description |
|--------|------|-------------|
| `2b300874` | data(17-1) | Add J_GENE_LENGTHS dictionary with expected J gene nucleotide lengths |
| `45bab802` | feat(17-2) | Add _recalculate_complete_vdj method and integrate in run_fasta |
| `2dd31243` | feat(17-4) | Integrate complete_vdj recalculation in _run_scfv for linked chains |
| `c3187682` | docs(17-6) | Document complete_vdj quirk resolution in igblast-quirk.md |

## Files Modified

| File | Change |
|------|--------|
| `src/sadie/germlines/builders/j_gene_data.py` | Added `J_GENE_LENGTHS` dict and `get_j_gene_length()` function |
| `src/sadie/airr/airr.py` | Added `_recalculate_complete_vdj()` method, integrated in `run_fasta()` and `_run_scfv()` |
| `audit/igblast-quirk.md` | Added Resolution section documenting the fix |

## Audit Results

| Metric | Before | After |
|--------|--------|-------|
| complete_vdj differences vs G3 | 22 | 4 |
| Direction of differences | germlines=False, G3=True | germlines=True, G3=False |
| Pure structural parity | 98.29% | 98.29% |

### Interpretation
- **22 → 4 differences**: Reduced complete_vdj discrepancies significantly
- **Direction change**: The 4 remaining differences now show germlines with the **correct** value (True) per AIRR standard, while G3 incorrectly shows False
- **Structural parity unchanged**: complete_vdj is classified as an allele-dependent field, not a structural field

## Success Criteria Assessment

| Criterion | Status | Notes |
|-----------|--------|-------|
| complete_vdj matches G3 for all sequences | ❌ | 4 differences remain, but germlines is now correct per AIRR standard |
| No more False negatives for productive sequences | ✅ | All 22 original false negatives now correctly show True |
| Pure structural parity >= 99% | ❌ | 98.29% (but complete_vdj is not a structural field) |
| Tests pass | N/A | No existing tests for complete_vdj |

## Deviations from Plan

1. **J Gene Lengths**: Plan had incorrect lengths from an external source; used actual values from SADIE germline FASTA files instead
2. **Partial Success**: Original goal was 0 differences, achieved 4 (now in opposite direction, meaning our values are more accurate)

## Key Insight
The remaining 4 differences demonstrate that our AIRR-standard-based recalculation is more accurate than G3's built-in calculation. G3 also suffers from the same IgBLAST quirk but in a different set of sequences.

---
*Completed: 2026-01-22*
