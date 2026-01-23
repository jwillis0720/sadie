# IgBLAST `complete_vdj` Quirk

## Summary

IgBLAST produces different `complete_vdj` values when selecting between nearly-identical V gene alleles, even when all alignment coordinates are identical.

## Discovery Context

During the AIRR backend parity audit comparing the germlines module vs G3 backend, 22 sequences showed `complete_vdj=False` in germlines but `complete_vdj=True` in G3.

## Root Cause

| Database | Allele Selected | complete_vdj | 
|----------|-----------------|--------------|
| Germlines (IMGT) | IGHV4-31*13 | **False** |
| G3 (older IMGT) | IGHV4-31*03 | True |

The germlines database contains `IGHV4-31*13` (a newer IMGT allele) which is not present in G3. When IgBLAST selects *13 as the best match, it marks `complete_vdj=False`, but when *13 is unavailable and *03 is selected instead, `complete_vdj=True`.

## Evidence: Identical Alignment Coordinates

Both alleles produce **identical** alignment coordinates:

```
v_sequence_start:  1   (both)
v_sequence_end:    296 (both)
v_germline_start:  1   (both)
v_germline_end:    298 (both)
j_germline_start:  6   (both)
j_germline_end:    62  (both)
```

## Evidence: Identical Allele Properties

| Property | IGHV4-31*03 | IGHV4-31*13 |
|----------|-------------|-------------|
| Sequence length | 299 nt | 299 nt |
| ndm.imgt FR3_end | 291 | 291 |
| All region boundaries | identical | identical |

The only sequence difference is at position 295:
- *03: `...TGTGCGAGAGA`
- *13: `...TGTGCCAGAGA`

## Test Sequence Analysis

Query sequence ends with: `...TGTGCCAGAGGGGGCCGG...`

- *13 is a better match (exact match up to CDR3 junction)
- *03 has a mismatch at position 295 (G vs C)

IgBLAST correctly selects *13 as the better match when available, but incorrectly marks `complete_vdj=False`.

## Conclusion

This appears to be an IgBLAST internal logic quirk where `complete_vdj` determination is affected by factors beyond just alignment coordinates. The flag is **not reliable** for comparing annotations across databases with different allele sets.

## Recommendation

When auditing backend parity:
1. Treat `complete_vdj` as an allele-dependent field
2. Focus on alignment coordinates (v_germline_start/end, j_germline_start/end) for structural comparison
3. The actual annotation quality is identical - only the derived flag differs

## Affected Sequences (22 total)

All 22 sequences with `complete_vdj` discrepancy share these characteristics:
- Productive annotation (productive=True in both backends)
- Complete CDR3, FWR4, and junction annotations
- Identical alignment positions
- Different V allele selection due to database differences

## Resolution (Phase 17)

**Status:** Fixed in commit `2dd31243` (2026-01-22)

### Solution

Post-processing recalculation in `Airr.run_fasta()` and `Airr._run_scfv()` using AIRR standard definition:

```python
def _recalculate_complete_vdj(self, result: AirrTable) -> AirrTable:
    """
    AIRR Definition: True if alignment spans entire V(D)J region:
    - v_germline_start == 1 (starts at V gene beginning)
    - j_germline_end == expected_j_length (extends to J gene end)
    """
    v_complete = v_germline_start == 1
    j_complete = j_germline_end == get_j_gene_length(j_call)
    return v_complete and j_complete
```

### Implementation

1. **J Gene Lengths:** Added `J_GENE_LENGTHS` dictionary to `j_gene_data.py` with expected nucleotide lengths for all human J alleles (sourced from SADIE germline FASTA files)

2. **Recalculation Method:** Added `_recalculate_complete_vdj()` method to `Airr` class that:
   - Checks `v_germline_start == 1` (V alignment complete)
   - Checks `j_germline_end == expected_j_length` (J alignment complete)
   - Returns True only if both conditions are met

3. **Integration:** Called after IgBLAST execution in both:
   - `run_fasta()` for standard annotations
   - `_run_scfv()` for linked H+L chain annotations

### Results

| Metric | Before | After |
|--------|--------|-------|
| complete_vdj differences vs G3 | 22 | 4 |
| Direction of differences | germlines=False, G3=True | germlines=True, G3=False |

The 4 remaining differences now show germlines with the **correct** value (True) per AIRR standard, while G3 shows the incorrect value (False). These are sequences where:
- `v_germline_start == 1` ✓
- `j_germline_end == expected_j_length` ✓
- Should be `complete_vdj=True` per AIRR standard

The fix ensures SADIE produces reliable, standard-compliant `complete_vdj` values regardless of which allele IgBLAST selects.
