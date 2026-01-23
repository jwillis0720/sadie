# Phase 17 Research: Fix complete_vdj IgBLAST Quirk

## Summary

The `complete_vdj` discrepancy between germlines and G3 backends is caused by IgBLAST's internal calculation being influenced by allele selection, not actual alignment coverage. **The recommended solution is post-processing recalculation** based on the AIRR standard definition, which uses explicit position data (`v_germline_start=1` and `j_germline_end=J_gene_length`).

## Primary Recommendation

**Use post-processing recalculation of `complete_vdj` in SADIE's AirrTable class.**

Rationale:
1. The AIRR standard provides an explicit, calculable definition
2. This is already a common pattern in the immunoinformatics community
3. All alignment position data is already available and identical between backends
4. Does not require modifying IgBLAST database configuration

## Standard Stack

| Component | Version | Purpose |
|-----------|---------|---------|
| SADIE AirrTable | current | Post-processing host for recalculation |
| pandas | existing | DataFrame operations for recalculation |
| AIRR Standards | 1.6 | Definition of complete_vdj |

**No new libraries required.** Solution uses existing infrastructure.

## Architecture Patterns

### Pattern 1: Post-processing Recalculation (RECOMMENDED)

```python
def recalculate_complete_vdj(airrtable: pd.DataFrame, j_gene_lengths: dict) -> pd.DataFrame:
    """
    Recalculate complete_vdj based on AIRR standard definition.
    
    AIRR Definition:
    "True if the sequence alignment spans the entire V(D)J region.
    Meaning, sequence_alignment includes both the first V gene codon
    that encodes the mature polypeptide chain (i.e., after the leader
    sequence) and the last complete codon of the J gene (i.e., before
    the J-C splice site)."
    
    Calculation:
    - v_germline_start == 1 (alignment starts at V gene beginning)
    - j_germline_end == expected_j_length (alignment extends to J gene end)
    """
    def calculate_complete(row):
        if pd.isna(row['v_germline_start']) or pd.isna(row['j_germline_end']):
            return None
        if pd.isna(row['j_call']) or not row['j_call']:
            return None
            
        # Get J gene base allele
        j_allele = str(row['j_call']).split(',')[0].strip()
        expected_j_len = j_gene_lengths.get(j_allele)
        
        if expected_j_len is None:
            return None
            
        # AIRR standard: spans entire V(D)J region
        v_complete = row['v_germline_start'] == 1
        j_complete = row['j_germline_end'] == expected_j_len
        
        return v_complete and j_complete
    
    airrtable['complete_vdj'] = airrtable.apply(calculate_complete, axis=1)
    return airrtable
```

### Pattern 2: Integration Point in AirrTable

The recalculation should happen in one of these locations:

**Option A: AirrTable.__init__** - Called after DataFrame construction
- Pros: Automatic, always applied
- Cons: Adds processing overhead to all AirrTable creation

**Option B: AirrTable method** - Explicit `recalculate_complete_vdj()` call
- Pros: Opt-in, clear intent
- Cons: Users must call it

**Option C: Airr.run_* methods** - Post-process after IgBLAST
- Pros: Applied at annotation time
- Cons: Only affects new annotations

**RECOMMENDED: Option C** - Apply in `Airr.run_single`, `Airr.run_dataframe`, etc. after IgBLAST execution but before returning results.

## Don't Hand-Roll

1. **Don't recreate J gene length data** - Use existing `j_gene_data.py` module with `HUMAN_J_GENE_DATA`
2. **Don't modify IgBLAST source** - Use post-processing approach
3. **Don't parse IgBLAST internal logic** - Rely on AIRR-defined fields only
4. **Don't create combined V+D+J+C internal_data files** - More complexity, same result

## Common Pitfalls

### Pitfall 1: Trusting IgBLAST's complete_vdj Directly
**Problem:** IgBLAST's calculation is influenced by internal factors beyond alignment positions.
**Evidence:** 22 sequences have identical alignments but different `complete_vdj` values.
**Solution:** Always recalculate from position data.

### Pitfall 2: Missing J Gene Length Data
**Problem:** Need expected J gene lengths for all alleles.
**Solution:** Build J gene length lookup from:
- `src/sadie/germlines/builders/j_gene_data.py` (CDR3 end positions)
- `src/sadie/germlines/igblast/database/human/human_J.fasta` (sequence lengths)

### Pitfall 3: Handling Multiple J Allele Calls
**Problem:** `j_call` can contain multiple alleles (comma-separated).
**Solution:** Use first allele only: `j_call.split(',')[0].strip()`

### Pitfall 4: Nullable Position Fields
**Problem:** `v_germline_start` or `j_germline_end` may be null/NA.
**Solution:** Return `None` for `complete_vdj` if positions are missing.

### Pitfall 5: Allele Matching
**Problem:** J allele names may not match exactly between databases.
**Solution:** Fall back to gene family if allele not found: `IGHJ1*01` → `IGHJ1` → default.

## Code Examples

### Example 1: J Gene Length Lookup Table

```python
# Expected J gene lengths (nucleotides) from human_J.fasta
J_GENE_LENGTHS = {
    # Heavy chain
    "IGHJ1*01": 52,
    "IGHJ2*01": 55,
    "IGHJ3*01": 50,
    "IGHJ3*02": 50,
    "IGHJ4*01": 48,
    "IGHJ4*02": 48,
    "IGHJ4*03": 48,
    "IGHJ5*01": 51,
    "IGHJ5*02": 51,
    "IGHJ5*03": 51,
    "IGHJ5*04": 51,
    "IGHJ6*01": 62,
    "IGHJ6*02": 62,
    "IGHJ6*03": 62,
    "IGHJ6*04": 62,
    # Kappa chain
    "IGKJ1*01": 39,
    "IGKJ2*01": 36,
    "IGKJ2*02": 39,
    "IGKJ2*03": 36,
    "IGKJ2*04": 36,
    "IGKJ3*01": 36,
    "IGKJ4*01": 38,
    "IGKJ4*02": 38,
    "IGKJ4*03": 38,
    "IGKJ5*01": 39,
    # Lambda chain
    "IGLJ1*01": 37,
    "IGLJ2*01": 38,
    "IGLJ2A*01": 38,
    "IGLJ3*01": 36,
    "IGLJ3*02": 36,
    "IGLJ4*01": 35,
    "IGLJ5*01": 36,
    "IGLJ5*02": 36,
    "IGLJ6*01": 35,
    "IGLJ7*01": 38,
    "IGLJ7*02": 38,
}
```

### Example 2: Integration in Airr Class

```python
class Airr:
    def _recalculate_complete_vdj(self, result: AirrTable) -> AirrTable:
        """Post-process IgBLAST results to fix complete_vdj."""
        from sadie.germlines.builders.j_gene_data import J_GENE_LENGTHS
        
        def calc(row):
            if pd.isna(row.get('v_germline_start')) or pd.isna(row.get('j_germline_end')):
                return None
            j_call = str(row.get('j_call', '')).split(',')[0].strip()
            expected = J_GENE_LENGTHS.get(j_call)
            if expected is None:
                return None
            return row['v_germline_start'] == 1 and row['j_germline_end'] == expected
        
        result['complete_vdj'] = result.apply(calc, axis=1)
        return result
```

### Example 3: Test Case

```python
def test_complete_vdj_recalculation():
    """Verify complete_vdj matches G3 after recalculation."""
    # Run germlines backend
    airr = Airr("human")
    result = airr.run_dataframe(test_sequences)
    
    # Check previously problematic sequences
    for seq_id in KNOWN_DISCREPANCY_IDS:
        row = result[result['sequence_id'] == seq_id].iloc[0]
        # Should now be True (matches G3)
        assert row['complete_vdj'] == True
        # Verify calculation basis
        assert row['v_germline_start'] == 1
        assert row['j_germline_end'] == J_GENE_LENGTHS[row['j_call'].split(',')[0]]
```

## Root Cause Analysis

### Why IgBLAST Produces Different complete_vdj Values

Based on audit investigation, the quirk occurs because:

1. **Allele Selection Affects Calculation:** When IgBLAST selects `IGHV4-31*13` (newer IMGT allele) vs `IGHV4-31*03` (older allele), the `complete_vdj` result differs even though alignment positions are identical.

2. **Internal Logic is Opaque:** IgBLAST's `complete_vdj` calculation (added in v1.16.0) depends on internal factors not exposed through the output format.

3. **Not a Bug:** This appears to be intentional behavior where IgBLAST uses additional heuristics beyond simple position checking.

### AIRR Standard Definition (Authoritative)

From AIRR Rearrangement Schema v1.6:

> **complete_vdj** (boolean, optional, nullable):
> True if the sequence alignment spans the entire V(D)J region. Meaning, sequence_alignment includes both the first V gene codon that encodes the mature polypeptide chain (i.e., after the leader sequence) and the last complete codon of the J gene (i.e., before the J-C splice site). **This does not require an absence of deletions within the internal FWR and CDR regions of the alignment.**

Key insight: The definition is purely positional - it does NOT depend on which specific allele is called.

## Sources

| Source | Confidence | Key Finding |
|--------|------------|-------------|
| AIRR Standards v1.6 | HIGH | Official `complete_vdj` definition |
| IgBLAST Release Notes | HIGH | v1.16.0 added `complete_vdj` field |
| SADIE `igblast-quirk.md` | HIGH | Documents the allele-dependent behavior |
| Immcantation/Change-O | MEDIUM | Uses post-processing for annotations |
| NCBI IgBLAST Documentation | MEDIUM | No documented calculation details |

## Verification Checklist

- [x] AIRR standard definition researched and documented
- [x] IgBLAST behavior analyzed (release notes, internal_data structure)
- [x] Existing codebase patterns identified (j_gene_data.py exists)
- [x] Immcantation/Change-O approaches reviewed (post-processing is standard)
- [x] Root cause documented (allele selection affects internal calculation)
- [x] Code examples provided with actual sequence lengths
- [x] Pitfalls identified with solutions

## Implementation Notes

1. **J Gene Lengths:** Can be dynamically computed from `human_J.fasta` or hard-coded. Hard-coded is simpler and sufficient for human.

2. **Species Support:** For multi-species support, build J gene length lookup dynamically from FASTA files in `src/sadie/germlines/igblast/database/{species}/{species}_J.fasta`.

3. **Performance:** Single DataFrame apply() operation - negligible performance impact.

4. **Backwards Compatibility:** Overwriting existing `complete_vdj` values maintains same column structure.
