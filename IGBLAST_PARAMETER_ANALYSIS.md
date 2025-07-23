# IgBLAST Parameter Analysis and CDR3 Detection

## Summary

After extensive testing, we found that **parameter differences between SADIE and the NCBI web interface do NOT explain the CDR3 detection failure** for sequences like Seq5. The issue lies deeper in the CDR3 detection algorithm itself, not in the alignment parameters.

## Web IgBLAST Default Parameters

Based on analysis of the NCBI IgBLAST web interface (https://www.ncbi.nlm.nih.gov/igblast/), the default parameters are:

| Parameter | Web Default | SADIE Default |
|-----------|-------------|---------------|
| **Organism** | Human | Human |
| **V gene mismatch penalty** | -1 | -1 |
| **D gene mismatch penalty** | -2 | -1 |
| **J gene mismatch penalty** | -2 | -2 |
| **Min D gene matches** | 5 | 5 |
| **Number of V alignments** | 1 | 3 |
| **Number of D alignments** | 1 | 3 |
| **Number of J alignments** | 1 | 3 |
| **Extend align 5' end** | False | True |
| **Extend align 3' end** | False | True |
| **Allow V(D)J overlap** | False | False |
| **Domain system** | IMGT | IMGT |
| **Output format** | 7 (tabular) | 19 (AIRR) |

## Implementation Changes

### 1. Added Configurable Parameters to Airr Class

```python
def __init__(
    self,
    reference_name: str,
    # ... existing parameters ...
    # New IgBLAST parameters
    num_alignments_v: Optional[int] = None,
    num_alignments_d: Optional[int] = None,
    num_alignments_j: Optional[int] = None,
    extend_align5end: Optional[bool] = None,
    extend_align3end: Optional[bool] = None,
    min_d_match: Optional[int] = None,
    word_size: Optional[int] = None,
    gap_open: Optional[int] = None,
    gap_extend: Optional[int] = None,
):
```

### 2. Usage Example

```python
# Using web-like defaults
airr_web = Airr(
    "human",
    d_gene_penalty=-2,
    num_alignments_v=1,
    num_alignments_d=1,
    num_alignments_j=1,
    extend_align5end=False,
    extend_align3end=False,
)

# Using more permissive settings
airr_permissive = Airr(
    "human",
    v_gene_penalty=-1,
    d_gene_penalty=-4,
    j_gene_penalty=-3,
    allow_vdj_overlap=True,
    word_size=4,
    gap_open=4,
    gap_extend=1,
)
```

## Test Results

Testing on `cdr3_known_bugs.fasta` with 5 sequences:

| Configuration | CDR3 Detection Rate | Failed Sequences |
|---------------|-------------------|------------------|
| SADIE Defaults | 4/5 (80%) | Seq5 |
| Web Defaults | 4/5 (80%) | Seq5 |
| Permissive Settings | 4/5 (80%) | Seq5 |

## Key Findings

1. **Parameter changes do not improve CDR3 detection** for problematic sequences
2. **The issue is algorithmic**, not parameter-based
3. **Seq5 fails consistently** regardless of parameter settings

### Why Seq5 Fails

The junction structure of Seq5:
```
V gene ends:    TGTCAGCAGTATAATAACTGG
N-nucleotides:  TCGGGT (6 nt insertion)
J gene starts:  GGGACCAAGGTGGAAATCAAAC (starting at J position 17)
```

IgBLAST correctly identifies:
- V gene: IGKV3-15*02
- J gene: IGKJ1*01
- Junction nucleotides

But fails to determine CDR3 boundaries because:
1. The J gene alignment starts at position 17 instead of position 1
2. The auxiliary file CDR3 end position (6) cannot be properly mapped
3. No fallback mechanism exists for partial J alignments

## Recommendations

### 1. Short-term Workarounds

Users can now adjust parameters if needed:
```python
# Try different penalty settings
airr = Airr("human", j_gene_penalty=-3, v_gene_penalty=-1)
```

### 2. Long-term Solutions

The real solution requires modifying the CDR3 detection algorithm to:
- Handle partial J gene alignments
- Implement fallback CDR3 boundary detection
- Use heuristic methods when standard detection fails

### 3. Best Practices

- Always check for CDR3 detection failures in your analysis
- Consider using multiple annotation tools for problematic sequences
- Flag sequences without CDR3 for manual review

## Conclusion

While SADIE now supports customizable IgBLAST parameters matching the web interface, the CDR3 detection issue for sequences like Seq5 remains. This is a fundamental limitation of the current IgBLAST CDR3 detection algorithm when dealing with sequences that have:

1. Significant N-nucleotide insertions
2. Partial J gene alignments
3. Non-standard junction structures

The web interface likely uses additional post-processing or a different CDR3 detection algorithm that can handle these edge cases better than the command-line version.
