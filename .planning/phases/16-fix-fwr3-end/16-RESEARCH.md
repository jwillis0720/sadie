# Phase 16 Research: Fix ndm.imgt FWR3 End Position

## Summary

Phase 16 fixes the IgBLAST internal data file (`ndm.imgt`) generation to use the correct FWR3 end position (IMGT position 312) instead of the full V gene sequence length. This is a data pipeline fix, not a library integration task.

**Primary Recommendation:** Modify `build_internal_data.py` to output `regions['FR3'][1]` (the ungapped FWR3 end position) as column 11 instead of `seq_len` (full V gene length). Regenerate `human.ndm.imgt` and validate with audit.

---

## Standard Stack

| Component | Library/Tool | Version | Purpose |
|-----------|-------------|---------|---------|
| Sequence Parsing | BioPython | >=1.80 | Parse FASTA sequences |
| Data Validation | pandas | >=1.5.0 | Audit comparisons |
| IgBLAST | NCBI IgBLAST | >=1.20.0 | V gene annotation |
| Python | Python | >=3.10 | Runtime |

**Note:** This phase requires no new libraries. All tools are already in use in the codebase.

---

## Architecture Patterns

### NDM File Format (IgBLAST Internal Data)

The `.ndm.imgt` file is a 13-column tab-separated file specifying V gene region boundaries:

| Column | Field | Description | Example |
|--------|-------|-------------|---------|
| 1 | gene_name | V gene allele name | IGHV1-69*01 |
| 2 | fwr1_start | FWR1 start position (1-based) | 1 |
| 3 | fwr1_end | FWR1 end position | 75 |
| 4 | cdr1_start | CDR1 start position | 76 |
| 5 | cdr1_end | CDR1 end position | 99 |
| 6 | fwr2_start | FWR2 start position | 100 |
| 7 | fwr2_end | FWR2 end position | 150 |
| 8 | cdr2_start | CDR2 start position | 151 |
| 9 | cdr2_end | CDR2 end position | 174 |
| 10 | fwr3_start | FWR3 start position | 175 |
| **11** | **fwr3_end** | **FWR3 end position (NOT seq_len)** | **288** |
| 12 | chain_type | Chain type (VH, VK, VL) | VH |
| 13 | frame_flag | Coding frame start (0-based) | 0 |

**Critical Insight:** Column 11 MUST be the FWR3 end position, NOT the full V gene sequence length. Using `seq_len` causes FWR3 to extend into CDR3 territory.

### IMGT Numbering System

IMGT unique numbering for V-DOMAIN defines conserved positions across species:

| Region | IMGT Positions | Nucleotide Positions (Gapped) |
|--------|----------------|-------------------------------|
| FR1 | 1-26 | 1-78 |
| CDR1 | 27-38 | 79-114 |
| FR2 | 39-55 | 115-165 |
| CDR2 | 56-65 | 166-195 |
| FR3 | 66-104 | 196-312 |
| CDR3 | 105+ | After 312 |

**Key Position:** IMGT position 312 marks the END of FR3 in gapped coordinates. The ungapped position varies by sequence (typically 285-297 for human V genes).

### Gapped-to-Ungapped Conversion

IMGT gapped sequences use `.` characters to maintain alignment. To convert:

```python
def calculate_ungapped_positions(gapped_seq: str) -> Dict[str, tuple]:
    """Convert IMGT gapped positions to ungapped positions."""
    ungapped_pos = 0
    pos_map = {}
    
    for gapped_pos, char in enumerate(gapped_seq, 1):
        if char not in ".-":
            ungapped_pos += 1
            pos_map[gapped_pos] = ungapped_pos
    
    # Find ungapped position for IMGT region boundaries
    regions = {}
    for region_name, (gapped_start, gapped_end) in IMGT_V_REGIONS.items():
        start_pos = None
        end_pos = None
        
        for g_pos in range(gapped_start, min(gapped_end + 1, len(gapped_seq) + 1)):
            if g_pos in pos_map:
                if start_pos is None:
                    start_pos = pos_map[g_pos]
                end_pos = pos_map[g_pos]
        
        if start_pos and end_pos:
            regions[region_name] = (start_pos, end_pos)
    
    return regions, ungapped_pos  # ungapped_pos is seq_len
```

---

## Don't Hand-Roll

| Problem | Existing Solution | Location |
|---------|------------------|----------|
| IMGT region definitions | `IMGT_V_REGIONS` constant | `build_internal_data.py` |
| Gapped→ungapped conversion | `calculate_ungapped_positions()` | `build_internal_data.py` |
| NDM entry generation | `generate_ndm_entry()` | `build_internal_data.py` |
| IMGT FASTA parsing | BioPython SeqIO | Standard library |
| Audit comparison | `audit.py` script | `audit/audit.py` |

The codebase already has the infrastructure. The fix is a one-line change in `generate_ndm_entry()`.

---

## Common Pitfalls

### 1. Confusing seq_len with FWR3 end
**Problem:** Using full V gene sequence length for column 11 causes FWR3 to extend past position 312.
**Impact:** FWR3 ~8nt too long, CDR3 ~8nt too short, all junction calculations offset.
**Solution:** Use `regions['FR3'][1]` not `seq_len`.

### 2. Off-by-one in position mapping
**Problem:** IMGT positions are 1-based, but code may assume 0-based.
**Solution:** Ensure `enumerate(gapped_seq, 1)` to start at position 1.

### 3. Incomplete V gene sequences
**Problem:** Some V genes may be truncated before position 312.
**Impact:** FWR3 end will be < 312 in gapped coordinates.
**Handling:** The code handles this by finding the last non-gap position <= 312.

### 4. Different CDR2 lengths
**Problem:** CDR1/CDR2 can have insertions at specific IMGT positions.
**Impact:** Ungapped FWR3 end varies (285-297 for human VH).
**Note:** This is expected behavior, not a bug.

### 5. Missing regeneration after code fix
**Problem:** Fixing the code but not regenerating the data file.
**Solution:** Always run `python build_internal_data.py human` after code changes.

---

## Code Examples

### The Bug (Before)

```python
# WRONG: Using seq_len for column 11
entry = (
    f"{gene_name}\t"
    f"{regions['FR1'][0]}\t{regions['FR1'][1]}\t"
    f"{regions['CDR1'][0]}\t{regions['CDR1'][1]}\t"
    f"{regions['FR2'][0]}\t{regions['FR2'][1]}\t"
    f"{regions['CDR2'][0]}\t{regions['CDR2'][1]}\t"
    f"{regions['FR3'][0]}\t{seq_len}\t"  # BUG: seq_len not FR3 end!
    f"{chain_type}\t0"
)
```

### The Fix (After)

```python
# CORRECT: Using regions['FR3'][1] for column 11
entry = (
    f"{gene_name}\t"
    f"{regions['FR1'][0]}\t{regions['FR1'][1]}\t"
    f"{regions['CDR1'][0]}\t{regions['CDR1'][1]}\t"
    f"{regions['FR2'][0]}\t{regions['FR2'][1]}\t"
    f"{regions['CDR2'][0]}\t{regions['CDR2'][1]}\t"
    f"{regions['FR3'][0]}\t{regions['FR3'][1]}\t"  # FIX: Use FR3 end
    f"{chain_type}\t0"
)
```

### Validation Script

```python
# Compare ndm.imgt values with reference
import pandas as pd

def validate_ndm_values(generated_path, reference_path, gene="IGHV1-69*01"):
    """Validate that generated ndm matches reference for key gene."""
    gen_df = pd.read_csv(generated_path, sep='\t', header=None)
    ref_df = pd.read_csv(reference_path, sep='\t', header=None)
    
    gen_row = gen_df[gen_df[0] == gene].iloc[0]
    ref_row = ref_df[ref_df[0] == gene].iloc[0]
    
    # Column 11 (index 10) is FWR3 end
    gen_fwr3_end = gen_row[10]
    ref_fwr3_end = ref_row[10]
    
    print(f"Generated FWR3 end: {gen_fwr3_end}")
    print(f"Reference FWR3 end: {ref_fwr3_end}")
    assert gen_fwr3_end == ref_fwr3_end, f"Mismatch: {gen_fwr3_end} != {ref_fwr3_end}"
```

---

## Sources

| Source | Confidence | Notes |
|--------|------------|-------|
| NCBI IgBLAST setup docs | HIGH | Official ndm file format documentation |
| IMGT Scientific Chart | HIGH | Authoritative IMGT numbering reference |
| `tests/data/fixtures/reference/igblast_internal/human.ndm.imgt` | HIGH | Reference values from NCBI IgBLAST distribution |
| Existing `build_internal_data.py` | HIGH | Current implementation to modify |
| Audit results (STATE.md) | HIGH | Empirical parity measurements |

---

## Implementation Notes

### Files to Modify
1. `src/sadie/germlines/scripts/build_internal_data.py` — Fix column 11 generation
2. `src/sadie/germlines/igblast/Ig/internal_data/human/human.ndm.imgt` — Regenerate

### Verification Steps
1. Run `python build_internal_data.py human` to regenerate
2. Compare IGHV1-69*01 column 11: should be 288, not 296
3. Run `audit/audit.py` to measure parity improvement
4. Target: Parity increase from ~77.6% toward 95%+

### Expected Impact
- FWR3 will be ~8 nucleotides shorter (correct length)
- CDR3 will be ~8 nucleotides longer (correct length)
- Junction boundaries will align with G3 backend
- Parity improvement: 77.6% → ~86-90%+ expected

---

*Research completed: 2026-01-22*
*Phase status: Implementation-ready*
