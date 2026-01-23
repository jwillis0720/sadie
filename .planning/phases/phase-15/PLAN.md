# Phase 15: J Gene Matching & CDR3 Annotation Fix — Plan

**Created:** 2026-01-22  
**Status:** Ready for execution

---

## Overview

Fix the `AuxFileBuilder` class in `src/sadie/germlines/builders/aux.py` to generate J gene entries in the correct IgBLAST auxiliary file format. This will enable J gene matching and CDR3/junction annotation.

**Bug Location:** `src/sadie/germlines/builders/aux.py` → `_create_aux_entry()` method

**Root Cause:** J gene entries are generated with 3 columns (FWR4 bounds) instead of 5 columns (IgBLAST aux format).

---

## Execution Plan

### Plan 1: Fix J Gene Aux Format (Core Fix)

**Goal:** Generate J gene entries with correct 5-column IgBLAST format

**Tasks:**

#### Task 1.1: Create J Gene Reference Data Module

Create a new module with reference data for known J gene alleles.

**File:** `src/sadie/germlines/builders/j_gene_data.py`

**Content:**
```python
"""
J Gene Reference Data
=====================

Contains CDR3 end positions and reading frames for J gene alleles.
Data sourced from IMGT reference and G3 aux files.
"""

# Human J gene reference data
# Format: {allele: (reading_frame, cdr3_end, is_functional)}
HUMAN_J_GENE_DATA = {
    # Heavy chain (JH)
    "IGHJ1*01": (0, 17, 1),
    "IGHJ2*01": (1, 18, 1),
    "IGHJ3*01": (1, 15, 1),
    "IGHJ3*02": (1, 15, 1),
    "IGHJ4*01": (2, 13, 1),
    "IGHJ4*02": (2, 13, 1),
    "IGHJ4*03": (2, 13, 1),
    "IGHJ5*01": (2, 16, 1),
    "IGHJ5*02": (2, 16, 1),
    "IGHJ5*03": (2, 16, 1),
    "IGHJ5*04": (2, 16, 1),
    "IGHJ6*01": (2, 28, 1),
    "IGHJ6*02": (2, 28, 0),
    "IGHJ6*03": (2, 28, 0),
    "IGHJ6*04": (2, 28, 1),
    # Kappa chain (JK)
    "IGKJ1*01": (1, 6, 1),
    "IGKJ2*01": (2, 7, 1),
    "IGKJ2*02": (1, 6, 1),
    "IGKJ2*03": (2, 7, 1),
    "IGKJ2*04": (2, 7, 1),
    "IGKJ3*01": (1, 6, 1),
    "IGKJ4*01": (1, 6, 1),
    "IGKJ4*02": (1, 6, 1),
    "IGKJ4*03": (1, 6, 1),
    "IGKJ5*01": (1, 6, 1),
    # Lambda chain (JL)
    "IGLJ1*01": (1, 6, 1),
    "IGLJ2*01": (1, 6, 1),
    "IGLJ2A*01": (1, 6, 1),
    "IGLJ3*01": (1, 6, 1),
    "IGLJ3*02": (1, 6, 1),
    "IGLJ4*01": (1, 6, 1),
    "IGLJ5*01": (1, 6, 1),
    "IGLJ5*02": (1, 6, 1),
    "IGLJ6*01": (1, 6, 1),
    "IGLJ7*01": (1, 6, 1),
    "IGLJ7*02": (1, 6, 1),
}

# Chain type mapping
CHAIN_TYPE_MAP = {
    "H": "JH",
    "K": "JK",
    "L": "JL",
}

def get_j_gene_data(allele_name: str, chain: str) -> tuple:
    """
    Get J gene reference data for an allele.
    
    Parameters
    ----------
    allele_name : str
        Full allele name (e.g., "IGHJ1*01")
    chain : str
        Chain type (H, K, or L)
    
    Returns
    -------
    tuple
        (reading_frame, chain_type, cdr3_end, is_functional)
        Returns default values if allele not found.
    """
    chain_type = CHAIN_TYPE_MAP.get(chain, f"J{chain}")
    
    # Check known reference data
    if allele_name in HUMAN_J_GENE_DATA:
        rf, cdr3_end, is_func = HUMAN_J_GENE_DATA[allele_name]
        return (rf, chain_type, cdr3_end, is_func)
    
    # Default fallback values based on chain type
    # These defaults are based on most common values
    defaults = {
        "H": (1, "JH", 15, 1),  # Most IGHJ are RF1, CDR3 ~15
        "K": (1, "JK", 6, 1),   # Most IGKJ are RF1, CDR3 6
        "L": (1, "JL", 6, 1),   # Most IGLJ are RF1, CDR3 6
    }
    
    return defaults.get(chain, (1, chain_type, 10, 1))
```

**Verification:** File compiles without errors.

---

#### Task 1.2: Modify AuxFileBuilder Class

Update the `AuxFileBuilder` class to:
1. Only generate J gene entries (remove V genes)
2. Use 5-column format for J genes

**File:** `src/sadie/germlines/builders/aux.py`

**Changes:**

1. **Update imports** (top of file):
```python
from sadie.germlines.builders.j_gene_data import get_j_gene_data
```

2. **Remove V genes from processing** (modify `build_for_species` method):
```python
def build_for_species(
    self,
    species: str,
    source_dir: Path,
    output_file: Path
) -> None:
    """
    Build auxiliary file for species.
    
    IgBLAST auxiliary files only contain J gene data.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Building auxiliary file for {species}")

    aux_lines = []

    # Process J segments only (IgBLAST aux files only need J genes)
    for chain in CHAINS:
        lines = self._process_segment(
            species,
            chain,
            "J",  # Only J segments
            source_dir
        )
        aux_lines.extend(lines)

    # Write auxiliary file
    if aux_lines:
        output_file.write_text("\n".join(aux_lines) + "\n")
        logger.info(
            f"Wrote {len(aux_lines)} entries to {output_file}"
        )
    else:
        logger.warning(f"No auxiliary entries generated for {species}")
```

3. **Rewrite J gene entry creation** (replace `_create_aux_entry` method):
```python
def _create_aux_entry(
    self,
    record,
    chain: str,
    segment: str
) -> Optional[str]:
    """
    Create auxiliary file entry for a J gene sequence.

    IgBLAST aux format for J genes (5 columns, tab-separated):
    <gene_name>\t<reading_frame>\t<chain_type>\t<cdr3_end>\t<is_functional>

    Parameters
    ----------
    record : SeqRecord
        Sequence record
    chain : str
        Chain type (H, K, L)
    segment : str
        Segment type (must be "J")

    Returns
    -------
    str or None
        Auxiliary file line for J gene, None for other segments
    """
    if segment != "J":
        # Only J genes go in aux file
        return None

    gene_name = record.id
    
    # Get reference data for this J gene
    reading_frame, chain_type, cdr3_end, is_functional = get_j_gene_data(
        gene_name, chain
    )

    return f"{gene_name}\t{reading_frame}\t{chain_type}\t{cdr3_end}\t{is_functional}"
```

4. **Remove unused methods** (optional cleanup):
   - `_parse_imgt_regions()` - No longer needed for J genes
   - `_build_position_map()` - No longer needed for J genes
   
   (Keep these if V gene aux generation might be needed in future, but mark as unused)

5. **Update module constants** (at top of file):
```python
# Only J segments are needed for aux files
SEGMENTS = ["J"]
```

**Verification:** 
- File compiles without errors
- Import of j_gene_data module works

---

#### Task 1.3: Rebuild Human Aux File

Rebuild the auxiliary file for human species.

**Commands:**
```bash
# In audit/ sandbox, test the rebuild
cd /Users/tmsincomb/sadie
python -c "
from pathlib import Path
from sadie.germlines.builders.aux import AuxFileBuilder

builder = AuxFileBuilder()
builder.build_for_species(
    'human',
    source_dir=Path('src/sadie/germlines/normalized/human/gapped'),
    output_file=Path('audit/test_human_gl.aux')
)

# Verify format
with open('audit/test_human_gl.aux') as f:
    lines = f.readlines()
    print(f'Total entries: {len(lines)}')
    print('First 5 entries:')
    for line in lines[:5]:
        print(line.strip())
        cols = line.strip().split('\t')
        print(f'  Columns: {len(cols)}')
"
```

**Expected Output:**
- ~35 entries (J genes only)
- Each line has 5 tab-separated columns
- Format: `IGHJ1*01	0	JH	17	1`

**Verification:**
- Count columns per line = 5
- All entries are J genes (start with IG*J)
- No V gene entries present

---

### Plan 2: Integration & Validation

**Goal:** Validate the fix enables J gene matching and CDR3 annotation

**Tasks:**

#### Task 2.1: Deploy Fixed Aux File

Copy the validated aux file to the production location.

**Command:**
```bash
# After validation, copy to production location
cp audit/test_human_gl.aux src/sadie/germlines/igblast/aux_db/human_gl.aux
```

**Verification:** File exists and has correct format.

---

#### Task 2.2: Run IgBLAST Validation Test

Create and run a validation script in the audit sandbox.

**File:** `audit/validate_j_gene_fix.py`

**Content:**
```python
"""
Validate J Gene Fix
===================

Tests that J gene matching and CDR3 annotation work after aux file fix.
"""

import pandas as pd
from pathlib import Path

# Use germlines backend
from sadie.airr import Airr

def main():
    # Test sequence (known productive heavy chain)
    test_fasta = Path("audit/test_sequences.fasta")
    
    # Create test FASTA if not exists
    if not test_fasta.exists():
        test_fasta.write_text(
            ">test_seq1\n"
            "CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTC"
            "TCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGGCATGCACTGGGTCCGCCAGGCT"
            "CCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTAC"
            "GCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTAT"
            "CTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCGAGAGATCGA"
            "CGGTTTGCTTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG\n"
        )
    
    # Run annotation with germlines backend
    airr = Airr(
        "test_seq1",
        species="human",
        backend="germlines"  # Use germlines module
    )
    
    result = airr.run_single(str(test_fasta.read_text().split('\n')[1]))
    
    # Check critical fields
    print("=== J Gene Fix Validation ===")
    print(f"v_call: {result.get('v_call', 'N/A')}")
    print(f"d_call: {result.get('d_call', 'N/A')}")
    print(f"j_call: {result.get('j_call', 'N/A')}")  # Should NOT be NaN
    print(f"junction: {result.get('junction', 'N/A')}")  # Should NOT be NaN
    print(f"junction_aa: {result.get('junction_aa', 'N/A')}")
    print(f"cdr3: {result.get('cdr3', 'N/A')}")  # Should NOT be NaN
    print(f"cdr3_aa: {result.get('cdr3_aa', 'N/A')}")
    print(f"fwr4: {result.get('fwr4', 'N/A')}")  # Should NOT be NaN
    print(f"complete_vdj: {result.get('complete_vdj', 'N/A')}")  # Should be True
    
    # Validation checks
    success = True
    if pd.isna(result.get('j_call')) or result.get('j_call') is None:
        print("\n❌ FAIL: j_call is NaN/None")
        success = False
    else:
        print("\n✓ PASS: j_call is populated")
    
    if pd.isna(result.get('junction')) or result.get('junction') is None:
        print("❌ FAIL: junction is NaN/None")
        success = False
    else:
        print("✓ PASS: junction is populated")
    
    if pd.isna(result.get('cdr3')) or result.get('cdr3') is None:
        print("❌ FAIL: cdr3 is NaN/None")
        success = False
    else:
        print("✓ PASS: cdr3 is populated")
    
    print(f"\n=== Overall: {'SUCCESS' if success else 'FAILURE'} ===")
    return success

if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)
```

**Command:**
```bash
cd /Users/tmsincomb/sadie
python audit/validate_j_gene_fix.py
```

**Expected Output:**
- j_call: IGHJ4*02 (or similar)
- junction: populated sequence
- cdr3: populated sequence
- fwr4: populated sequence
- complete_vdj: True

**Verification:** All critical fields are populated (not NaN/None).

---

#### Task 2.3: Run Full Audit Comparison

Compare germlines vs G3 backend on the audit dataset.

**Command:**
```bash
cd /Users/tmsincomb/sadie
python audit/audit.py  # Or run audit notebook
```

**Verification:**
- J gene match rate approaches 100%
- CDR3/junction annotation rate approaches 100%
- complete_vdj rate approaches 100%

---

## Success Criteria

| Criterion | Target | Validation Method |
|-----------|--------|-------------------|
| J gene aux format | 5 columns | Inspect aux file |
| J gene entries only | No V genes | Count line types |
| j_call populated | >95% | Audit comparison |
| junction populated | >95% | Audit comparison |
| cdr3 populated | >95% | Audit comparison |
| fwr4 populated | >95% | Audit comparison |
| complete_vdj = True | >95% | Audit comparison |
| G3 parity | >95% match | Audit comparison |

---

## Risk Mitigation

1. **Fallback Values:** If allele not in reference data, use sensible defaults
2. **Species Support:** Currently focused on human; extend pattern for other species
3. **Backup:** Original aux file can be restored from git

---

## Execution Order

1. **Plan 1, Task 1.1** — Create j_gene_data.py module
2. **Plan 1, Task 1.2** — Modify aux.py builder
3. **Plan 1, Task 1.3** — Rebuild aux file to sandbox
4. **Plan 2, Task 2.1** — Deploy to production location
5. **Plan 2, Task 2.2** — Run validation test
6. **Plan 2, Task 2.3** — Run full audit comparison

---

## Files Modified

| File | Change Type | Description |
|------|-------------|-------------|
| `src/sadie/germlines/builders/j_gene_data.py` | **NEW** | J gene reference data |
| `src/sadie/germlines/builders/aux.py` | MODIFY | Fix J gene aux format |
| `src/sadie/germlines/igblast/aux_db/human_gl.aux` | REGENERATE | Rebuilt aux file |
| `audit/validate_j_gene_fix.py` | **NEW** | Validation script |

---

## Dependencies

- No external dependencies
- No pip installs required
- Uses existing biopython (SeqIO)

---

*Plan ready for execution.*
