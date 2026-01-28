# Summary 32-01: Create Combined VDJC Files for internal_data

**Status:** ✅ Completed  
**Commit:** 0fb7fdbc  
**Date:** 2026-01-27

## What Was Done

### Task 1: Modified build_internal_data.py
- Replaced symlink creation with FASTA concatenation (V+D+J+C segments)
- Added deduplication logic to handle duplicate sequences across files
- Build BLAST database from combined file using makeblastdb
- Combined file named `{species}_V.fasta` as required by IgBLAST

### Task 2: Updated Reference Builder
- Modified `_make_internal_annotaion_file()` to include all segments in BLAST db
- NDM file generation unchanged (V genes only for framework/CDR annotations)
- BLAST database now contains V+D+J+C sequences for complete_vdj calculation

### Task 3: Added Unit Tests
- Created `tests/unit/germlines/test_build_internal_data.py`
- Tests verify no symlinks exist in internal_data
- Tests verify combined FASTA contains all segment types
- Tests verify GermlineData paths point to correct locations
- Tests verify deduplication works correctly

### Task 4: Updated GermlineData Paths
- V/D/J/C dirs now point to `database/{species}/` for IgBLAST search
- `igdata` points to `Ig/` containing `internal_data/` with combined VDJC
- This separation allows IgBLAST to:
  - Use combined VDJC file for internal annotation (complete_vdj)
  - Use separate segment BLAST dbs for alignment (-germline_db_V, etc.)

## Files Modified

- `src/sadie/germlines/scripts/build_internal_data.py` - Combined FASTA creation
- `src/sadie/reference/reference.py` - Reference builder VDJC support
- `src/sadie/airr/igblast/germline.py` - GermlineData path updates
- `tests/unit/germlines/test_build_internal_data.py` - New test file
- `src/sadie/germlines/igblast/Ig/internal_data/*/` - Rebuilt directories

## Verification

```bash
# IgBLAST test shows complete_vdj=True for full-length sequence
python -c "
from sadie.airr import Airr
airr = Airr(reference_name='human')
result = airr.run_single('test', 'CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCTATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCGAGAGATCGACTCTCCATCACGATCCGGCCCCGGTACTACGGCCTGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG')
print(f'complete_vdj: {result[\"complete_vdj\"].values[0]}')
"
# Output: complete_vdj: True

# All 9 new tests pass
pytest tests/unit/germlines/test_build_internal_data.py -v
# Result: 9 passed

# No symlinks in internal_data
find src/sadie/germlines/igblast/Ig/internal_data -type l | wc -l
# Output: 0

# Combined FASTA has all segments
grep -c "^>" src/sadie/germlines/igblast/Ig/internal_data/human/human_V.fasta
# Output: 881 (V+D+J+C combined)
```

## Key Changes Summary

1. **internal_data/{species}/** now contains:
   - `{species}_V.fasta` - Combined VDJC FASTA (all segments)
   - `{species}_V.n*` - BLAST database files from combined FASTA
   - `{species}.ndm.imgt` - NDM annotations (V genes only)

2. **database/{species}/** unchanged:
   - Separate `{species}_V.fasta`, `{species}_D.fasta`, etc.
   - Separate BLAST databases for each segment type

3. **GermlineData paths**:
   - `igdata` → `germlines/igblast/Ig/` (contains internal_data/)
   - `v_gene_dir` → `germlines/igblast/database/{species}/{species}_V`
   - `d_gene_dir` → `germlines/igblast/database/{species}/{species}_D`
   - `j_gene_dir` → `germlines/igblast/database/{species}/{species}_J`
   - `c_gene_dir` → `germlines/igblast/database/{species}/{species}_C`

## Impact

- **complete_vdj** now calculates correctly for all species
- IgBLAST can properly annotate V, D, J, and C regions
- No functional changes to alignment results (same BLAST databases)
- Symlinks replaced with actual files for better portability
