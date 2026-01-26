# Phase 25: Macaque Germlines Integration — Executable Plan

## Goal Statement

Build macaque germline databases from VDJbase (V/D/J genes) + IMGT (C genes) to enable 6 currently skipped tests.

**Key Finding:** VDJbase has NO C genes for macaque. IMGT GENE-DB provides Macaca mulatta C genes (IGHA, IGHD, IGHE, IGHG, IGHM, IGKC, IGLC). Both sources are required.

## Context Summary

- **6 tests skipped** with `@skip_no_macaque` decorator checking for `internal_data/macaque/`
- **Single code change needed:** `vdjbase.py` SPECIES_MAP: `"Rhesus Macaque": "macaque"` (currently `rhesus_macaque`)
- **Second code change:** `download_imgt.py` needs `macaque` alias for `Macaca mulatta`
- **Build produces:** BLAST databases, NDM files, aux files in `igblast/` structure
- **Tests validate:** C gene annotation, IGL sequences, 5'/3' prime extensions

---

## Plan 1: Provider Configuration Updates

**Objective:** Configure VDJbase and IMGT providers to use `macaque` as the canonical species name.

### Task 1.1: Update VDJbase SPECIES_MAP

**Action:** Change species mapping from `rhesus_macaque` to `macaque`

**File:** `src/sadie/germlines/providers/vdjbase.py`

**Code Change (line ~44):**
```python
# Before:
SPECIES_MAP = {
    "Human": "human",
    "Rhesus Macaque": "rhesus_macaque",
    "Mouse": "mouse",
}

# After:
SPECIES_MAP = {
    "Human": "human",
    "Rhesus Macaque": "macaque",  # Changed from rhesus_macaque
    "Mouse": "mouse",
}
```

**Expected Output:** VDJbase provider saves macaque data to `sources/vdjbase/macaque/`

**Verification:**
```bash
grep -n '"Rhesus Macaque"' src/sadie/germlines/providers/vdjbase.py
# Should show: "Rhesus Macaque": "macaque"
```

### Task 1.2: Add Macaque Alias to IMGT Download Script

**Action:** Add `macaque` as an alias for `Macaca mulatta` in both mapping dictionaries

**File:** `src/sadie/germlines/scripts/download_imgt.py`

**Code Change #1 (SPECIES_GENEDB_MAP, ~line 71):**
```python
SPECIES_GENEDB_MAP = {
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
    # ... existing entries ...
    "rhesus_macaque": "Macaca mulatta",
    "macaque": "Macaca mulatta",  # ADD THIS LINE - alias for canonical name
    # ... rest of dict ...
}
```

**Code Change #2 (SPECIES_MAP, ~line 91):**
```python
SPECIES_MAP = {
    "human": "Homo_sapiens",
    "mouse": "Mus_musculus",
    # ... existing entries ...
    "rhesus_macaque": "Macaca_mulatta",
    "macaque": "Macaca_mulatta",  # ADD THIS LINE - alias for canonical name
    # ... rest of dict ...
}
```

**Expected Output:** IMGT downloader saves macaque data to `sources/imgt/macaque/`

**Verification:**
```bash
grep -n '"macaque"' src/sadie/germlines/scripts/download_imgt.py
# Should show two lines with "macaque": "Macaca..." mappings
```

---

## Plan 2: Download Germline Data

**Objective:** Fetch macaque germline sequences from IMGT (for C genes) and VDJbase (for V/D/J genes).

**Depends on:** Plan 1 (provider configuration)

### Task 2.1: Download IMGT Macaque Data

**Action:** Download IMGT reference sequences for Macaca mulatta (includes C genes)

**Commands:**
```bash
cd /Users/tmsincomb/sadie
python src/sadie/germlines/scripts/download_imgt.py --species macaque
```

**Expected Outputs:**
- `src/sadie/germlines/sources/imgt/macaque/IGHV.fasta`
- `src/sadie/germlines/sources/imgt/macaque/IGHV_gapped.fasta`
- `src/sadie/germlines/sources/imgt/macaque/IGHD.fasta`
- `src/sadie/germlines/sources/imgt/macaque/IGHJ.fasta`
- `src/sadie/germlines/sources/imgt/macaque/IGHC.fasta` (C genes!)
- `src/sadie/germlines/sources/imgt/macaque/IGK*.fasta`
- `src/sadie/germlines/sources/imgt/macaque/IGL*.fasta`

**Verification:**
```bash
ls -la src/sadie/germlines/sources/imgt/macaque/
# Should list FASTA files including IGHC.fasta, IGKC.fasta, IGLC.fasta
```

### Task 2.2: Download VDJbase Macaque Data

**Action:** Download VDJbase reference sequences for Rhesus Macaque (V/D/J genes only)

**Commands:**
```bash
cd /Users/tmsincomb/sadie
python -c "
from sadie.germlines.providers.vdjbase import VDJbaseProvider
provider = VDJbaseProvider()
provider.download(['macaque'])
"
```

**Expected Outputs:**
- `src/sadie/germlines/sources/vdjbase/macaque/IGHV.fasta`
- `src/sadie/germlines/sources/vdjbase/macaque/IGHD.fasta`
- `src/sadie/germlines/sources/vdjbase/macaque/IGHJ.fasta`
- `src/sadie/germlines/sources/vdjbase/macaque/IGKV.fasta`
- `src/sadie/germlines/sources/vdjbase/macaque/IGKJ.fasta`
- `src/sadie/germlines/sources/vdjbase/macaque/IGLV.fasta`
- `src/sadie/germlines/sources/vdjbase/macaque/IGLJ.fasta`

**Verification:**
```bash
ls -la src/sadie/germlines/sources/vdjbase/macaque/
# Should list V, D, J FASTA files (NO C files - VDJbase lacks C genes)
```

---

## Plan 3: Build Germline Databases

**Objective:** Generate IgBLAST databases, internal_data, and aux files for macaque.

**Depends on:** Plan 2 (downloaded data)

### Task 3.1: Run Germline Pipeline Force Rebuild

**Action:** Normalize and build IgBLAST BLAST databases for macaque

**Commands:**
```bash
cd /Users/tmsincomb/sadie
python -c "
from pathlib import Path
from sadie.germlines.pipeline import GermlinePipeline

pipeline = GermlinePipeline(Path('src/sadie/germlines'))
pipeline.force_rebuild('macaque')
"
```

**Expected Outputs:**
- `src/sadie/germlines/normalized/macaque/gapped/*.fasta`
- `src/sadie/germlines/normalized/macaque/ungapped/*.fasta`
- `src/sadie/germlines/igblast/database/macaque/macaque_V.n*` (BLAST db files)
- `src/sadie/germlines/igblast/database/macaque/macaque_D.n*`
- `src/sadie/germlines/igblast/database/macaque/macaque_J.n*`
- `src/sadie/germlines/igblast/database/macaque/macaque_C.n*`
- `src/sadie/germlines/igblast/aux_db/macaque_gl.aux`

**Verification:**
```bash
ls -la src/sadie/germlines/igblast/database/macaque/
ls -la src/sadie/germlines/igblast/aux_db/macaque*
# Should show BLAST database files and aux file
```

### Task 3.2: Generate Internal Data (NDM Files)

**Action:** Create internal_data directory with symlinks and NDM files

**Commands:**
```bash
cd /Users/tmsincomb/sadie
python src/sadie/germlines/scripts/build_internal_data.py macaque
```

**Expected Outputs:**
- `src/sadie/germlines/igblast/Ig/internal_data/macaque/` directory
- `src/sadie/germlines/igblast/Ig/internal_data/macaque/macaque.ndm.imgt`
- Symlinks to BLAST database files

**Verification:**
```bash
ls -la src/sadie/germlines/igblast/Ig/internal_data/macaque/
cat src/sadie/germlines/igblast/Ig/internal_data/macaque/macaque.ndm.imgt | head -5
# Should show NDM entries with FR/CDR positions
```

### Task 3.3: Verify GermlineData Resolution

**Action:** Confirm GermlineData can resolve macaque species

**Commands:**
```bash
cd /Users/tmsincomb/sadie
python -c "
from sadie.airr.igblast.germline import GermlineData
gd = GermlineData('macaque')
print(f'Database dir: {gd.database_dir}')
print(f'Aux path: {gd.aux_path}')
print(f'Internal data: {gd.internal_data_dir}')
print('SUCCESS: GermlineData resolves macaque')
"
```

**Expected Output:**
```
Database dir: .../igblast/database/macaque
Aux path: .../igblast/aux_db/macaque_gl.aux
Internal data: .../igblast/Ig/internal_data/macaque
SUCCESS: GermlineData resolves macaque
```

---

## Plan 4: Enable Tests

**Objective:** Remove skip markers and verify all 6 tests pass.

**Depends on:** Plan 3 (databases built)

### Task 4.1: Verify Test Infrastructure Detects Macaque

**Action:** Confirm `_macaque_available()` returns True

**Commands:**
```bash
cd /Users/tmsincomb/sadie
python -c "
from sadie.germlines import get_germlines_base_dir
macaque_path = get_germlines_base_dir() / 'igblast' / 'Ig' / 'internal_data' / 'macaque'
print(f'Path: {macaque_path}')
print(f'Exists: {macaque_path.exists()}')
assert macaque_path.exists(), 'macaque internal_data not found!'
print('SUCCESS: Macaque germlines available')
"
```

**Expected Output:**
```
Path: .../igblast/Ig/internal_data/macaque
Exists: True
SUCCESS: Macaque germlines available
```

### Task 4.2: Run Macaque Tests (with skip markers)

**Action:** Run tests to confirm they now execute (not skipped)

**Commands:**
```bash
cd /Users/tmsincomb/sadie
pytest tests/unit/airr/test_airr.py::test_five_and_three_prime_extension \
       tests/unit/airr/test_airr.py::test_hard_igl_seqs \
       tests/unit/airr/test_airr.py::test_hard_igl_seqs_linked \
       tests/unit/airr/test_airr.py::test_airr_constant_region_macaque \
       tests/unit/airr/test_methods.py::test_run_five_prime_buffer \
       tests/unit/airr/test_methods.py::test_run_three_prime_buffer \
       -v 2>&1 | grep -E "(PASSED|FAILED|SKIPPED|ERROR)"
```

**Expected Output:** All 6 tests show PASSED (not SKIPPED)

**Verification:** If any test fails, investigate and fix. Common issues:
- Missing C genes → Verify IMGT download included IGHC.fasta
- NDM parse errors → Check gapped sequences in `sources/imgt/macaque/`
- Annotation failures → Run with `-s` flag to see IgBLAST output

### Task 4.3: Run Full Test Suite

**Action:** Verify no regressions in other tests

**Commands:**
```bash
cd /Users/tmsincomb/sadie
pytest tests/unit/airr/ tests/unit/germlines/ -v --tb=short
```

**Expected Output:** All tests pass, no new failures introduced.

---

## Dependencies Between Plans

```
Plan 1 (Provider Config) 
    ↓
Plan 2 (Download Data) ← Requires Plan 1 for correct directory naming
    ↓
Plan 3 (Build Databases) ← Requires Plan 2 for source data
    ↓
Plan 4 (Enable Tests) ← Requires Plan 3 for functional databases
```

---

## Success Criteria (Phase Level)

1. ✅ `GermlineData("macaque")` resolves without error
2. ✅ All 6 macaque tests pass (not skipped):
   - `test_five_and_three_prime_extension`
   - `test_hard_igl_seqs`
   - `test_hard_igl_seqs_linked`
   - `test_airr_constant_region_macaque`
   - `test_run_five_prime_buffer`
   - `test_run_three_prime_buffer`
3. ✅ Macaque annotation produces valid AIRR output with C gene calls
4. ✅ No regressions in human/mouse tests

---

## Files Modified

### Code Changes (2 files)
- `src/sadie/germlines/providers/vdjbase.py` — SPECIES_MAP: `"Rhesus Macaque": "macaque"`
- `src/sadie/germlines/scripts/download_imgt.py` — Add `"macaque"` alias in both mapping dicts

### Generated Files (created by build process)
- `src/sadie/germlines/sources/imgt/macaque/` — IMGT FASTA files
- `src/sadie/germlines/sources/vdjbase/macaque/` — VDJbase FASTA files
- `src/sadie/germlines/normalized/macaque/` — Normalized sequences
- `src/sadie/germlines/igblast/database/macaque/` — BLAST databases
- `src/sadie/germlines/igblast/aux_db/macaque_gl.aux` — Aux file
- `src/sadie/germlines/igblast/Ig/internal_data/macaque/` — NDM files + symlinks

---

## Risk Assessment

### Low Risk
- **VDJbase API availability:** API has been stable; download includes retry logic
- **IMGT data availability:** IMGT GENE-DB is stable; C genes confirmed present

### Medium Risk
- **Test fixture updates may be needed:** VDJbase may have different allele names than original fixtures
  - **Mitigation:** Run regenerate scripts if test assertions fail on gene names

### Mitigated
- **No C genes from VDJbase:** Confirmed by research; using IMGT GENE-DB as source
- **Species name mismatch:** Single line change in `vdjbase.py` resolves

---

## Verification Checklist

```bash
# Pre-build checks
[ ] vdjbase.py SPECIES_MAP updated: "Rhesus Macaque": "macaque"
[ ] download_imgt.py has "macaque" in both SPECIES_MAP and SPECIES_GENEDB_MAP

# Data download checks
[ ] sources/imgt/macaque/ exists with IGHC.fasta
[ ] sources/vdjbase/macaque/ exists with V/D/J files

# Build checks
[ ] igblast/database/macaque/ has BLAST db files
[ ] igblast/aux_db/macaque_gl.aux exists
[ ] igblast/Ig/internal_data/macaque/macaque.ndm.imgt exists

# Runtime checks
[ ] GermlineData("macaque") resolves without error
[ ] All 6 tests pass (not skipped)
[ ] No regressions in human/mouse tests
```
