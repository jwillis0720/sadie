# Phase 25 Research: Macaque Germlines Integration

## Summary

**Primary Recommendation:** Build macaque germline databases by combining VDJbase (V/D/J genes) with IMGT GENE-DB (C genes). VDJbase does NOT provide C genes for macaque - this is the critical finding that changes the implementation strategy.

**Key Change Required:** Update VDJbase provider species mapping from `"Rhesus Macaque": "rhesus_macaque"` to `"Rhesus Macaque": "macaque"`.

---

## Standard Stack

| Component | Version/Source | Purpose |
|-----------|---------------|---------|
| VDJbase API | `https://vdjbase.org/admin/api/` | V, D, J genes for macaque |
| IMGT GENE-DB | `IMGTGENEDB-ReferenceSequences.fasta-nt-*` | C genes for macaque |
| GermlinePipeline | `src/sadie/germlines/pipeline.py` | Orchestrates database builds |
| build_internal_data.py | `src/sadie/germlines/scripts/` | Generates NDM files + symlinks |

---

## Architecture Patterns

### Data Flow
```
1. VDJbase API → sources/vdjbase/macaque/ (V, D, J genes)
2. IMGT GENE-DB → sources/imgt/macaque/ or merged (C genes only)
3. GermlineManager → normalized/macaque/ (merged, gapped)
4. BlastDBBuilder → igblast/database/macaque/
5. AuxFileBuilder → igblast/aux_db/macaque_gl.aux
6. build_internal_data.py → igblast/Ig/internal_data/macaque/
```

### Provider Priority
The existing provider priority (`['custom', 'ogrdb', 'vdjbase', 'imgt']`) ensures:
- VDJbase V/D/J genes are used when available
- IMGT C genes fill the gap since VDJbase lacks C genes

---

## VDJbase API Details

### Species Availability
```bash
curl "https://vdjbase.org/admin/api/repseq/species"
# Returns: ["Rhesus Macaque", "Human"]
```

### Macaque Datasets
```bash
curl "https://vdjbase.org/admin/api/repseq/ref_seqs/Rhesus%20Macaque"
# Returns: [{"dataset": "IGH"}, {"dataset": "IGK"}, {"dataset": "IGL"}]
```

### Sequence Retrieval (with pagination)
```bash
curl "https://vdjbase.org/admin/api/repseq/sequences/Rhesus%20Macaque/IGH?per_page=500"
# Returns: V=2147, D=524, J=27, C=0 (NO C GENES!)
```

### Critical Finding: No C Genes in VDJbase
VDJbase repseq endpoint returns gene types: `['IGHD', 'IGHJ', 'IGHV']` only.
**C genes are NOT available in VDJbase for Rhesus Macaque.**

---

## C Gene Source: IMGT GENE-DB

### IMGT Has Macaque C Genes
```bash
# Macaca mulatta C genes available:
IGHA*01 through IGHA*11 (IgA isotypes)
IGHD*01 (IgD)
IGHE*01 (IgE)
IGHG*01-04 (IgG isotypes)
IGHM*01 (IgM)
IGKC (Kappa constant)
IGLC (Lambda constant)
```

### Download C Genes
The existing `download_imgt.py` already handles C gene downloads from GENE-DB:
```python
# In download_imgt.py:
SPECIES_GENEDB_MAP = {
    "rhesus_macaque": "Macaca mulatta",  # ← Already defined
    ...
}
```

### Integration Strategy
1. Download VDJbase macaque → `sources/vdjbase/macaque/`
2. Download IMGT macaque → `sources/imgt/macaque/` (includes C genes)
3. GermlineManager merges both sources during normalization

---

## Don't Hand-Roll

### Use Existing Components
1. **VDJbaseProvider.download()** - Already handles pagination, error recovery
2. **IMGTDownloader._download_c_genes_from_genedb()** - Already fetches C genes
3. **GermlinePipeline.update()** - Change detection + rebuild orchestration
4. **build_internal_data.py** - NDM file generation for IgBLAST

### Key Code Changes (Minimal)
Only ONE line needs to change in `vdjbase.py`:
```python
# Before:
SPECIES_MAP = {
    "Rhesus Macaque": "rhesus_macaque",
}

# After:
SPECIES_MAP = {
    "Rhesus Macaque": "macaque",  # ← Single change!
}
```

---

## Common Pitfalls

### 1. Missing C Genes Will Break Tests
- `test_airr_constant_region_macaque` requires C gene annotation
- Without C genes, the test fails with empty c_call results
- **Solution:** Ensure IMGT C genes are downloaded before VDJbase build

### 2. Species Name Mismatch
- VDJbase API uses "Rhesus Macaque"
- Codebase expects "macaque"
- Tests check for `internal_data/macaque/` not `internal_data/rhesus_macaque/`
- **Solution:** Map VDJbase name to `macaque` in SPECIES_MAP

### 3. Gapping Required for VDJbase Sequences
- VDJbase provides UNGAPPED sequences
- IgBLAST aux files require gapped sequences for FR/CDR positions
- **Solution:** GapperService auto-gaps using human IMGT templates (already implemented)

### 4. NDM Generation Requires Gapped V Genes
- `build_internal_data.py` reads `sources/imgt/{species}/IG*V_gapped.fasta`
- Must ensure gapped V genes exist from IMGT or auto-gapping works

---

## Code Examples

### 1. Download Macaque Data
```bash
# Download VDJbase macaque (V, D, J)
cd /Users/tmsincomb/sadie
python -m sadie.germlines.providers.vdjbase --download macaque

# Download IMGT macaque (includes C genes)
python src/sadie/germlines/scripts/download_imgt.py --species rhesus_macaque
```

### 2. Build Macaque Databases
```bash
# Update/rebuild macaque germlines
python -c "
from pathlib import Path
from sadie.germlines.pipeline import GermlinePipeline
pipeline = GermlinePipeline(Path('src/sadie/germlines'))
pipeline.force_rebuild('macaque')
"

# Generate internal_data
python src/sadie/germlines/scripts/build_internal_data.py macaque
```

### 3. Verify Macaque Availability
```python
from sadie.airr import Airr

# This should work after build:
airr = Airr("macaque")
results = airr.run_fasta("test.fasta")
```

---

## Sources

| Source | Confidence | Notes |
|--------|------------|-------|
| VDJbase API (live query) | HIGH | Verified `/repseq/sequences` returns no C genes |
| IMGT GENE-DB (live query) | HIGH | Confirmed Macaca mulatta IGHA, IGHD, IGHE, IGHG, IGHM available |
| `vdjbase.py` | HIGH | Current implementation reviewed |
| `download_imgt.py` | HIGH | C gene download logic confirmed |
| `pipeline.py` | HIGH | Build orchestration reviewed |
| VDJbase REST API Docs | HIGH | https://wordpress.vdjbase.org/index.php/vdjbase_help/using-the-vdjbase-rest-api/ |
| IMGT Rhesus Macaque Loci | MEDIUM | PMC8950363 paper confirms C gene structure |

---

## Quality Gate Verification

### Pre-Implementation Checklist
- [ ] VDJbase API accessible (verified working)
- [ ] IMGT GENE-DB accessible (verified working)
- [ ] `SPECIES_MAP` change identified in `vdjbase.py` line ~19
- [ ] 6 skipped tests identified in `test_airr.py` and `test_methods.py`

### Post-Implementation Checklist
- [ ] `GermlineData("macaque")` resolves without error
- [ ] `internal_data/macaque/` directory exists with NDM files
- [ ] C genes present in `normalized/macaque/gapped/IGHC.fasta`
- [ ] All 6 macaque tests pass (not skipped)
- [ ] No regressions in human/mouse tests

---

## Implementation Sequence

1. **Update SPECIES_MAP** in `vdjbase.py`:
   ```python
   "Rhesus Macaque": "macaque"
   ```

2. **Download IMGT macaque data** (for C genes):
   ```bash
   python src/sadie/germlines/scripts/download_imgt.py --species rhesus_macaque
   ```

3. **Download VDJbase macaque data** (V, D, J):
   ```python
   from sadie.germlines.providers.vdjbase import VDJbaseProvider
   provider = VDJbaseProvider()
   provider.download(["macaque"])
   ```

4. **Build germline databases**:
   ```python
   from sadie.germlines.pipeline import GermlinePipeline
   pipeline = GermlinePipeline(Path("src/sadie/germlines"))
   pipeline.force_rebuild("macaque")
   ```

5. **Generate internal_data**:
   ```bash
   python src/sadie/germlines/scripts/build_internal_data.py macaque
   ```

6. **Remove skip decorators** from tests

7. **Run tests** to verify
