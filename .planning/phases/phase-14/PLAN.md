# Phase 14: C Region Data Integration

**Goal:** Add C gene constant region data to germlines module sources and IgBLAST databases

**Wave Structure:** 2 waves of parallel work

---

## Plan 14.1: C Gene Data Sources

**Objective:** Update IMGT downloader to fetch C genes from GENE-DB during pipeline build

### Task 14.1.1: Update IMGTDownloader for GENE-DB
**File:** `src/sadie/germlines/scripts/download_imgt.py`

Update the downloader to fetch C genes from IMGT/GENE-DB:

1. Add GENE-DB URL constant:
```python
GENEDB_URL = "https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP"
```

2. Add method `download_c_genes(species: str)`:
   - Download bulk GENE-DB FASTA (cache locally)
   - Parse FASTA headers to filter by species
   - Filter for C genes (names matching `IG[HKL][ACDEFGM]*`)
   - Group by chain (H, K, L)
   - Write to `sources/imgt/{species}/IG{H|K|L}C.fasta`

3. Update main `download()` method to call `download_c_genes()` after V/D/J

**GENE-DB FASTA header format:**
```
>ACCESSION|GENE*ALLELE|SPECIES|FUNCTIONALITY|REGION|...
```
Example: `>M87789|IGHG1*01|Homo sapiens|F|CH1+H+CH2+CH3+CHS|...`

**Filter criteria:**
- Field 3 (species) = "Homo sapiens" (for human)
- Field 2 (gene) matches `IG[HKL][ACDEFGM]*` (C gene pattern)

**Must-haves:**
- Downloads run during normal `IMGTProvider.download()` call
- C gene FASTAs created alongside V/D/J FASTAs
- Works for any species in GENE-DB (human, mouse, etc.)

**Verification:**
```bash
# Run IMGT download (includes C genes now)
python -c "
from pathlib import Path
from sadie.germlines.providers.imgt import IMGTProvider
provider = IMGTProvider(data_dir=Path('src/sadie/germlines/sources/imgt'))
provider.download(['human'], force=True)
"

# Verify C gene files created
ls -la src/sadie/germlines/sources/imgt/human/IG*C.fasta
```

### Task 14.1.2: Update OGRDBDownloader for C Genes
**File:** `src/sadie/germlines/scripts/download_ogrdb.py`

Update the OGRDB downloader to fetch C genes:

1. OGRDB published 105 IGHC sequences (Jana et al. 2025)
2. Check OGRDB API/Zenodo for C gene germline sets
3. Add C gene download logic to existing downloader

**OGRDB C gene source:**
- Check: `https://ogrdb.airr-community.org/` for IGHC germline set
- Or Zenodo archive for C gene sequences

**Must-haves:**
- Downloads C genes during normal `OGRDBProvider.download()` call
- IGHC FASTA created at `sources/ogrdb/human/IGHC.fasta`
- Note: OGRDB only has heavy chain C genes (IGHC), not IGKC/IGLC

**Verification:**
```bash
# Run OGRDB download (includes C genes now)
python -c "
from pathlib import Path
from sadie.germlines.providers.ogrdb import OGRDBProvider
provider = OGRDBProvider(data_dir=Path('src/sadie/germlines/sources/ogrdb'))
provider.download(['human'])
"

# Verify C gene files created
ls -la src/sadie/germlines/sources/ogrdb/human/IG*C.fasta
```

### Task 14.1.3: Update VDJbaseProvider for C Genes
**File:** `src/sadie/germlines/providers/vdjbase.py`

Update VDJbase provider to fetch C genes via existing API:

1. The VDJbase API already supports pagination (existing code)
2. Add "C" segment to the download logic
3. Check if VDJbase API returns C gene sequences

**Implementation:**
1. Update `_download_chain()` to handle C segments
2. Update segment grouping logic to recognize C genes
3. Write `sources/vdjbase/human/IG{H|K|L}C.fasta`

**Must-haves:**
- Downloads C genes during normal `VDJbaseProvider.download()` call
- C gene FASTAs created alongside V/D/J FASTAs

**Verification:**
```bash
# Run VDJbase download (includes C genes now)
python -c "
from pathlib import Path
from sadie.germlines.providers.vdjbase import VDJbaseProvider
provider = VDJbaseProvider(data_dir=Path('src/sadie/germlines/sources/vdjbase'))
provider.download(['human'])
"

# Verify C gene files created
ls -la src/sadie/germlines/sources/vdjbase/human/IG*C.fasta
```

### Task 14.1.4: Update GermlineGene Model
**File:** `src/sadie/germlines/models.py`

Modify `validate_segment` to accept "C" as a valid segment:

```python
@field_validator("segment")
@classmethod
def validate_segment(cls, v: str) -> str:
    """Validate segment is V, D, J, or C."""
    v = v.upper()
    if v not in ["V", "D", "J", "C"]:
        raise ValueError(f"Segment must be V, D, J, or C, got: {v}")
    return v
```

**Must-haves:**
- Add "C" to the list of valid segments
- Update docstring to reflect C segment support

**Verification:**
```python
from sadie.germlines.models import GermlineGene
gene = GermlineGene(name="IGHG1*01", species="human", segment="C", chain="H", sequence="ACGT", source="imgt")
assert gene.segment == "C"
```

---

## Plan 14.2: Pipeline and Builder Updates

**Objective:** Update the pipeline to process C segments and build C gene BLAST databases

### Task 14.2.1: Update IMGT Provider for C Genes
**File:** `src/sadie/germlines/providers/imgt.py`

Modify `fetch_genes` and `fetch_gene_by_name` to handle C segment:

1. Add "C" to segment iteration in `fetch_gene_by_name`:
```python
for segment in ["V", "D", "J", "C"]:
```

2. Add C gene FASTA path support in `get_fasta_path`:
   - Handle `IG{H|K|L}C.fasta` pattern

**Must-haves:**
- `fetch_genes("human", "C", "H")` returns IGHC genes
- `fetch_genes("human", "C", "K")` returns IGKC genes  
- `fetch_genes("human", "C", "L")` returns IGLC genes

**Verification:**
```python
from sadie.germlines.providers.imgt import IMGTProvider
provider = IMGTProvider(data_dir=Path("src/sadie/germlines/sources/imgt"))
genes = provider.fetch_genes("human", "C", "H")
assert len(genes) > 0
```

### Task 14.2.2: Update Pipeline Constants
**File:** `src/sadie/germlines/pipeline.py`

Modify the SEGMENTS constant and normalization logic:

1. Change `SEGMENTS = ["V", "D", "J"]` to `SEGMENTS = ["V", "D", "J", "C"]`

2. Update `_write_gapped_fasta` to handle C genes (C genes don't have IMGT gapping):
```python
gapped_records = [
    SeqRecord(
        seq=Seq(str(gene.sequence_gapped or gene.sequence)),
        id=gene.name,
        description=f"source={gene.source}"
    )
    for gene in genes
    if gene.sequence_gapped or segment in ["D", "C"]  # D and C segments may be ungapped
]
```

**Must-haves:**
- Normalized FASTA files created: `normalized/human/ungapped/IG{H|K|L}C.fasta`
- Pipeline processes C genes without errors

**Verification:**
```bash
# Force rebuild
python -c "from sadie.germlines.pipeline import GermlinePipeline; GermlinePipeline(Path('src/sadie/germlines')).force_rebuild('human')"

# Check normalized files
ls -la src/sadie/germlines/normalized/human/ungapped/IG*C.fasta
```

### Task 14.2.3: Update BLAST Builder for C Segment
**File:** `src/sadie/germlines/builders/blast.py`

1. Change `SEGMENTS = ["V", "D", "J"]` to `SEGMENTS = ["V", "D", "J", "C"]`

**Must-haves:**
- `human_C.fasta` and associated BLAST files created in `igblast/database/human/`
- BLAST database files: `human_C.nhr`, `human_C.nin`, `human_C.nsq`

**Verification:**
```bash
ls -la src/sadie/germlines/igblast/database/human/human_C.*
```

### Task 14.2.4: Copy BLAST DB to internal_data
**File:** `src/sadie/germlines/scripts/build_internal_data.py`

The current `build_internal_data.py` script needs to also copy C gene databases to the IgBLAST internal_data directory:

1. Add C segment to the segments copied
2. Ensure `human_C.*` files are copied to `igblast/Ig/internal_data/human/`

**Alternative:** Modify `GermlineData` in `germline.py` to point to `igblast/database/human/` for C genes if different from V/D/J location.

**Must-haves:**
- C gene BLAST databases accessible at path expected by IgBLAST
- No "C gene directory not found" warning

**Verification:**
```bash
ls -la src/sadie/germlines/igblast/Ig/internal_data/human/human_C.*
```

---

## Plan 14.3: Verification and Audit

**Objective:** Validate C gene integration and parity improvement

### Task 14.3.1: End-to-End Test
Run the full pipeline and verify IgBLAST can use C genes:

```bash
# Rebuild everything
cd /Users/tmsincomb/sadie
python -c "
from pathlib import Path
from sadie.germlines.pipeline import GermlinePipeline
pipeline = GermlinePipeline(Path('src/sadie/germlines'))
pipeline.force_rebuild('human')
"

# Verify no warnings
python -c "
from sadie.airr.igblast.germline import GermlineData
gd = GermlineData('human')
print(f'V: {gd.v_gene_dir}')
print(f'D: {gd.d_gene_dir}')
print(f'J: {gd.j_gene_dir}')
print(f'C: {gd.c_gene_dir}')
"
```

**Success criteria:**
- No warnings about missing C gene directory
- All paths resolve correctly

### Task 14.3.2: Re-run Audit Notebook
Re-execute the audit notebook and verify:

1. All 129 columns present (matching G3)
2. CDR3/junction fields populated
3. `complete_vdj` flag matches G3
4. Parity approaches or reaches 100%

```bash
# Run audit comparison
jupyter nbconvert --execute audit/audit.ipynb --to notebook --inplace
```

---

## Dependency Graph

```
Plan 14.1 (Wave 1)              Plan 14.2 (Wave 2)              Plan 14.3 (Wave 3)
├── Task 14.1.1 (IMGT C)    →   ├── Task 14.2.1 (provider)  →  ├── Task 14.3.1 (e2e test)
├── Task 14.1.2 (OGRDB C)   →   ├── Task 14.2.2 (pipeline)  →  └── Task 14.3.2 (audit)
├── Task 14.1.3 (VDJbase C) →   ├── Task 14.2.3 (blast)
└── Task 14.1.4 (model)     →   └── Task 14.2.4 (internal_data)
```

**Wave 1:** Foundation (Tasks 14.1.1-14.1.4) - Can run in parallel
**Wave 2:** Integration (Tasks 14.2.1-14.2.4) - Depends on Wave 1
**Wave 3:** Verification (Tasks 14.3.1-14.3.2) - Depends on Wave 2

---

## Files to Modify Summary

| File | Change Type | Description |
|------|-------------|-------------|
| `scripts/download_imgt.py` | MODIFY | Add GENE-DB URL, download C genes during build |
| `scripts/download_ogrdb.py` | MODIFY | Add C gene download from OGRDB |
| `providers/vdjbase.py` | MODIFY | Add C segment to API download logic |
| `models.py` | MODIFY | Add "C" to valid segments |
| `providers/imgt.py` | MODIFY | Handle C segment in fetch_genes |
| `providers/ogrdb.py` | MODIFY | Handle C segment in fetch_genes |
| `pipeline.py` | MODIFY | Add "C" to SEGMENTS constant |
| `builders/blast.py` | MODIFY | Add "C" to SEGMENTS constant |
| `scripts/build_internal_data.py` | MODIFY | Copy C gene databases |

---

## Risk Assessment

| Risk | Mitigation |
|------|------------|
| C genes lack IMGT gapping | Use ungapped sequences only (acceptable for IgBLAST) |
| GENE-DB download fails | Cache file locally, retry with exponential backoff |
| Species name mismatch | Map common names to IMGT format (human → "Homo sapiens") |
| IgBLAST path resolution | Test GermlineData after changes |

---

## Acceptance Criteria

1. ✓ No "C gene directory not found" warnings
2. ✓ All 129 columns present in germlines output (matching G3)  
3. ✓ CDR3/junction fields populated for productive sequences
4. ✓ `complete_vdj` flag matches G3 backend
5. ✓ Parity approaches 100%

---

*Created: 2026-01-22*
*Phase: 14 of Milestone v1.1 Audit*
