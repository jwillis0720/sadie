# Phase 32 Research: Fix IgBLAST internal_data to Use Combined VDJC File

## Summary

The germlines module's `internal_data/` directory uses symlinks to separate V/D/J/C files, while IgBLAST expects a combined VDJC file (named `{species}_V.fasta`) for proper `complete_vdj` calculation. This mismatch causes incorrect `complete_vdj` values for all species, leading to the Phase 17 workaround that only works for human.

## Primary Recommendation

**Modify `build_internal_data.py` to create combined VDJC files (not symlinks) matching G3's working structure.** This allows IgBLAST to correctly calculate `complete_vdj` for all species without post-processing workarounds.

## Current Implementation Analysis

### 1. Symlink Creation Script

**File:** `src/sadie/germlines/scripts/build_internal_data.py`

The script creates symlinks from `internal_data/{species}/` to `database/{species}/`:

```python
def build_internal_data(species: str, germlines_root: Path) -> bool:
    database_dir = germlines_root / "igblast" / "database" / species
    internal_data_dir = germlines_root / "igblast" / "Ig" / "internal_data" / species

    # Symlink BLAST database files
    for db_file in database_dir.glob(f"{species}_*"):
        link_path = internal_data_dir / db_file.name
        rel_path = os.path.relpath(db_file, internal_data_dir)
        link_path.symlink_to(rel_path)
```

**Current Result (Broken):**
```
internal_data/human/
├── human.ndm.imgt          (actual file - NDM annotation)
├── human_V.* -> ../../../database/human/human_V.* (symlinks, V-only, 684 seqs)
├── human_D.* -> ../../../database/human/human_D.* (symlinks)
├── human_J.* -> ../../../database/human/human_J.* (symlinks)
└── human_C.* -> ../../../database/human/human_C.* (symlinks)
```

### 2. G3's Working Structure

**Path:** `src/sadie/airr/data/germlines/Ig/internal_data/human/`

G3 has a **combined VDJC file** named `human_V.fasta` with ALL segments:

```
# human_V.fasta content (479 sequences total):
>IGHA1*01
GCATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCAGCACC...
>IGHA2*01
...
>IGHD1-1*01
GGTACAACTGGAACGAC
...
>IGHJ1*01
GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG
...
>IGHV1-2*01
...
```

**G3 Structure (Working):**
```
internal_data/human/
├── human.ndm.imgt          (NDM annotation)
└── human_V.* (BLAST db from combined VDJC file, 479 seqs)
                            NO human_D.*, human_J.*, human_C.*
```

### 3. Phase 17 Workaround

**Files:**
- `src/sadie/airr/airr.py` - `_recalculate_complete_vdj()` method
- `src/sadie/germlines/builders/j_gene_data.py` - `J_GENE_LENGTHS` dictionary

The workaround post-processes IgBLAST results to fix `complete_vdj`:

```python
# j_gene_data.py - Only has HUMAN J alleles (34 total)
J_GENE_LENGTHS = {
    "IGHJ1*01": 52, "IGHJ2*01": 53, "IGHJ3*01": 50, ...
    "IGKJ1*01": 38, "IGKJ2*01": 39, ...
    "IGLJ1*01": 38, "IGLJ2*01": 38, ...
}

# airr.py - Recalculates complete_vdj using hardcoded lengths
def _recalculate_complete_vdj(self, result: AirrTable) -> AirrTable:
    """Called after IgBLAST to fix complete_vdj quirk"""
    expected_j_len = get_j_gene_length(j_allele)  # Only works for human!
    return v_complete and j_complete
```

**Problem:** This workaround only supports human. All other 28+ species fail because their J alleles aren't in `J_GENE_LENGTHS`.

## Scripts to Modify

### 1. Primary: `build_internal_data.py`

**Current behavior:** Creates symlinks to separate V/D/J/C files
**Required change:** Concatenate V+D+J+C FASTAs → single `{species}_V.fasta` → build BLAST db

```python
def build_internal_data(species: str, germlines_root: Path) -> bool:
    database_dir = germlines_root / "igblast" / "database" / species
    internal_data_dir = germlines_root / "igblast" / "Ig" / "internal_data" / species

    # NEW: Concatenate all segment FASTAs into combined file
    combined_fasta = internal_data_dir / f"{species}_V.fasta"
    with open(combined_fasta, "w") as out:
        for segment in ["V", "D", "J", "C"]:
            fasta_file = database_dir / f"{species}_{segment}.fasta"
            if fasta_file.exists():
                with open(fasta_file) as f:
                    out.write(f.read())

    # NEW: Build BLAST db from combined file (actual files, not symlinks)
    write_blast_db(combined_fasta, internal_data_dir / f"{species}_V")

    # KEEP: Generate NDM file (unchanged)
    # ...
```

### 2. Remove Workaround: `airr.py`

**File:** `src/sadie/airr/airr.py`

Remove:
- `_recalculate_complete_vdj()` method (~40 lines)
- Call in `run_fasta()`: `result = self._recalculate_complete_vdj(result)`
- Call in `_run_scfv()`: `_heavy_airr = self._recalculate_complete_vdj(_heavy_airr)`

### 3. Remove Workaround Data: `j_gene_data.py`

**File:** `src/sadie/germlines/builders/j_gene_data.py`

Remove:
- `J_GENE_LENGTHS` dictionary (~34 entries)
- `get_j_gene_length()` function

Keep:
- `HUMAN_J_GENE_DATA` (used for aux file generation)
- `get_j_gene_data()` (used for aux file generation)

### 4. Reference Builder: `reference.py`

**File:** `src/sadie/reference/reference.py`

The `_make_internal_annotaion_file()` method already creates combined files correctly for the reference builder workflow. It creates V-gene filtered BLAST databases:

```python
def _make_internal_annotaion_file(self, outpath: Path):
    # Already creates combined files for custom references
    # May need update to include D/J/C segments
```

**Note:** This method currently only includes V genes. Update to include all segments.

## Species Count

**Total species in `database/`:** 29 directories
- alpaca, atlantic_cod, atlantic_salmon, camel, cat, channel_catfish, chicken
- cow, cynomolgus, dog, ferret, goat, gorilla, horse, human, lemur
- macaque, mink, mouse, mouse_c57bl6j, orangutan_bornean, orangutan_sumatran
- pig, platypus, rabbit, rainbow_trout, rat, rhesus_macaque, sheep, zebrafish

**Current `internal_data/`:** Only 5 species built
- chicken, human, macaque, mouse, rhesus_macaque

## Dependencies and Code Paths

### Files that Read internal_data

| File | Usage |
|------|-------|
| `germline.py` | `GermlineData` class finds BLAST dbs in `internal_data/{species}/` |
| `igblast.py` | `igdata` property points to `Ig/` containing `internal_data/` |
| `reference.py` | `_make_internal_annotaion_file()` creates internal_data for custom references |

### GermlineData Path Logic

```python
# germline.py - looks for databases in internal_data
internal_data_species = germlines_igblast / "Ig" / "internal_data" / name
if internal_data_species.exists():
    self.blast_dir = internal_data_species / f"{name}_"  # e.g., human_
    self.v_gene_dir = Path(self.blast_dir.__str__() + "V")  # human_V
    # Currently expects human_D, human_J, human_C to exist
    # After fix: only human_V will exist (combined)
```

**Impact:** Code expects separate D/J/C paths. After fix, only `_V` path will exist but will contain all segments.

### Tests for internal_data

| Test File | What it Checks |
|-----------|----------------|
| `test_reference_integration.py` | Internal data file structure, NDM files |
| `test_multi_species.py` | Species internal_data directories exist |
| `test_germline_data_legacy.py` | GermlineData path resolution |

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Breaking existing species that use symlinks | Rebuild all species with combined files |
| D/J/C database path references in code | Audit all `_D`, `_J`, `_C` path usages - may need to keep database/ structure for search |
| Test fixtures expect separate files | Update test fixtures to expect combined structure |
| Reference builder creates V-only internal files | Update `_make_internal_annotaion_file()` to include all segments |

## Recommended Implementation Approach

### Step 1: Modify `build_internal_data.py`

1. Remove symlink creation logic
2. Add FASTA concatenation: V + D + J + C → `{species}_V.fasta`
3. Build BLAST database from combined file
4. Keep NDM file generation unchanged

### Step 2: Rebuild All Species

```bash
# Run for all 29 species
python build_internal_data.py alpaca atlantic_cod ... zebrafish
```

### Step 3: Remove Phase 17 Workaround

1. Delete `J_GENE_LENGTHS` from `j_gene_data.py`
2. Delete `get_j_gene_length()` function
3. Delete `_recalculate_complete_vdj()` from `airr.py`
4. Remove workaround calls from `run_fasta()` and `_run_scfv()`

### Step 4: Update Reference Builder

Modify `_make_internal_annotaion_file()` in `reference.py` to include D/J/C segments in the combined file.

### Step 5: Update Tests

1. Fix any tests expecting separate V/D/J/C files in internal_data
2. Add test verifying combined file contains all segments
3. Add test verifying `complete_vdj` works for non-human species

## Verification Checklist

- [ ] No symlinks in any `internal_data/{species}/` directory
- [ ] Each species has `{species}_V.fasta` with ALL segments (V+D+J+C)
- [ ] Each species has `{species}_V.*` BLAST database files
- [ ] `complete_vdj` works correctly for human (current test cases)
- [ ] `complete_vdj` works correctly for mouse, macaque, chicken
- [ ] Phase 17 workaround code completely removed
- [ ] All existing tests pass
- [ ] Reference builder creates correct structure for custom databases

## Sources

| Source | Confidence | Key Finding |
|--------|------------|-------------|
| G3 `human_V.fasta` | HIGH | 479 sequences, all segments combined |
| `build_internal_data.py` | HIGH | Creates symlinks to separate files |
| `audit/igblast-quirk.md` | HIGH | Documents Phase 17 workaround |
| `j_gene_data.py` | HIGH | Only human J alleles (34 entries) |
| IgBLAST documentation | MEDIUM | Expects combined file in internal_data |

---
*Research completed: 2026-01-27*
