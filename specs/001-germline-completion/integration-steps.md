# Integration Steps: Germlines Module → Sadie AIRR

## Current State

### What Exists

**Germlines Module (NEW)**
```
src/sadie/germlines/
├── __init__.py          # get_germline_genes(), get_gene_by_name(), GermlineManager
├── manager.py           # Priority-based gene selection
├── pipeline.py          # sources → normalized → igblast
├── providers/           # IMGT, OGRDB, VDJbase, Custom
├── builders/            # BlastDBBuilder, AuxFileBuilder, HMMBuilder, GapperService
├── sources/             # Raw FASTA files
├── normalized/          # Merged gapped/ungapped
└── igblast/             # BLAST databases, aux files, organism.yaml
```

**AIRR Module (EXISTING)**
```
src/sadie/airr/
├── airr.py              # Airr class - main entry point
├── igblast/
│   ├── germline.py      # GermlineData class - paths to BLAST DBs
│   └── igblast.py       # IgBLASTN wrapper
└── data/germlines/      # OLD location for germline data
```

### Integration Points

| Component | Current Source | New Source | Status |
|-----------|---------------|------------|--------|
| GermlineData paths | `airr/data/germlines/` | `germlines/igblast/` | ⚠️ Flag exists, not wired |
| Reference system | G3 API | `germlines.get_gene_by_name()` | ⚠️ Needs adapter |
| HMM builder | External | `germlines.builders.HMMBuilder` | ⚠️ Needs wiring |

---

## Integration Options

### Option A: Path Switching (Minimal Change)

Modify `GermlineData.__init__` to check feature flag and use new paths.

**Pros**: Minimal code change, backward compatible
**Cons**: Two parallel data sources, potential drift

```python
# In src/sadie/airr/igblast/germline.py
def __init__(self, name, receptor="Ig", database_dir=None, scheme="imgt"):
    self.name = name

    if database_dir:
        self.base_dir = Path(database_dir).absolute()
    elif _use_germlines_module():
        self.base_dir = _get_germlines_igblast_dir()
    else:
        self.base_dir = Path(__file__).parent / "../data/germlines/"

    # Rest of path setup...
```

**Changes Required**:
1. `germline.py`: Add conditional in `__init__`
2. Ensure `germlines/igblast/` structure matches expected paths

### Option B: Delegate to Germlines Module

Have `GermlineData` delegate entirely to germlines module when flag is true.

**Pros**: Single source of truth, cleaner long-term
**Cons**: More code change, needs careful API matching

```python
# In src/sadie/airr/igblast/germline.py
class GermlineData:
    def __init__(self, name, receptor="Ig", database_dir=None, scheme="imgt"):
        self.name = name

        if _use_germlines_module() and database_dir is None:
            self._init_from_germlines_module(name, receptor, scheme)
        else:
            self._init_legacy(name, receptor, database_dir, scheme)

    def _init_from_germlines_module(self, name, receptor, scheme):
        from sadie.germlines import get_germlines_base_dir
        base = get_germlines_base_dir() / "igblast"

        self.base_dir = base
        self.blast_dir = base / "database" / name / f"{name}_"
        self.v_gene_dir = Path(str(self.blast_dir) + "V")
        self.d_gene_dir = Path(str(self.blast_dir) + "D")
        self.j_gene_dir = Path(str(self.blast_dir) + "J")
        self.aux_path = base / "aux_db" / f"{name}_gl.aux"
        self.igdata = base
```

### Option C: Facade Pattern (Recommended)

Create a facade that presents unified API, hiding implementation details.

**Pros**: Clean separation, easy to swap implementations
**Cons**: Additional abstraction layer

```python
# New file: src/sadie/airr/igblast/germline_facade.py
class GermlineFacade:
    """Unified interface for germline data access."""

    def __init__(self, species: str, receptor: str = "Ig"):
        self.species = species
        self.receptor = receptor
        self._use_new = _use_germlines_module()

    @property
    def v_gene_db(self) -> Path:
        if self._use_new:
            from sadie.germlines import get_germlines_base_dir
            return get_germlines_base_dir() / "igblast/database" / self.species / f"{self.species}_V"
        else:
            return self._legacy_path("V")

    # Similar for d_gene_db, j_gene_db, aux_path...
```

---

## Directory Structure Mapping

### Current (airr/data/germlines/)
```
airr/data/germlines/
├── Ig/
│   ├── blastdb/
│   │   └── human/
│   │       └── human_V.nhr, .nin, .nsq
│   └── internal_data/
│       └── human/
└── aux_db/
    └── imgt/
        └── human_gl.aux
```

### New (germlines/igblast/)
```
germlines/igblast/
├── database/
│   └── human/
│       └── human_V.nhr, .nin, .nsq, human_V.fasta
├── aux_db/
│   └── human_gl.aux
└── internal_data/
    └── organism.yaml
```

### Path Translation Table

| GermlineData Attribute | Old Path | New Path |
|------------------------|----------|----------|
| `base_dir` | `airr/data/germlines/` | `germlines/igblast/` |
| `blast_dir` | `{base}/Ig/blastdb/{species}/{species}_` | `{base}/database/{species}/{species}_` |
| `v_gene_dir` | `{blast_dir}V` | `{blast_dir}V` |
| `d_gene_dir` | `{blast_dir}D` | `{blast_dir}D` |
| `j_gene_dir` | `{blast_dir}J` | `{blast_dir}J` |
| `aux_path` | `{base}/aux_db/{scheme}/{species}_gl.aux` | `{base}/aux_db/{species}_gl.aux` |
| `igdata` | `{base}/Ig/` | `{base}/` |

**Key Differences**:
1. No `Ig/` subdirectory in new structure
2. `blastdb/` renamed to `database/`
3. `aux_db/` no longer has scheme subdirectory
4. `internal_data/` uses `organism.yaml` instead of species subdirs

---

## Step-by-Step Integration Plan

### Phase 1: Verify Data Compatibility

1. **Check BLAST DB format matches**
   ```bash
   # Compare file outputs
   diff <(ls germlines/igblast/database/human/) <(ls airr/data/germlines/Ig/blastdb/human/)
   ```

2. **Verify aux file format**
   ```bash
   head -5 germlines/igblast/aux_db/human_gl.aux
   head -5 airr/data/germlines/aux_db/imgt/human_gl.aux
   ```

3. **Test IgBLAST with new paths manually**
   ```bash
   igblastn -germline_db_V germlines/igblast/database/human/human_V \
            -auxiliary_data germlines/igblast/aux_db/human_gl.aux \
            -query test.fasta
   ```

### Phase 2: Implement Path Switching

1. **Update `germline.py`** with Option A or B
2. **Add integration test**
   ```python
   def test_germline_data_uses_new_paths():
       import os
       os.environ["SADIE_USE_GERMLINES_MODULE"] = "true"
       gd = GermlineData("human")
       assert "germlines/igblast" in str(gd.base_dir)
   ```

### Phase 3: Verify Full Pipeline

1. **Run AIRR with new germlines**
   ```python
   from sadie.airr import Airr
   airr = Airr("human")
   result = airr.run_single("test", "CAGGTGCAGCTGGTGCAGTCTGGGGCT...")
   ```

2. **Compare results with old germlines**
   ```python
   # With old
   os.environ["SADIE_USE_GERMLINES_MODULE"] = "false"
   result_old = Airr("human").run_single("test", seq)

   # With new
   os.environ["SADIE_USE_GERMLINES_MODULE"] = "true"
   result_new = Airr("human").run_single("test", seq)

   # Compare V gene calls, CDR3, etc.
   assert result_old["v_call"] == result_new["v_call"]
   ```

### Phase 4: Migration & Cleanup

1. **Deprecation notice** (after validation period)
2. **Remove old data directory**
3. **Remove feature flag code**
4. **Update documentation**

---

## Risks & Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Path structure mismatch | IgBLAST fails | Verify structure in Phase 1 |
| Aux file format differs | Wrong CDR annotations | Compare formats, regenerate if needed |
| Missing species data | Runtime errors | Add data validation on startup |
| Performance regression | Slower annotation | Benchmark before/after |

---

## Testing Checklist

- [ ] GermlineData resolves correct paths with flag=true
- [ ] GermlineData resolves correct paths with flag=false
- [ ] IgBLAST runs successfully with new BLAST DBs
- [ ] Aux file provides correct CDR/FWR boundaries
- [ ] Airr.run_single() produces same results
- [ ] Airr.run_fasta() produces same results
- [ ] Custom germlines are picked up
- [ ] Priority ordering affects results

---

## Quick Reference: Code Locations

| File | Purpose | Changes Needed |
|------|---------|----------------|
| `airr/igblast/germline.py` | GermlineData class | Add path switching logic |
| `airr/igblast/igblast.py` | IgBLASTN wrapper | None (uses GermlineData) |
| `airr/airr.py` | Main Airr class | None (uses GermlineData) |
| `germlines/pipeline.py` | Build BLAST DBs | Ensure output matches expected structure |
| `germlines/builders/blast.py` | makeblastdb wrapper | Verify output path/naming |
| `germlines/builders/aux.py` | Aux file generator | Verify format matches IgBLAST expectations |
