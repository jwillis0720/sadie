# Quickstart: Germline Integration

**Feature**: 002-germline-integration
**Audience**: Developers integrating germlines module with SADIE
**Time**: 30-45 minutes

## Prerequisites

- ✅ Germlines module installed and populated (from `001-germline-completion`)
- ✅ Python 3.10+ environment
- ✅ SADIE development environment set up
- ✅ Understanding of AIRR/renumbering workflows

## Overview

This quickstart walks through integrating the germlines module into three SADIE components:
1. **IgBLAST** (AIRR annotation) - ~10 min
2. **HMM/Renumbering** - ~15 min
3. **Reference System** - ~15 min

## Step 1: Verify Germlines Module (5 min)

###1.1 Check Installation

```bash
cd /path/to/sadie
python -c "from sadie.germlines import get_germlines_base_dir; print(get_germlines_base_dir())"
```

**Expected output**: `/path/to/sadie/src/sadie/germlines`

### 1.2 Verify Data Populated

```bash
ls src/sadie/germlines/normalized/human/ungapped/
```

**Expected files**:
- `IGHV.fasta`, `IGHD.fasta`, `IGHJ.fasta`
- `IGKV.fasta`, `IGKJ.fasta`
- `IGLV.fasta`, `IGLJ.fasta`

### 1.3 Build Databases (if not already built)

```python
from sadie.germlines import update_databases

# Build for human (takes ~2-3 minutes)
update_databases("human", force=True)
```

**Verifies**:
- ✅ BLAST databases created in `germlines/igblast/database/human/`
- ✅ Auxiliary files created in `germlines/igblast/aux_db/`
- ✅ HMM models ready for generation

## Step 2: IgBLAST Integration (10 min)

### 2.1 Update GermlineData Class

**File**: `src/sadie/airr/igblast/germline.py`

**Changes**:
```python
def __init__(self, name: str, receptor: str = "Ig", database_dir: Optional[str | Path] = None, scheme: str = "imgt"):
    self.name = name

    # NEW: Check feature flag
    if database_dir:
        self.base_dir = Path(database_dir).absolute()
    elif _use_germlines_module():  # NEW CONDITION
        # Use germlines module paths
        germlines_igblast = _get_germlines_igblast_dir()
        self.base_dir = germlines_igblast
        self.blast_dir = germlines_igblast / "database" / name / f"{name}_"
        self.v_gene_dir = Path(self.blast_dir.__str__() + "V")
        # ... (D, J, C genes similar)
        self.aux_path = germlines_igblast / "aux_db" / f"{name}_gl.aux"
        self.igdata = germlines_igblast / "internal_data"
    else:
        # Legacy G3 paths (existing code)
        self.base_dir = Path(__file__).absolute().parent / "../data/germlines/"
        # ... (existing path setup)
```

### 2.2 Test IgBLAST Integration

```python
import os
from sadie.airr.igblast.germline import GermlineData

# Enable germlines module
os.environ["SADIE_USE_GERMLINES_MODULE"] = "true"

# Test path switching
gd = GermlineData("human")
print(f"Base dir: {gd.base_dir}")
print(f"V gene DB: {gd.v_gene_dir}")
print(f"Aux file: {gd.aux_path}")

# Verify files exist
assert gd.v_gene_dir.with_suffix(".nhr").exists(), "V database not found"
assert gd.aux_path.exists(), "Auxiliary file not found"

print("✅ IgBLAST integration working!")
```

### 2.3 Test with AIRR Annotation

```python
from sadie.airr import Airr

# Create AIRR instance (will use germlines paths)
airr = Airr(species="human")

# Run on test sequence
test_seq = "CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGG..."
result = airr.run([test_seq])

print(f"V gene call: {result['v_call'][0]}")
print("✅ AIRR annotation with germlines backend working!")
```

## Step 3: HMM/Renumbering Integration (15 min)

### 3.1 Create LocalHMMBuilder

**File**: `src/sadie/germlines/renumbering_integration.py` (already created if following along)

**Key class**:
```python
from sadie.germlines.renumbering_integration import LocalHMMBuilder

builder = LocalHMMBuilder()
hmm = builder.get_hmm(species="human", chain="H")
print(f"HMM name: {hmm.name}")  # b'human_H'
```

### 3.2 Update HMMER Class

**File**: `src/sadie/renumbering/aligners/hmmer.py`

**Add at top**:
```python
def _use_local_hmm_builder() -> bool:
    """Check if germlines module should be used for HMM building."""
    try:
        from sadie.germlines.renumbering_integration import use_local_hmm_builder
        return use_local_hmm_builder()
    except ImportError:
        return False
```

**Update get_hmm_models()** (see contract in contracts/integration-api.md for full implementation)

### 3.3 Test Renumbering Integration

```python
import os
from sadie.renumbering.aligners.hmmer import HMMER

# Enable germlines
os.environ["SADIE_USE_GERMLINES_MODULE"] = "true"

# Test HMM loading
hmmer = HMMER(species="human", chains="H")
print(f"Loaded {len(hmmer.hmms)} HMM models")

# Test renumbering
test_seq = "QVQLVQSGAEVKKPGASVKVSCKASGYTFT..."
result = hmmer.align_sequences([test_seq])

print(f"Alignment score: {result[0]['score']}")
print("✅ Renumbering with germlines backend working!")
```

## Step 4: Reference System Integration (15 min)

### 4.1 Create G3 Adapter

**File**: `src/sadie/germlines/g3_adapter.py` (already created if following along)

**Test adapter**:
```python
from sadie.germlines.g3_adapter import GermlineToG3Adapter
from sadie.germlines import get_gene_by_name

adapter = GermlineToG3Adapter()

# Get gene from germlines
gene = get_gene_by_name("human", "IGHV1-69*01")

# Transform to G3 format
g3_dict = adapter.to_g3_format(gene)

print(f"Gene: {g3_dict['gene']}")
print(f"Source: {g3_dict['source']}")
print(f"Functional: {g3_dict['imgt']['imgt_functional']}")
print("✅ G3 adapter working!")
```

### 4.2 Update Reference Class

**File**: `src/sadie/reference/reference.py`

**Add parameter to `__init__`**:
```python
class Reference:
    def __init__(self, endpoint: str = _endpoint, use_germlines: bool = False):
        """
        Initialize reference object.

        Args:
            endpoint: G3 API endpoint (ignored if use_germlines=True)
            use_germlines: Use local germlines module instead of G3 API
        """
        self.data = []
        self.use_germlines = use_germlines
        if not use_germlines:
            self.endpoint = endpoint
        else:
            # Import germlines components
            from sadie.germlines import get_manager
            from sadie.germlines.g3_adapter import GermlineToG3Adapter
            self.germline_manager = get_manager()
            self.g3_adapter = GermlineToG3Adapter()
```

**Update `_get_gene()`**:
```python
def _get_gene(self, gene: GeneEntry) -> Dict[str, str]:
    if self.use_germlines:
        # Query germlines module
        from sadie.germlines import get_gene_by_name
        germline_gene = get_gene_by_name(gene.species, gene.gene)
        if not germline_gene:
            raise G3Error(f"Gene {gene.gene} not found in germlines database")
        # Transform to G3 format
        return self.g3_adapter.to_g3_format(germline_gene)
    else:
        # Use G3 API (existing code)
        # ... existing implementation
```

### 4.3 Test Reference Integration

```python
from sadie.reference import Reference
from sadie.reference.models import GeneEntry

# Create reference with germlines backend
ref = Reference(use_germlines=True)

# Test gene lookup
gene_entry = GeneEntry(species="human", gene="IGHV1-69*01", source="imgt")
ref.add_gene(species="human", gene="IGHV1-69*01", source="imgt")

# Get dataframe
df = ref.get_dataframe()
print(f"Retrieved {len(df)} genes")
print(f"Columns: {df.columns.tolist()}")
print("✅ Reference system with germlines backend working!")
```

## Step 5: Testing & Validation (10 min)

### 5.1 Run Existing Tests (Backwards Compatibility)

```bash
# Should all pass with feature flag disabled
export SADIE_USE_GERMLINES_MODULE=false
pytest tests/unit/airr/test_igblast.py -v
pytest tests/unit/renumbering/test_hmmer.py -v
```

### 5.2 Run Integration Tests

```bash
# Should all pass with feature flag enabled
export SADIE_USE_GERMLINES_MODULE=true
pytest tests/unit/germlines/ -v
```

### 5.3 Validate Results Match

```python
import os
from sadie.airr import Airr

test_seq = "CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGG..."

# Run with G3
os.environ["SADIE_USE_GERMLINES_MODULE"] = "false"
airr_g3 = Airr(species="human")
result_g3 = airr_g3.run([test_seq])

# Run with germlines
os.environ["SADIE_USE_GERMLINES_MODULE"] = "true"
airr_local = Airr(species="human")
result_local = airr_local.run([test_seq])

# Compare
print(f"G3 V gene: {result_g3['v_call'][0]}")
print(f"Local V gene: {result_local['v_call'][0]}")
assert result_g3['v_call'][0] == result_local['v_call'][0], "Mismatch!"
print("✅ Results match between backends!")
```

## Troubleshooting

### Issue: FileNotFoundError when loading databases

**Solution**:
```python
from sadie.germlines import update_databases
update_databases("human", force=True)
```

### Issue: HMM building fails

**Symptoms**: "Insufficient sequences" warning

**Solution**: Check gapped sequences exist:
```bash
ls src/sadie/germlines/normalized/human/gapped/
```

If missing, regenerate gapped sequences:
```python
from sadie.germlines.pipeline import GermlinePipeline
pipeline = GermlinePipeline()
pipeline.normalize("human")
```

### Issue: G3 format validation errors

**Solution**: Check GermlineGene has all required fields:
```python
gene = get_gene_by_name("human", "IGHV1-69*01")
print(f"Has gapped? {gene.sequence_gapped is not None}")
print(f"Has regions? {gene.regions is not None}")
```

If missing, might need to rebuild normalized data.

### Issue: Tests fail with germlines enabled

**Debug steps**:
1. Check feature flag: `echo $SADIE_USE_GERMLINES_MODULE`
2. Verify paths: Print `GermlineData` paths in test
3. Check database files exist
4. Compare with G3 output for same input

## Next Steps

After completing this quickstart:

1. **Create mirrored test suite** (`tests/unit/germlines/`)
2. **Performance benchmarking** (germlines vs G3)
3. **Documentation updates** (user-facing docs)
4. **Migration guide** for existing users

## Summary

You've successfully integrated the germlines module with:
- ✅ IgBLAST (AIRR annotation)
- ✅ HMM generation (renumbering)
- ✅ Reference system (gene queries)

All integrations maintain backwards compatibility and can be toggled via the `SADIE_USE_GERMLINES_MODULE` environment variable.

