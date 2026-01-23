# Phase 11: IMGT Gapped Fix - Execution Plan

## Overview

**Goal**: Fix IMGT provider to load gapped sequences from `*_gapped.fasta` files, enabling HMM generation for all species with gapped data

**Root Cause**: IMGT provider only reads `IGHV.fasta` (ungapped) and ignores `IGHV_gapped.fasta` files. Human works because its main file happens to contain dots (gapped), but rabbit/chicken/etc fail because their main files are ungapped.

**Solution**: Add the same gapped file loading pattern that OGRDB provider uses.

---

## Task Breakdown

### T073: Add `_get_gapped_fasta_path()` method to IMGT provider

**File**: `src/sadie/germlines/providers/imgt.py`

**Changes**:
```python
def _get_gapped_fasta_path(
    self,
    species: str,
    segment: str,
    chain: str
) -> Path:
    """
    Get path to gapped FASTA file.

    Parameters
    ----------
    species : str
        Species name
    segment : str
        Segment type
    chain : str
        Chain type

    Returns
    -------
    Path
        Path to gapped FASTA file
    """
    return self.data_dir / species / f"IG{chain}{segment}_gapped.fasta"
```

**Location**: After `_create_imgt_gene()` method (around line 180)

**Acceptance Criteria**:
- Method exists and returns correct path
- Path format matches existing gapped files: `sources/imgt/{species}/IG{chain}{segment}_gapped.fasta`

---

### T074: Add `_load_gapped_sequences()` method to IMGT provider

**File**: `src/sadie/germlines/providers/imgt.py`

**Dependencies**: T073 (for context, not code dependency)

**Changes**:
```python
def _load_gapped_sequences(self, fasta_path: Path) -> Dict[str, str]:
    """
    Load gapped sequences from FASTA file.

    Parameters
    ----------
    fasta_path : Path
        Path to gapped FASTA file

    Returns
    -------
    Dict[str, str]
        Mapping of gene name to gapped sequence
    """
    gapped = {}
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            # IMGT format: >ACCESSION|GENE_NAME|...
            header_parts = record.id.split("|")
            gene_name = header_parts[1] if len(header_parts) > 1 else header_parts[0]
            gapped[gene_name] = str(record.seq).upper()
    except Exception as e:
        logger.warning(f"Failed to load gapped sequences from {fasta_path}: {e}")
    return gapped
```

**Location**: After `_get_gapped_fasta_path()` method

**Note**: IMGT header format is `>ACCESSION|GENE_NAME|...` (gene name at index 1), different from OGRDB `>GENE_NAME|...` (gene name at index 0)

**Acceptance Criteria**:
- Method loads sequences from gapped fasta files
- Returns dict mapping gene_name -> gapped_sequence
- Handles missing/invalid files gracefully with warning

---

### T075: Update `fetch_genes()` to merge gapped sequences

**File**: `src/sadie/germlines/providers/imgt.py`

**Dependencies**: T073, T074

**Changes to `fetch_genes()`**:

1. Add import at top: `from typing import Dict` (if not present)

2. Update `fetch_genes()`:
```python
def fetch_genes(
    self,
    species: str,
    segment: str,
    chain: str
) -> List[GermlineGene]:
    fasta_path = self.get_fasta_path(species, segment, chain)
    gapped_path = self._get_gapped_fasta_path(species, segment, chain)

    # Guard: file doesn't exist
    if not fasta_path.exists():
        logger.debug(f"No IMGT file: {fasta_path}")
        logger.info(
            "Run download script or add FASTA manually. "
            "See sources/imgt/README.md"
        )
        return []

    # Load gapped sequences if available
    gapped_sequences: Dict[str, str] = {}
    if gapped_path.exists():
        gapped_sequences = self._load_gapped_sequences(gapped_path)
        logger.debug(f"Loaded {len(gapped_sequences)} gapped sequences from {gapped_path}")

    logger.info(f"Loading IMGT FASTA: {fasta_path}")

    genes = self._parse_imgt_fasta(fasta_path, species, segment, chain, gapped_sequences)

    logger.info(f"Loaded {len(genes)} IMGT genes")

    return genes
```

3. Update `_parse_imgt_fasta()` signature and body:
```python
def _parse_imgt_fasta(
    self,
    fasta_path: Path,
    species: str,
    segment: str,
    chain: str,
    gapped_sequences: Optional[Dict[str, str]] = None
) -> List[GermlineGene]:
    genes = []
    gapped_sequences = gapped_sequences or {}
    
    # ... existing parsing code ...
    
    for record in records:
        gene = self._create_imgt_gene(record, species, segment, chain, gapped_sequences)
        if gene:
            genes.append(gene)
    
    return genes
```

4. Update `_create_imgt_gene()` to use gapped_sequences dict:
```python
def _create_imgt_gene(
    self,
    record,
    species: str,
    segment: str,
    chain: str,
    gapped_sequences: Optional[Dict[str, str]] = None
) -> Optional[GermlineGene]:
    # ... existing header parsing ...
    gapped_sequences = gapped_sequences or {}
    
    # Get sequence
    sequence_raw = str(record.seq).upper()
    is_gapped = "." in sequence_raw
    
    if is_gapped:
        # Main file has gapped sequences (e.g., human)
        sequence_gapped = sequence_raw
        sequence_ungapped = sequence_raw.replace(".", "").replace("-", "")
    else:
        # Main file is ungapped, look up from gapped file
        sequence_ungapped = sequence_raw
        sequence_gapped = gapped_sequences.get(gene_name)
    
    # ... rest of gene creation ...
```

**Acceptance Criteria**:
- `fetch_genes()` loads from both ungapped and gapped files
- `_parse_imgt_fasta()` accepts gapped_sequences parameter
- `_create_imgt_gene()` uses gapped_sequences dict when main sequence is ungapped
- Human (already gapped in main file) continues to work
- Rabbit/chicken now get gapped sequences from `*_gapped.fasta` files

---

### T076: Test rabbit HMM generation now works

**Dependencies**: T073, T074, T075

**Test Location**: Create or update test in `tests/unit/germlines/test_multi_species.py`

**Test Code**:
```python
def test_rabbit_hmm_generation():
    """Test that rabbit HMM generation works after gapped fix."""
    from sadie.germlines import GermlineManager
    from sadie.germlines.renumbering_integration import LocalHMMBuilder
    
    manager = GermlineManager(providers=["imgt"])
    builder = LocalHMMBuilder(manager)
    
    # This should now work (was failing before)
    hmm_path = builder.get_hmm("rabbit", "H")
    assert hmm_path.exists(), "Rabbit H chain HMM should be generated"
    
    # Also test K and L chains
    hmm_path_k = builder.get_hmm("rabbit", "K")
    assert hmm_path_k.exists(), "Rabbit K chain HMM should be generated"
    
    hmm_path_l = builder.get_hmm("rabbit", "L")
    assert hmm_path_l.exists(), "Rabbit L chain HMM should be generated"
```

**Acceptance Criteria**:
- Rabbit HMM generation succeeds for H, K, L chains
- No exceptions thrown during HMM building

---

### T077: Test chicken HMM generation now works

**Dependencies**: T073, T074, T075

**Test Location**: `tests/unit/germlines/test_multi_species.py`

**Test Code**:
```python
def test_chicken_hmm_generation():
    """Test that chicken HMM generation works after gapped fix."""
    from sadie.germlines import GermlineManager
    from sadie.germlines.renumbering_integration import LocalHMMBuilder
    
    manager = GermlineManager(providers=["imgt"])
    builder = LocalHMMBuilder(manager)
    
    # Chicken only has H and L chains (no kappa)
    hmm_path = builder.get_hmm("chicken", "H")
    assert hmm_path.exists(), "Chicken H chain HMM should be generated"
    
    hmm_path_l = builder.get_hmm("chicken", "L")
    assert hmm_path_l.exists(), "Chicken L chain HMM should be generated"
```

**Acceptance Criteria**:
- Chicken HMM generation succeeds for H, L chains
- No exceptions thrown during HMM building

---

### T078: Verify all 29 species have gapped sequences loaded

**Dependencies**: T073, T074, T075

**Test Location**: `tests/unit/germlines/test_multi_species.py`

**Test Code**:
```python
def test_all_species_have_gapped_sequences():
    """Verify all 29 species have gapped sequences loaded from IMGT."""
    from sadie.germlines import GermlineManager
    
    manager = GermlineManager(providers=["imgt"])
    
    # All 29 species that have BLAST databases
    species_list = [
        "human", "mouse", "mouse_c57bl6j", "rat",
        "rhesus_macaque", "cynomolgus", "gorilla", "orangutan_sumatran", 
        "orangutan_bornean", "lemur",
        "dog", "cat", "ferret", "mink",
        "rabbit", "pig", "cow", "sheep", "goat", "horse", "alpaca", "camel",
        "chicken",
        "zebrafish", "atlantic_salmon", "rainbow_trout", "atlantic_cod", "channel_catfish",
        "platypus"
    ]
    
    species_with_gapped = []
    for species in species_list:
        genes = manager.get_genes(species, segment="V", chain="H")
        if genes and any(g.sequence_gapped for g in genes):
            species_with_gapped.append(species)
    
    # Most species should have gapped sequences for V genes
    assert len(species_with_gapped) >= 25, f"Expected at least 25 species with gapped V genes, got {len(species_with_gapped)}"
```

**Acceptance Criteria**:
- At least 25 of 29 species have gapped sequences for V genes
- GermlineGene.sequence_gapped is populated

---

### T079: Run full test suite to ensure no regressions

**Dependencies**: T076, T077, T078

**Commands**:
```bash
# Run all germlines tests
pytest tests/unit/germlines/ -v

# Run renumbering tests that use germlines
pytest tests/unit/test_renumbering*.py -v -k germline

# Run AIRR tests
pytest tests/unit/test_airr*.py -v

# Run full test suite
pytest tests/ -v --ignore=tests/integration
```

**Acceptance Criteria**:
- All existing tests pass
- No regressions in germlines, renumbering, or AIRR modules
- New multi-species tests pass

---

## Execution Order

```
T073 (add _get_gapped_fasta_path)
    ↓
T074 (add _load_gapped_sequences)
    ↓
T075 (update fetch_genes)
    ↓
T076, T077, T078 (parallel: rabbit, chicken, all species tests)
    ↓
T079 (full test suite)
```

## Estimated Effort

| Task | Effort | Risk |
|------|--------|------|
| T073 | 5 min | Low |
| T074 | 10 min | Low |
| T075 | 20 min | Medium |
| T076 | 5 min | Low |
| T077 | 5 min | Low |
| T078 | 10 min | Low |
| T079 | 15 min | Low |

**Total**: ~70 minutes

## Rollback Plan

If issues arise, the changes are isolated to:
1. `src/sadie/germlines/providers/imgt.py` - revert to previous version
2. New tests - can be removed

No database migrations or external dependencies involved.

---
*Created: 2026-01-21*
