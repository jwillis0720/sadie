# Phase 20: Integration Foundation

## Goal
Route reference.yml processing through germlines module with explicit source selection

## Context Analysis

### Current State
1. `Reference` class already has `use_germlines` parameter in `__init__()` (line 41)
2. `References.from_yaml()` does NOT have `use_germlines` parameter (line 313)
3. `GermlineToG3Adapter.to_g3_format()` does NOT generate `_id` field
4. `_id` field is required for deduplication in `References.get_dataframe()` (lines 280-282)
5. YAML structure explicitly defines source: `name: source: species: [genes]`
6. GermlineManager uses priority-based lookup, but we need explicit source (no fallback)

### Target Files
- `src/sadie/reference/reference.py`
- `src/sadie/germlines/g3_adapter.py`

---

## Task 1: Add `_id` Field Generation to G3 Adapter

**File**: `src/sadie/germlines/g3_adapter.py`

**Requirement**: INT-03 — Generate synthetic `_id` field using deterministic hash

### Changes

#### 1.1 Add hashlib import (top of file)
```python
import hashlib
```

#### 1.2 Add `_generate_id()` method to `GermlineToG3Adapter` class
```python
def _generate_id(self, source: str, species: str, gene_name: str) -> str:
    """
    Generate deterministic _id for gene deduplication.

    Uses SHA-256 hash of source:species:gene to match G3 API behavior.

    Parameters
    ----------
    source : str
        Provider source (e.g., "imgt", "ogrdb")
    species : str
        Species name (e.g., "human")
    gene_name : str
        Gene name (e.g., "IGHV1-69*01")

    Returns
    -------
    str
        Hex digest of hash (first 24 chars for readability)
    """
    key = f"{source}:{species}:{gene_name}"
    return hashlib.sha256(key.encode()).hexdigest()[:24]
```

#### 1.3 Modify `to_g3_format()` to include `_id` field
Add after building base structure (line ~68, before return):
```python
# Generate deterministic _id for deduplication
g3_dict["_id"] = self._generate_id(gene.source, gene.species, gene.name)
```

### Verification
- Unit test: adapter generates consistent `_id` for same gene
- Unit test: different genes produce different `_id` values

---

## Task 2: Add `use_germlines` Parameter to `References.from_yaml()`

**File**: `src/sadie/reference/reference.py`

**Requirement**: INT-01 — Add `use_germlines=True` parameter to `References.from_yaml()`

### Changes

#### 2.1 Modify `from_yaml()` signature and docstring (line 313)
```python
@staticmethod
def from_yaml(yaml_path: Optional[Path] = None, use_germlines: bool = False) -> "References":
    """Parse a yaml file into a references file object

    Parameters
    ----------
    yaml_path : Path
        Path to yaml file
    use_germlines : bool, optional
        If True, use local germlines module instead of G3 API. Defaults to False.

    Returns
    -------
    Reference - Reference Object
    """
```

#### 2.2 Pass `use_germlines` to `Reference()` constructor (line 325)
Change:
```python
reference_object = Reference()
```
To:
```python
reference_object = Reference(use_germlines=use_germlines)
```

### Verification
- `References.from_yaml(use_germlines=True)` creates Reference objects with germlines enabled
- `References.from_yaml(use_germlines=False)` still works (backward compatible)

---

## Task 3: Route Explicit Source to GermlineManager

**File**: `src/sadie/reference/reference.py`

**Requirement**: INT-02 — Route source selection through GermlineManager (explicit source, no priority)

### Changes

#### 3.1 Modify `_get_gene()` germlines path (around line 149-160)
Change the germlines branch to use explicit provider:
```python
# Use germlines module if enabled
if self.use_germlines:
    from sadie.germlines import GermlineManager
    from sadie.germlines.g3_adapter import GermlineToG3Adapter

    # Create manager with explicit source (no priority fallback)
    manager = GermlineManager(providers=[gene.source])
    germline_gene = manager.get_gene_by_name(gene.gene, gene.species)

    if not germline_gene:
        raise G3Error(f"Gene {gene.gene} not found in {gene.source} database for species {gene.species}")

    # Transform to G3 format
    adapter = GermlineToG3Adapter()
    g3_dict = adapter.to_g3_format(germline_gene)
    logger.debug(f"Retrieved {gene.gene} from germlines module ({gene.source})")
    return g3_dict
```

#### 3.2 Modify `_get_genes()` germlines path (around line 198-215)
Change the germlines branch to use explicit provider:
```python
# Use germlines module if enabled
if self.use_germlines:
    from sadie.germlines import GermlineManager
    from sadie.germlines.g3_adapter import GermlineToG3Adapter

    # Create manager with explicit source (no priority fallback)
    manager = GermlineManager(providers=[genes.source])
    adapter = GermlineToG3Adapter()

    results = []
    for gene_name in genes.genes:
        germline_gene = manager.get_gene_by_name(gene_name, genes.species)
        if germline_gene:
            g3_dict = adapter.to_g3_format(germline_gene)
            results.append(g3_dict)
        else:
            logger.warning(f"Gene {gene_name} not found in {genes.source} database for {genes.species}")

    logger.debug(f"Retrieved {len(results)} genes from germlines module ({genes.source})")
    return results
```

#### 3.3 Remove redundant imports from `__init__` (cleanup)
Remove the now-unnecessary early imports in `__init__()` (lines 58-61):
```python
if use_germlines:
    # Germlines components imported in _get_gene/_get_genes when needed
    self._endpoint = endpoint
```

### Verification
- Source from YAML explicitly passed to GermlineManager (not using priority)
- When `source=imgt`, only IMGT provider is queried
- Error message includes source name when gene not found

---

## Execution Order

1. **Task 1** (g3_adapter.py) — No dependencies, can be done first
2. **Task 2** (reference.py) — Depends on Task 1 for `_id` field
3. **Task 3** (reference.py) — Depends on Task 1 and 2

## Wave Assignment
- **Wave 1**: Task 1 (isolated adapter change)
- **Wave 2**: Tasks 2 + 3 (reference.py changes, can be done together)

---

## Test Strategy

### Unit Tests
1. `test_g3_adapter_generates_id`: Verify `_id` field exists and is deterministic
2. `test_references_from_yaml_use_germlines`: Verify parameter flows through
3. `test_explicit_source_selection`: Verify only specified provider is queried

### Integration Tests
1. `test_full_yaml_processing_with_germlines`: Load reference.yml with `use_germlines=True`
2. `test_backward_compatibility`: Ensure `use_germlines=False` still works with G3 API

---

## Success Criteria Checklist

- [ ] `References.from_yaml(use_germlines=True)` loads genes from germlines module
- [ ] Source field from YAML explicitly passed to GermlineManager (not using priority fallback)
- [ ] All returned gene dicts contain `_id` field (hash of `source:species:gene`)
- [ ] Downstream code using `_id` for deduplication/indexing works correctly
- [ ] G3 API path still works with `use_germlines=False` for backwards compatibility
