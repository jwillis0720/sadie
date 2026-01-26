# Architecture Research: Reference Module Unification

## Executive Summary

The Reference module already has foundational integration with Germlines module via the `use_germlines` flag. The core adapter (`GermlineToG3Adapter`) exists and transforms data correctly. The primary gap is propagating the `use_germlines` flag through the YAML-based initialization path and expanding source validation.

**Recommendation**: Extend existing integration points rather than creating new abstractions.

---

## Current Architecture

### Reference Module (`src/sadie/reference/`)

```
reference/
├── reference.py    # Reference, References classes
├── models.py       # GeneEntry, GeneEntries (Pydantic validation)
├── yaml.py         # YamlRef - parses reference.yml
├── util.py         # Blast DB utilities
└── data/
    └── reference.yml  # Gene definitions per species/source
```

**Key Classes:**
- `Reference` - Single reference database (list of genes)
- `References` - Collection of named Reference objects
- `GeneEntry` / `GeneEntries` - Validated gene lookup models
- `YamlRef` - YAML parser with gene lookup methods

### Germlines Module (`src/sadie/germlines/`)

```
germlines/
├── __init__.py     # Public API: get_germline_genes, get_gene_by_name
├── manager.py      # GermlineManager - priority-based lookup
├── models.py       # GermlineGene, ProviderMetadata
├── g3_adapter.py   # GermlineToG3Adapter - format conversion
├── pipeline.py     # Database build pipeline
├── providers/
│   ├── base.py     # GermlineProvider ABC
│   ├── imgt.py     # IMGT provider
│   ├── ogrdb.py    # OGRDB provider
│   ├── vdjbase.py  # VDJbase provider
│   └── custom.py   # Custom provider
└── sources/        # Raw germline data by provider
```

**Key Classes:**
- `GermlineManager` - Multi-provider lookup with priority (custom > ogrdb > vdjbase > imgt)
- `GermlineGene` - Unified gene model across all providers
- `GermlineToG3Adapter` - Transforms GermlineGene → G3 API dict format

---

## Integration Points

### 1. Already Integrated (Working)

| Component | Location | Status |
|-----------|----------|--------|
| `Reference.__init__(use_germlines)` | reference.py:46-59 | ✅ Complete |
| `Reference._get_gene()` germlines path | reference.py:149-162 | ✅ Complete |
| `Reference._get_genes()` germlines path | reference.py:184-203 | ✅ Complete |
| `GermlineToG3Adapter.to_g3_format()` | g3_adapter.py:46-77 | ✅ Complete |

**Existing integration code in `Reference.__init__()`:**
```python
def __init__(self, endpoint: str = _endpoint, use_germlines: bool = False):
    if not use_germlines:
        self.endpoint = endpoint
    else:
        from sadie.germlines import get_manager
        from sadie.germlines.g3_adapter import GermlineToG3Adapter
        self.germline_manager = get_manager()
        self.g3_adapter = GermlineToG3Adapter()
```

### 2. Gap: References.from_yaml() Doesn't Propagate Flag

**Current code (reference.py:313-339):**
```python
@staticmethod
def from_yaml(yaml_path: Optional[Path] = None) -> "References":
    yaml_ref_object = YamlRef(yaml_path)
    references_object = References()
    for name in yaml_ref:
        reference_object = Reference()  # ← No use_germlines parameter!
        for source in yaml_ref.get(name):
            for species in yaml_ref.get(name).get(source):
                reference_object.add_genes(species, source, list_of_genes)
        references_object.add_reference(name, reference_object)
    return references_object
```

**Required change:**
```python
@staticmethod
def from_yaml(yaml_path: Optional[Path] = None, use_germlines: bool = False) -> "References":
    # ...
    for name in yaml_ref:
        reference_object = Reference(use_germlines=use_germlines)  # Pass flag
```

### 3. Gap: Source Validation Restricts to imgt/custom

**Current code (models.py:28-31):**
```python
@field_validator("source")
@classmethod
def check_source(cls, v: str) -> str:
    if v not in ["imgt", "custom"]:
        raise ValueError(f"{v} is not a valid source")
```

**Required change** (for multi-source support):
```python
VALID_SOURCES = {"imgt", "ogrdb", "vdjbase", "custom"}

@field_validator("source")
@classmethod
def check_source(cls, v: str) -> str:
    if v not in VALID_SOURCES:
        raise ValueError(f"{v} is not a valid source, choices: {VALID_SOURCES}")
```

### 4. Gap: G3 Adapter Missing Fields

The `GermlineToG3Adapter` produces output missing some fields expected by downstream consumers:

**Currently produced:**
- `source`, `common`, `latin`, `gene`, `label`, `gene_segment`, `receptor`, `sequence`, `species`
- `imgt.*` fields (partial)

**Potentially missing:**
- `_id` - MongoDB-style ID used by from_dataframe()
- `gene_curation_source` - Curation metadata
- `chimera` - Chimeric gene flag

---

## Data Flow

### Current Flow (G3 API)

```
┌─────────────┐     ┌─────────────┐     ┌──────────────┐     ┌──────────┐
│ reference.yml│ ──▶ │   YamlRef   │ ──▶ │ References   │ ──▶ │  G3 API  │
│ (YAML file) │     │  (parser)   │     │ .from_yaml() │     │ (remote) │
└─────────────┘     └─────────────┘     └──────────────┘     └──────────┘
                                               │
                                               ▼
                                        ┌──────────────┐
                                        │  Reference   │
                                        │ .add_genes() │
                                        └──────────────┘
                                               │
                                               ▼
                                        ┌──────────────┐
                                        │  _get_genes  │ ──▶ HTTP GET to G3
                                        └──────────────┘
```

### Target Flow (Germlines Module)

```
┌─────────────┐     ┌─────────────┐     ┌──────────────────────┐
│ reference.yml│ ──▶ │   YamlRef   │ ──▶ │     References       │
│ (YAML file) │     │  (parser)   │     │ .from_yaml(          │
└─────────────┘     └─────────────┘     │   use_germlines=True)│
                                        └──────────────────────┘
                                               │
                                               ▼
                                        ┌──────────────────────┐
                                        │ Reference(           │
                                        │   use_germlines=True)│
                                        └──────────────────────┘
                                               │
                                               ▼
                                        ┌──────────────────────┐
                                        │ _get_genes()         │
                                        │   if use_germlines:  │
                                        │     GermlineManager  │──▶ Local Data
                                        │       .get_gene_by_  │
                                        │        name()        │
                                        │     GermlineToG3     │
                                        │       Adapter        │
                                        └──────────────────────┘
```

**Key change:** The `_get_genes()` method already branches on `self.use_germlines`. We just need to ensure the flag reaches it through `from_yaml()`.

---

## New Components

### No New Components Required

The existing architecture is sufficient. We need only:

1. **Parameter propagation** - Add `use_germlines` param to `from_yaml()`
2. **Validation expansion** - Update `VALID_SOURCES` in models.py
3. **Adapter enhancement** - Add missing fields to G3 adapter output

### Optional Enhancement: Provider Selection

If fine-grained provider control is needed, consider:

```python
def from_yaml(
    yaml_path: Optional[Path] = None,
    use_germlines: bool = False,
    providers: Optional[List[str]] = None  # ["imgt", "ogrdb"]
) -> "References":
```

This would pass `providers` through to `GermlineManager(providers=providers)`.

---

## Build Order

### Phase 1: Core Integration (Critical Path)

| Step | File | Change | Effort |
|------|------|--------|--------|
| 1.1 | `models.py` | Expand `VALID_SOURCES` to include ogrdb, vdjbase | Small |
| 1.2 | `reference.py` | Add `use_germlines` param to `from_yaml()` | Small |
| 1.3 | `reference.py` | Propagate flag to `Reference()` instantiation | Small |
| 1.4 | Tests | Test `from_yaml(use_germlines=True)` | Medium |

**Dependencies:** 1.1 before 1.2 (source validation must pass first)

### Phase 2: Adapter Completeness

| Step | File | Change | Effort |
|------|------|--------|--------|
| 2.1 | `g3_adapter.py` | Add `_id` generation (UUID or hash) | Small |
| 2.2 | `g3_adapter.py` | Add `gene_curation_source` field | Small |
| 2.3 | `g3_adapter.py` | Add `chimera` field (default False) | Small |
| 2.4 | Tests | Verify adapter output matches G3 schema | Medium |

**Dependencies:** Can run in parallel with Phase 1

### Phase 3: CLI & API Surface

| Step | File | Change | Effort |
|------|------|--------|--------|
| 3.1 | CLI (if exists) | Add `--use-germlines` flag | Small |
| 3.2 | Documentation | Update usage examples | Small |

**Dependencies:** After Phase 1

### Dependency Graph

```
           ┌────────┐
           │ 1.1    │ (models validation)
           └────┬───┘
                │
           ┌────▼───┐     ┌────────┐
           │ 1.2    │     │ 2.1    │
           │ 1.3    │     │ 2.2    │ (adapter, parallel)
           └────┬───┘     │ 2.3    │
                │         └────┬───┘
           ┌────▼───┐          │
           │ 1.4    │◄─────────┤
           │ Tests  │     ┌────▼───┐
           └────┬───┘     │ 2.4    │
                │         └────┬───┘
           ┌────▼────────────▼───┐
           │      3.1, 3.2       │
           └─────────────────────┘
```

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| G3 adapter output incompatible with downstream | Medium | High | Add comprehensive schema tests |
| Missing genes in germlines vs G3 | Low | Medium | Fall back to G3 if gene not found locally |
| Performance regression (local slower than API) | Low | Low | Local should be faster; benchmark to confirm |

---

## Recommendations

1. **Start with Phase 1** - The core integration is straightforward and unblocks offline usage
2. **Add fallback to G3** - For genes missing from local germlines, optionally fall back to G3 API
3. **Schema validation tests** - Add tests comparing adapter output against known-good G3 responses
4. **Consider caching** - If G3 fallback is used, cache responses to local germlines format

---

## Appendix: Key Code References

### Reference._get_genes() with germlines (already exists)
```python
# reference.py:184-203
def _get_genes(self, genes: GeneEntries) -> List[Dict[str, str]]:
    if self.use_germlines:
        from sadie.germlines import get_gene_by_name
        results = []
        for gene_name in genes.genes:
            germline_gene = get_gene_by_name(gene_name, genes.species)
            if germline_gene:
                g3_dict = self.g3_adapter.to_g3_format(germline_gene)
                g3_dict["species"] = genes.species
                results.append(g3_dict)
        return results
    # ... G3 API path
```

### GermlineToG3Adapter.to_g3_format() signature
```python
# g3_adapter.py:46
def to_g3_format(self, gene: GermlineGene) -> Dict[str, Any]:
    """Convert GermlineGene to G3 API response format."""
```

### GermlineManager.get_gene_by_name() signature
```python
# manager.py:175
def get_gene_by_name(self, name: str, species: str) -> Optional[GermlineGene]:
    """Get specific gene by name (first provider wins)."""
```
