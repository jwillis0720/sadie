# Stack Research: Reference Module Unification

## Current Stack

### Reference Module (`src/sadie/reference/`)
| File | Purpose | Integration Status |
|------|---------|-------------------|
| `models.py` | Pydantic validators for GeneEntry/GeneEntries | ⚠️ Hardcoded to `["imgt", "custom"]` only |
| `reference.py` | Main Reference/References classes | ✅ Has `use_germlines` param with full integration |
| `yaml.py` | YamlRef parser for reference.yml | ❌ Does NOT pass `use_germlines` to Reference |
| `reference.yml` | Allele selection configuration | ⚠️ Schema only supports imgt/custom sources |

### Germlines Module (`src/sadie/germlines/`)
| Component | Purpose | Ready for Integration |
|-----------|---------|----------------------|
| `GermlineManager` | Priority-based provider lookup | ✅ Fully functional |
| `GermlineGene` | Unified gene model | ✅ Has all required fields |
| `GermlineToG3Adapter` | Converts GermlineGene → G3 dict | ✅ Already used in reference.py |
| 4 Providers | imgt, ogrdb, vdjbase, custom | ✅ All operational |
| Public API | `get_manager()`, `get_gene_by_name()` | ✅ Exposed in `__init__.py` |

### Existing Integration (Already Implemented)
```python
# reference.py lines 31-38
class Reference:
    def __init__(self, endpoint: str = _endpoint, use_germlines: bool = False):
        if use_germlines:
            from sadie.germlines import get_manager
            from sadie.germlines.g3_adapter import GermlineToG3Adapter
            self.germline_manager = get_manager()
            self.g3_adapter = GermlineToG3Adapter()
```

The `_get_gene()` and `_get_genes()` methods already have complete `if self.use_germlines:` branches that work correctly.

## Required Changes

### 1. Expand Source Validation (models.py)
**File**: `src/sadie/reference/models.py`

**Current** (lines 24-26, 44-46):
```python
@field_validator("source")
def check_source(cls, v: str) -> str:
    if v not in ["imgt", "custom"]:
        raise ValueError(f"{v} is not a valid source, choices are 'imgt' or 'custom'")
    return v
```

**Change to**:
```python
VALID_SOURCES = ["imgt", "ogrdb", "vdjbase", "custom"]

@field_validator("source")
def check_source(cls, v: str) -> str:
    if v not in VALID_SOURCES:
        raise ValueError(f"{v} is not a valid source, choices are {VALID_SOURCES}")
    return v
```

**Rationale**: This is the primary blocker - reference.yml cannot use ogrdb/vdjbase sources because validation fails before data is even fetched.

### 2. Enable Germlines Mode in References.from_yaml() (reference.py)
**File**: `src/sadie/reference/reference.py`

**Current** (line ~320):
```python
@staticmethod
def from_yaml(yaml_path: Optional[Path] = None) -> "References":
    # ...
    for name in yaml_ref:
        reference_object = Reference()  # Always uses G3 API
```

**Change to**:
```python
@staticmethod
def from_yaml(yaml_path: Optional[Path] = None, use_germlines: bool = True) -> "References":
    # ...
    for name in yaml_ref:
        reference_object = Reference(use_germlines=use_germlines)
```

**Rationale**: The Reference class already supports germlines integration, but `from_yaml()` never activates it. Default to `True` for new behavior while allowing `False` for backwards compatibility.

### 3. NO Changes Needed to yaml.py
The YamlRef class is data-agnostic - it simply parses YAML structure. No modifications required.

### 4. NO Changes Needed to g3_adapter.py
The adapter already correctly transforms GermlineGene → G3 dict format. Testing has validated this works.

### 5. NO Changes Needed to Germlines Module
All infrastructure is already in place and tested.

## Integration Approach

### Architecture Pattern: Adapter + Strategy
```
                    ┌─────────────────────────────────────┐
                    │         reference.yml               │
                    │   (now supports all 4 sources)      │
                    └─────────────────┬───────────────────┘
                                      │
                                      ▼
                    ┌─────────────────────────────────────┐
                    │            YamlRef                   │
                    │      (unchanged - data parser)       │
                    └─────────────────┬───────────────────┘
                                      │
                                      ▼
┌───────────────────────────────────────────────────────────────────────┐
│                    References.from_yaml()                             │
│                 use_germlines=True (new default)                      │
└───────────────────────────────┬───────────────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────────────┐
│                       Reference class                                 │
│                                                                       │
│   if use_germlines:              │    else (backwards compat):        │
│     GermlineManager              │      G3 REST API                   │
│        │                         │         │                          │
│        ▼                         │         ▼                          │
│   GermlineGene                   │    JSON response                   │
│        │                         │         │                          │
│        ▼                         │         │                          │
│   GermlineToG3Adapter            │         │                          │
│        │                         │         │                          │
│        └─────────────────────────┴─────────┘                          │
│                         │                                             │
│                         ▼                                             │
│              Unified G3-format dict                                   │
│              (same downstream processing)                             │
└───────────────────────────────────────────────────────────────────────┘
```

### Data Flow
1. User specifies genes in reference.yml with source: ogrdb/vdjbase/imgt/custom
2. `References.from_yaml()` parses YAML, creates Reference instances with `use_germlines=True`
3. `Reference.add_genes()` calls `_get_genes()` which queries GermlineManager
4. GermlineManager delegates to appropriate provider based on source
5. `GermlineToG3Adapter` transforms result to G3 dict format
6. Downstream IgBLAST database building proceeds unchanged

## Backwards Compatibility

### What MUST NOT Change
| Component | Reason |
|-----------|--------|
| G3 dict schema | All downstream code (IgBLAST building) expects this format |
| `Reference.get_dataframe()` | Used by multiple pipelines |
| `Reference.from_dataframe()` | Serialization/deserialization paths |
| YAML structure | Existing reference.yml files must continue to work |
| CLI interfaces | External scripts depend on current behavior |

### What CAN Change (with migration path)
| Change | Migration |
|--------|-----------|
| `from_yaml()` signature | Add optional `use_germlines=True` param with backwards-compatible default |
| Source validation | Expand allowed values (additive change) |
| Default behavior | New installs get germlines; explicit `use_germlines=False` preserves G3 API |

### Breaking Change Risk Assessment
- **NONE**: All changes are additive or behind feature flags
- The G3 API path remains fully functional
- Existing tests will pass without modification

## Recommendations

### Prescriptive Guidance

1. **Start with models.py** - This unblocks everything else with a 2-line change

2. **Add parameter to from_yaml()** - Make `use_germlines=True` the default, allowing:
   ```python
   # New behavior (default)
   refs = References.from_yaml()  # Uses germlines
   
   # Legacy behavior (explicit)
   refs = References.from_yaml(use_germlines=False)  # Uses G3 API
   ```

3. **Update reference.yml examples** - Show ogrdb/vdjbase sources in documentation

4. **DO NOT**:
   - Add new external dependencies
   - Change the G3 dict schema
   - Remove G3 API support
   - Modify GermlineManager or provider logic

### Implementation Order
```
1. models.py:VALID_SOURCES expansion     → 5 min
2. reference.py:from_yaml() parameter    → 10 min  
3. Test with existing reference.yml      → 5 min
4. Add ogrdb/vdjbase entries to YAML     → 15 min
5. Integration tests                     → 30 min
```

### No New External Dependencies Needed
All required libraries are already in `pyproject.toml`:
- `pydantic` (validation)
- `requests` (G3 API fallback)
- `pyyaml` (YAML parsing)
- `biopython` (sequence handling)

The germlines module is fully self-contained within the project.
