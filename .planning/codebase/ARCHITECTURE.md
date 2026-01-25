# SADIE Architecture

## Overview

SADIE (Sequence Analysis and Design for Immunology Engineering) is an immunoglobulin sequence annotation and analysis toolkit built on IgBLAST. The architecture centers on the **Reference Module (v1.2)** which unifies germline data management from multiple sources.

## Core Module Organization

```
src/sadie/
├── airr/          # AIRR-standard annotation (main entry point)
├── germlines/     # Multi-source germline database management (v1.2)
├── reference/     # YAML configuration → IgBLAST database builder
├── renumbering/   # Antibody numbering (IMGT, Kabat, etc.)
├── numbering/     # Legacy numbering (being deprecated)
├── receptor/      # Receptor utilities
├── cluster/       # Sequence clustering
├── typing/        # Type definitions (Species, Chain, Source)
└── utility/       # Shared utilities
```

## Reference Module Architecture (v1.2)

### Data Flow Pipeline

```
┌──────────────────────────────────────────────────────────────────────┐
│                     YAML Configuration (reference.yml)               │
│   name:                                                              │
│     source: [imgt, ogrdb, vdjbase, custom]                          │
│       species:                                                       │
│         - IGHV1-69*01                                                │
│         - IGHD3-3*01                                                 │
└──────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────┐
│                     References.from_yaml()                           │
│   src/sadie/reference/reference.py                                   │
│   - Parses YAML via YamlRef                                         │
│   - Creates Reference objects per name                               │
│   - Fetches genes via GermlineManager (or legacy G3 API)            │
└──────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────┐
│                     GermlineManager                                  │
│   src/sadie/germlines/manager.py                                    │
│   - Priority-based database lookup (custom > imgt > ogrdb > vdjbase)│
│   - Deduplication: first provider wins on name/sequence conflicts   │
│   - Returns GermlineGene objects                                     │
└──────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────┐
│                     GermlineToG3Adapter                              │
│   src/sadie/germlines/g3_adapter.py                                 │
│   - Transforms GermlineGene → G3 API format                         │
│   - Enables backward compatibility with existing Reference code      │
│   - Generates deterministic _id via SHA-256                         │
└──────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────┐
│                     References.make_airr_database()                  │
│   - _make_internal_annotation_file() → internal_data/{name}/*.ndm   │
│   - _make_igblast_ref_database() → blastdb/{name}/{name}_V,D,J      │
│   - _make_auxillary_file() → aux_db/{scheme}/{name}_gl.aux          │
└──────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────┐
│                     IgBLAST Database Structure                       │
│   output_path/                                                       │
│   ├── Ig/                                                           │
│   │   ├── blastdb/{name}/{name}_V, {name}_D, {name}_J              │
│   │   └── internal_data/{name}/{name}.ndm.imgt                     │
│   ├── aux_db/imgt/{name}_gl.aux                                    │
│   └── .references_dataframe.csv.gz                                  │
└──────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────┐
│                     Airr(database=<path>) or Airr(reference_name)    │
│   src/sadie/airr/airr.py                                            │
│   - Uses GermlineData to locate databases                           │
│   - Executes IgBLASTN with configured paths                         │
│   - Returns AirrTable with AIRR-standard annotations                │
└──────────────────────────────────────────────────────────────────────┘
```

### Provider System

The germlines module implements a pluggable provider architecture:

```
GermlineProvider (abstract base)
├── IMGTProvider      # IMGT reference database
├── OGRDBProvider     # OGRDB novel alleles
├── VDJbaseProvider   # VDJbase population-specific
└── CustomProvider    # User-defined sequences
```

**Provider Priority** (default): `["custom", "ogrdb", "vdjbase", "imgt"]`

**Deduplication Rules**:
1. Same gene name → first provider wins
2. Same exact sequence → first provider wins  
3. Novel gene → include from any provider

### Integration Points

#### AIRR ↔ Germlines
```python
# airr/igblast/germline.py
class GermlineData:
    def __init__(self, name, receptor, database_dir=None, scheme="imgt", prebuilt=False):
        # Option 1: Use prebuilt database (database= parameter)
        if prebuilt:
            paths = validate_prebuilt_database(database_dir, name)
            
        # Option 2: Use germlines module (SADIE_USE_GERMLINES_MODULE=true)
        elif _use_germlines_module():
            germlines_igblast = _get_germlines_igblast_dir()
            
        # Option 3: Legacy G3 API paths (deprecated)
        else:
            self._use_legacy_paths(name, receptor, scheme)
```

#### Reference ↔ Germlines (via Adapter)
```python
# reference/reference.py
class Reference:
    def _get_gene(self, gene: GeneEntry):
        if self.use_germlines:
            manager = GermlineManager(providers=[gene.source])
            germline_gene = manager.get_gene_by_name(gene.gene, gene.species)
            
            # Transform to G3 format for backward compatibility
            adapter = GermlineToG3Adapter()
            return adapter.to_g3_format(germline_gene)
```

#### Renumbering ↔ Germlines (HMM Building)
```python
# germlines/renumbering_integration.py
class LocalHMMBuilder:
    def get_hmm(self, species, chain, source="imgt"):
        # Build HMM from local germlines instead of G3 API
        vj_pairs = self._get_vj_alignment_pairs(species, chain, source)
        # Use pyhmmer for HMM construction
```

## Key Design Patterns

### 1. Adapter Pattern (G3 Compatibility)
`GermlineToG3Adapter` transforms GermlineGene objects into the legacy G3 API response format, enabling incremental migration without breaking existing code.

### 2. Strategy Pattern (Providers)
Each `GermlineProvider` subclass implements the same interface but retrieves data from different sources (IMGT, OGRDB, VDJbase, custom files).

### 3. Priority Chain
`GermlineManager` iterates providers in priority order, implementing a chain-of-responsibility for gene lookup with deduplication.

### 4. Feature Flag Pattern
```python
SADIE_USE_GERMLINES_MODULE=true   # Use local germlines (default)
SADIE_USE_GERMLINES_MODULE=false  # Use deprecated G3 API
```

### 5. Prebuilt Database Support (v1.2)
```python
# Build once
references = References.from_yaml("reference.yml", use_germlines=True)
references.make_airr_database("/path/to/database")

# Use anywhere
airr = Airr("human", database="/path/to/database")
```

## Module Dependencies

```
airr
├── reference.References (for custom databases)
├── germlines.GermlineManager (gene lookup)
└── igblast (blast database execution)

reference
├── germlines.GermlineManager (gene fetching)
├── germlines.g3_adapter.GermlineToG3Adapter (format conversion)
└── yaml (configuration parsing)

germlines
├── providers/ (IMGT, OGRDB, VDJbase, custom)
├── pipeline (normalize → build workflow)
└── renumbering_integration (HMM building)

renumbering
└── germlines.renumbering_integration.LocalHMMBuilder
```

## Validation Models

**Pydantic Models** in `reference/models.py`:
- `GeneEntry`: Single gene with species/gene/source validation
- `GeneEntries`: Multiple genes with consistent validation
- `VALID_SOURCES`: `["imgt", "ogrdb", "vdjbase", "custom"]`

**Germline Models** in `germlines/models.py`:
- `GermlineGene`: Core gene representation with regions and positions
- `ProviderMetadata`: Provider version and species availability

## Error Handling Strategy

1. **Strict Mode**: Raise errors when required data is missing (FR-014)
2. **Graceful Degradation**: Warn and continue when optional data is missing
3. **Fail-Fast**: Validate database structure upfront (`validate_prebuilt_database`)
4. **Deprecation Warnings**: G3 API usage triggers warnings with migration guidance
