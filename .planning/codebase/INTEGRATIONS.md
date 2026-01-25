# SADIE External Integrations

## IgBLAST Integration

### Overview
IgBLAST (Immunoglobulin BLAST) is the core annotation engine for AIRR analysis.

### Key Files
- `src/sadie/airr/airr.py` → `Airr` class wrapping IgBLAST
- `src/sadie/airr/igblast/igblast.py` → `IgBLASTN` execution wrapper
- `src/sadie/airr/igblast/germline.py` → `GermlineData` path management

### Integration Flow
```
Input FASTA → Airr.run_fasta() → IgBLASTN.run_file() → Parse Output → AirrTable
```

### Database Parameters
```python
Airr(
    reference_name="human",      # Species/reference name
    database=Path("/path/to/db") # Prebuilt database (v1.2 feature)
)
```

### Prebuilt Database Support (v1.2)
The `Airr(database=<path>)` parameter enables using prebuilt databases:
```python
# From sadie reference build output
airr = Airr("custom_ref", database="/data/my_germlines")
```

Expected structure validated by `validate_prebuilt_database()`:
```
database_path/
├── Ig/
│   ├── blastdb/{name}/{name}_V, _D, _J
│   └── internal_data/{name}/{name}.ndm.imgt
├── aux_db/imgt/{name}_gl.aux
└── .references_dataframe.csv.gz (optional)
```

---

## Germlines Module Data Sources

### Provider Architecture (v1.2)
Located in `src/sadie/germlines/providers/`:

| Provider | File | Data Source |
|----------|------|-------------|
| IMGT | `imgt.py` | IMGT/GENE-DB reference sequences |
| OGRDB | `ogrdb.py` | Open Germline Receptor Database |
| VDJbase | `vdjbase.py` | Population-specific genotypes |
| Custom | `custom.py` | User-provided FASTA files |

### Priority System
Default order (first wins on conflicts):
```python
GermlineManager.DEFAULT_PROVIDERS = ["custom", "ogrdb", "vdjbase", "imgt"]
```

### GermlineManager API
```python
from sadie.germlines import GermlineManager

# Default priority
manager = GermlineManager()

# Custom priority - explicit source
manager = GermlineManager(providers=["imgt"])

# Fetch genes
genes = manager.get_genes("human", "V", "H", functional_only=True)

# Get specific gene
gene = manager.get_gene_by_name("IGHV1-69*01", "human")
```

### Data Directory Structure
```
src/sadie/germlines/
├── sources/
│   ├── imgt/{species}/*.fasta
│   ├── ogrdb/{species}/*.fasta
│   ├── vdjbase/{species}/*.fasta
│   └── custom/{species}/*.fasta
├── normalized/{species}/
└── igblast/
    ├── Ig/internal_data/{species}/
    └── aux_db/{species}_gl.aux
```

---

## Reference Module Integration (v1.2 Feature)

### Core Classes
**File**: `src/sadie/reference/reference.py`

```python
class Reference:
    def __init__(self, endpoint: str = _endpoint, use_germlines: bool = False):
        """
        endpoint: G3 API URL (legacy, ignored if use_germlines=True)
        use_germlines: Use local germlines module instead of G3 API
        """
```

### use_germlines Parameter
When `use_germlines=True`:
1. Imports `GermlineManager` and `GermlineToG3Adapter`
2. Creates manager with explicit source (no fallback)
3. Fetches gene via `manager.get_gene_by_name()`
4. Transforms to G3 format via adapter

```python
# In Reference._get_gene():
if self.use_germlines:
    from sadie.germlines import GermlineManager
    from sadie.germlines.g3_adapter import GermlineToG3Adapter
    
    manager = GermlineManager(providers=[gene.source])
    germline_gene = manager.get_gene_by_name(gene.gene, gene.species)
    adapter = GermlineToG3Adapter()
    return adapter.to_g3_format(germline_gene)
```

### Valid Sources
**File**: `src/sadie/reference/models.py`
```python
VALID_SOURCES = ["imgt", "ogrdb", "vdjbase", "custom"]
```

Pydantic validation ensures only valid sources are accepted:
```python
@field_validator("source")
def check_source(cls, v: str) -> str:
    if v not in VALID_SOURCES:
        raise ValueError(f"{v} is not a valid source")
    return v
```

---

## G3 Adapter Integration

### Purpose
Transforms `GermlineGene` objects to G3 API response format for backward compatibility.

**File**: `src/sadie/germlines/g3_adapter.py`

### Transformation
```python
class GermlineToG3Adapter:
    def to_g3_format(self, gene: GermlineGene) -> Dict[str, Any]:
        """
        Input: GermlineGene from germlines module
        Output: Dict matching G3 API JSON structure
        """
        return {
            "_id": self._generate_id(gene.source, gene.species, gene.name),
            "source": gene.source,
            "common": gene.species,
            "gene": gene.name,
            "sequence": gene.sequence,
            "imgt": {
                "imgt_functional": gene.functionality,
                "sequence_gapped": gene.sequence_gapped,
                # ... CDR/FWR regions
            }
        }
```

---

## G3 API Integration (Legacy)

### Status
**DEPRECATED** - Being replaced by germlines module. Deprecation warning added.

### API Endpoint
```python
Reference._endpoint = "https://g3.jordanrwillis.com/api/v1/genes"
```

### Environment Control
```bash
# Force legacy G3 API (not recommended)
export SADIE_USE_GERMLINES_MODULE=false
```

### Deprecation Timeline
```python
# In germline.py
logger.warning(
    "G3 API is deprecated. Set SADIE_USE_GERMLINES_MODULE=true. "
    "G3 will be removed after 2026-06-01."
)
```

---

## CLI Commands and Integrations

### Entry Point
**File**: `src/sadie/app.py`
```python
@click.group()
def sadie(): pass
```

### Commands

#### `sadie airr`
Runs AIRR annotation pipeline.
```bash
sadie airr --name human input.fasta output.tsv.gz
```
Integrations: IgBLAST, GermlineData

#### `sadie renumbering`
Runs antibody numbering/renumbering.
```bash
sadie renumbering -q input.fasta -s imgt -o output.csv
```
Integrations: HMMER, pyhmmer

#### `sadie reference make` (Legacy)
Builds reference from G3 API.
```bash
sadie reference make --outpath ./germlines --reference reference.yml
```
Integrations: G3 API, BLAST makeblastdb

#### `sadie reference build` (v1.2)
Builds IgBLAST database with germlines support.
```bash
sadie reference build reference.yml --output ./db --use-germlines
```
Integrations: GermlineManager, References.from_yaml()

#### `sadie germlines populate`
Downloads germline data from providers.
```bash
sadie germlines populate -p imgt -s human
```
Integrations: IMGT/OGRDB/VDJbase providers

#### `sadie germlines status`
Shows germline database status.
```bash
sadie germlines status
```
Integrations: Provider metadata

#### `sadie make-all`
Comprehensive database generation.
```bash
sadie make-all --outpath ./data --species human
```
Integrations: References, IgBLAST, auxiliary files

---

## Renumbering Integration

### HMMER Integration
**File**: `src/sadie/renumbering/renumbering.py`

Uses pyhmmer for antibody chain numbering:
```python
class Renumbering:
    def __init__(
        self,
        scheme="imgt",           # Numbering scheme
        region_assign="imgt",    # CDR/FWR definition
        allowed_species=["human"],
        use_numbering_hmms=True  # Use germlines HMMs
    ):
        self.hmmer = HMMER(...)
```

### Germlines Integration (v1.2)
**File**: `src/sadie/germlines/renumbering_integration.py`

Provides HMM paths from germlines module for renumbering.

---

## Integration Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                         CLI (app.py)                            │
├─────────────────────────────────────────────────────────────────┤
│  sadie airr  │  sadie reference build  │  sadie germlines      │
└──────┬───────┴──────────┬──────────────┴──────────┬────────────┘
       │                  │                         │
       ▼                  ▼                         ▼
┌─────────────┐   ┌───────────────┐        ┌───────────────────┐
│  Airr Class │   │  References   │        │  GermlineManager  │
│  (airr.py)  │   │ (reference.py)│        │   (manager.py)    │
└──────┬──────┘   └───────┬───────┘        └─────────┬─────────┘
       │                  │                          │
       │         ┌────────┴────────┐                 │
       │         │ use_germlines?  │                 │
       │         └────────┬────────┘                 │
       │                  │                          │
       │    ┌─────────────┼─────────────┐            │
       │    │ True        │       False │            │
       │    ▼             │             ▼            │
       │ ┌──────────┐     │      ┌──────────┐        │
       │ │G3 Adapter│     │      │ G3 API   │        │
       │ └────┬─────┘     │      │(deprecated)       │
       │      │           │      └──────────┘        │
       │      └───────────┼──────────────────────────┘
       │                  │
       ▼                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                    GermlineData (germline.py)                   │
├─────────────────────────────────────────────────────────────────┤
│  - Validates prebuilt database structure                        │
│  - Provides paths to IgBLAST databases                          │
│  - Supports germlines/ or legacy G3 paths                       │
└─────────────────────────────────────────────────────────────────┘
       │
       ▼
┌─────────────────────────────────────────────────────────────────┐
│                      IgBLAST Execution                          │
├─────────────────────────────────────────────────────────────────┤
│  Input: FASTA sequences                                         │
│  Databases: V/D/J BLAST DBs + auxiliary files + internal data   │
│  Output: AIRR-formatted annotation tables                       │
└─────────────────────────────────────────────────────────────────┘
```

---

## Data Flow Summary

### Germlines → Reference → AIRR Pipeline (v1.2)
1. **User specifies**: `sadie reference build ref.yml --use-germlines`
2. **References.from_yaml()**: Parses YAML with `use_germlines=True`
3. **Reference._get_gene()**: Uses `GermlineManager` with explicit provider
4. **GermlineToG3Adapter**: Transforms to G3 format for compatibility
5. **References.make_airr_database()**: Builds IgBLAST databases
6. **Output**: Complete database structure at `--output` path
7. **Airr(database=path)**: Uses prebuilt database directly
