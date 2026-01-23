# G3 Backend Technology Stack

> Analysis focused on J gene and constant region segment position discovery in the G3 (non-germlines module) backend for AIRR annotation.

## Core Technologies

### Python Version
- **Minimum**: Python 3.10
- **Supported**: Python 3.10, 3.11, 3.12, 3.13
- **Defined in**: `pyproject.toml`

### Primary Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `biopython` | >=1.80 | Sequence parsing (SeqIO), SeqRecord handling |
| `pandas` | >=1.5 | Data manipulation, AIRR table handling |
| `pydantic` | >=2.0.0,<3.0.0 | Data validation, models |
| `semantic-version` | ^2.10.0 | IgBLAST version parsing |
| `requests` | ^2.32.0 | G3 API HTTP calls |
| `PyYAML` | ^6.0 | Reference YAML file parsing |

### IgBLAST Integration

- **External binary**: `igblastn` (required at runtime)
- **Version detection**: Automatic via semantic versioning
- **Platform binaries**: Pre-compiled for macOS/Linux in `src/sadie/airr/bin/`
- **Python wrapper**: `src/sadie/airr/igblast/igblast.py` → `IgBLASTN` class

## Feature Flag System

### Environment Variable
```bash
SADIE_USE_GERMLINES_MODULE=true|false
```

- **Default**: `true` (uses new germlines module)
- **When `false`**: Uses legacy G3 paths (deprecated, removal date: 2026-06-01)
- **Implementation**: `src/sadie/airr/igblast/germline.py` → `_use_germlines_module()`

## Data File Types

### 1. Auxiliary Files (`.aux`)

**Location (G3 legacy)**: `src/sadie/airr/data/germlines/aux_db/imgt/{species}_gl.aux`

**Format for J genes** (5 columns, tab-separated):
```
<gene_name>	<reading_frame>	<chain_type>	<cdr3_end>	<is_functional>
```

**Example**:
```
IGHJ1*01	0	JH	17	1
IGHJ2*01	1	JH	18	1
IGKJ1*01	1	JK	6	1
```

**Purpose**: CDR3 boundary calculation - tells IgBLAST where CDR3 ends relative to J gene alignment start.

### 2. NDM Files (`.ndm.imgt`)

**Location (G3 legacy)**: `src/sadie/airr/data/germlines/Ig/internal_data/{species}/{species}.ndm.imgt`

**Format for V genes** (12 columns, tab-separated):
```
<gene>	<FWR1_start>	<CDR1_start>	<CDR1_end>	<FWR2_start>	<FWR2_end>	<CDR2_start>	<CDR2_end>	<FWR3_start>	<FWR3_end>	<chain>	<flags>
```

**Example**:
```
IGHV1-18*01	1	75	76	99	100	150	151	174	175	288	VH	0
```

**Purpose**: V gene region boundaries (FWR1, CDR1, FWR2, CDR2, FWR3 positions).

### 3. BLAST Databases

**Location (G3 legacy)**: `src/sadie/airr/data/germlines/Ig/blastdb/{species}/`

**File patterns**:
- `{species}_V.*` - V gene database files
- `{species}_D.*` - D gene database files
- `{species}_J.*` - J gene database files
- `{species}_C.*` - C gene database files

**Extensions**: `.ndb`, `.nhi`, `.nhr`, `.nin`, `.nog`, `.nos`, `.not`, `.nsq`, `.ntf`, `.nto`, `.fasta`

### 4. Reference YAML

**Location**: `src/sadie/reference/data/reference.yml`

**Purpose**: Maps species names to gene lists for G3 API queries.

**Structure**:
```yaml
{species_name}:
  {source}:
    {species}:
      - IGHV1-2*02
      - IGKV1-33*01
      ...
```

## Key Classes and Files

### GermlineData Class
**File**: `src/sadie/airr/igblast/germline.py`

**Responsibility**: Manages paths to all germline data files for a species.

**Key properties**:
- `base_dir` - Base germline data directory
- `v_gene_dir` - V gene BLAST database prefix
- `d_gene_dir` - D gene BLAST database prefix
- `j_gene_dir` - J gene BLAST database prefix
- `c_gene_dir` - C gene BLAST database prefix
- `aux_path` - Auxiliary file path for CDR3 boundaries
- `igdata` - IGDATA environment variable path

### IgBLASTN Class
**File**: `src/sadie/airr/igblast/igblast.py`

**Responsibility**: Python wrapper for IgBLAST command-line execution.

**Key parameters for segment discovery**:
- `germline_db_V` - V gene database path
- `germline_db_D` - D gene database path
- `germline_db_J` - J gene database path
- `c_region_db` - C region database path
- `auxiliary_data` - Aux file for J gene CDR3 boundaries

### Airr Class
**File**: `src/sadie/airr/airr.py`

**Responsibility**: High-level API for AIRR annotation.

**Key initialization**:
1. Sets up GermlineData for species
2. Configures IgBLASTN with database paths
3. Manages penalty parameters for alignment
4. Handles adaptive penalty correction for failed annotations

## G3 vs Germlines Module Path Resolution

### G3 Legacy Paths
```python
base_dir = Path(__file__).parent / "../data/germlines/"
blast_dir = f"{base_dir}/{receptor}/blastdb/{name}/{name}_"
aux_path = base_dir / f"aux_db/{scheme}/{name}_gl.aux"
igdata = base_dir / f"{receptor}/"
```

### Germlines Module Paths
```python
germlines_igblast = get_germlines_base_dir() / "igblast"
blast_dir = germlines_igblast / "Ig" / "internal_data" / name / f"{name}_"
aux_path = germlines_igblast / "aux_db" / f"{name}_gl.aux"
igdata = germlines_igblast / "Ig"
```

## Species Support (G3 Legacy)

Built-in species in `src/sadie/airr/data/germlines/`:
- `human`
- `mouse`
- `rabbit`
- `rat`
- `dog`
- `macaque`
- `clk` (Custom library)
- `se09` (Custom library)

## Build/Test Tools

| Tool | Purpose |
|------|---------|
| `poetry` | Dependency management |
| `pytest` | Test framework |
| `mypy` | Static type checking |
| `pyright` | Additional type checking |
| `black` | Code formatting |
| `flake8` | Linting |
| `pre-commit` | Git hooks |

## Data Type Constants

**File**: `src/sadie/airr/airrtable/constants.py`

Defines pandas dtypes for:
- `IGBLAST_AIRR` - Core IgBLAST output columns (V, D, J gene data)
- `CONSTANTS_AIRR` - C region columns
- `OTHER_COLS` - Additional SADIE-specific columns
