# SADIE Technology Stack

## Python Version
- **Required**: Python >= 3.10
- **Supported**: 3.10, 3.11, 3.12, 3.13
- **Implementations**: CPython, PyPy

## Package Management
- **Build System**: Poetry with poetry-dynamic-versioning
- **Version Source**: Git tags (pattern: `^(?:test-)?v(?P<base>\d+\.\d+\.\d+)`)
- **Lock File**: `poetry.lock`
- **Package Name**: `sadie-antibody`

## Core Dependencies

### Bioinformatics
| Package | Version | Purpose |
|---------|---------|---------|
| `biopython` | >=1.80 | Sequence I/O, SeqRecord handling, FASTA parsing |
| `pyhmmer` | ^0.11.1 | HMMER3 for antibody numbering/renumbering |

### Data Processing
| Package | Version | Purpose |
|---------|---------|---------|
| `pandas` | >=1.5 | AIRR table manipulation, dataframes |
| `numpy` | * | Numerical operations |
| `pyarrow` | * | Feather file format support |
| `scipy` | ^1.11.0 | Scientific computing |
| `scikit-learn` | ^1.5.0 | ML utilities |

### Validation & Models
| Package | Version | Purpose |
|---------|---------|---------|
| `pydantic` | >=2.0.0,<3.0.0 | Data validation, Reference models |
| `PyYAML` | ^6.0 | Reference YAML parsing |

### CLI & Output
| Package | Version | Purpose |
|---------|---------|---------|
| `click` | >=8.0,<8.2 | CLI commands (`sadie airr`, `sadie reference build`) |
| `rich` | ^14.1.0 | Terminal output formatting |

### Network & Files
| Package | Version | Purpose |
|---------|---------|---------|
| `requests` | ^2.32.0 | G3 API calls (legacy), HTTP requests |
| `yarl` | ^1.9.0 | URL handling |
| `filetype` | ^1.2.0 | File type detection |
| `openpyxl` | ^3.1.0 | Excel file support |

### Utilities
| Package | Version | Purpose |
|---------|---------|---------|
| `python-Levenshtein` | ^0.27.0 | String similarity |
| `semantic-version` | ^2.10.0 | Version handling |
| `ipython` | ^8.18.0 | Interactive development |

## Development Dependencies

### Testing
| Package | Version | Purpose |
|---------|---------|---------|
| `pytest` | >=8.0.0 | Test framework |
| `pytest-cov` | ^6.2.1 | Coverage reporting |
| `coverage` | ^7.0 | Code coverage |
| `airr` | ^1.5.0 | AIRR schema validation in tests |

### Type Checking
| Package | Version | Purpose |
|---------|---------|---------|
| `mypy` | ^1.8.0 | Static type checking |
| `pyright` | ^1.1.350 | Type checking |
| `types-PyYAML` | ^6.0.12 | Type stubs |
| `pandas-stubs` | ^2.1.0 | Type stubs |
| `types-requests` | ^2.31.0 | Type stubs |

### Linting & Formatting
| Package | Version | Purpose |
|---------|---------|---------|
| `black` | ^24.0.0 | Code formatting (line-length: 120) |
| `flake8` | ^7.0.0 | Linting |
| `pre-commit` | ^3.6.0 | Git hooks |

### Documentation
| Package | Version | Purpose |
|---------|---------|---------|
| `mkdocs` | ^1.5.0 | Documentation site |
| `mkdocs-material` | ^9.5.0 | Material theme |
| `mkdocs-git-revision-date-plugin` | ^0.3.2 | Git dates in docs |

## Database Technologies

### IgBLAST
- **Binary**: Platform-specific `igblastn` (macOS/Linux)
- **Location**: `src/sadie/airr/bin/{platform}/igblastn`
- **Database Structure**:
  ```
  germlines/
  ├── Ig/
  │   ├── blastdb/{species}/{species}_V, _D, _J
  │   └── internal_data/{species}/{species}.ndm.imgt
  ├── aux_db/imgt/{species}_gl.aux
  └── .references_dataframe.csv.gz
  ```

### BLAST Databases
- **Format**: BLAST database files (.nhr, .nin, .nsq, etc.)
- **Generation**: `makeblastdb` from NCBI BLAST+
- **Used By**: IgBLAST for germline alignment

### Germline Databases (v1.2 Feature)
- **Location**: `src/sadie/germlines/`
- **Structure**:
  ```
  germlines/
  ├── sources/{provider}/{species}/*.fasta
  ├── normalized/{species}/
  └── igblast/
      ├── Ig/internal_data/{species}/
      └── aux_db/{species}_gl.aux
  ```

## Reference Module Dependencies (v1.2)

### Core Classes
- `Reference` → Pydantic validation, G3/Germlines adapter
- `References` → Multi-reference management
- `GermlineManager` → Multi-provider priority lookup
- `GermlineToG3Adapter` → Format transformation

### Validation Sources
```python
VALID_SOURCES = ["imgt", "ogrdb", "vdjbase", "custom"]
```

### Environment Variables
| Variable | Default | Purpose |
|----------|---------|---------|
| `SADIE_USE_GERMLINES_MODULE` | `true` | Use germlines vs G3 API |
| `IGDATA` | auto-detected | IgBLAST data directory |
| `TMPDIR` | system temp | Temporary file storage |

## Build & CI Configuration

### Build System
```toml
[build-system]
requires = ["poetry-core>=1.0.0", "poetry-dynamic-versioning>=1.0.0,<2.0.0"]
build-backend = "poetry_dynamic_versioning.backend"
```

### Code Style
- **Line Length**: 120 characters
- **Target Python**: 3.9, 3.10, 3.11, 3.12
- **Import Style**: isort with black profile

### Coverage Configuration
- **Source**: `sadie` package
- **Excluded**: tests, pyhmmer, scipy, numpy, pandas, sklearn, Bio
- **Reports**: HTML and XML
