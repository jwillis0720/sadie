# SADIE Directory and File Structure

## Project Root

```
/Users/tmsincomb/sadie/
├── src/sadie/           # Main source code
├── tests/               # Test suite
├── docs/                # Documentation (MkDocs)
├── examples/            # Usage examples
├── notebooks/           # Jupyter notebooks
├── scripts/             # Utility scripts
├── specs/               # Specification documents
├── reference.yml        # Default reference configuration
├── pyproject.toml       # Poetry project config
├── pytest.ini           # pytest configuration
├── mkdocs.yml           # MkDocs configuration
└── .planning/           # Planning documentation (generated)
```

## Source Code Structure

### Main Package (`src/sadie/`)

```
src/sadie/
├── __init__.py          # Package exports: airr, numbering, receptor, reference, renumbering
├── app.py               # Streamlit web application
│
├── airr/                # AIRR-standard annotation (PRIMARY ENTRY POINT)
│   ├── __init__.py
│   ├── airr.py          # Airr class - main annotation API
│   ├── methods.py       # Annotation helper methods
│   ├── exceptions.py    # Custom exceptions (BadDataSet, BadIgBLASTExe, etc.)
│   ├── airrtable/       # AirrTable and LinkedAirrTable classes
│   ├── igblast/         # IgBLAST wrapper
│   │   ├── __init__.py
│   │   ├── igblast.py   # IgBLASTN class - blast execution
│   │   └── germline.py  # GermlineData - database path management
│   ├── models/          # Data models
│   ├── bin/             # Platform-specific igblastn binaries
│   └── data/            # Default germline data (legacy)
│
├── reference/           # YAML → IgBLAST database builder
│   ├── __init__.py      # Exports: Reference, References, YamlRef
│   ├── reference.py     # Reference and References classes
│   ├── models.py        # Pydantic models (GeneEntry, GeneEntries)
│   ├── yaml.py          # YamlRef - YAML configuration parser
│   ├── genbank.py       # GenBank utilities
│   ├── settings.py      # Reference settings
│   ├── util.py          # Utilities (make_blast_db, write_fasta)
│   ├── bin/             # makeblastdb binaries
│   └── data/            # Default reference.yml
│
├── germlines/           # Multi-source germline database (v1.2)
│   ├── __init__.py      # Public API: get_germline_genes, GermlineManager, etc.
│   ├── manager.py       # GermlineManager - priority-based lookup
│   ├── models.py        # GermlineGene, ProviderMetadata
│   ├── pipeline.py      # GermlinePipeline - normalize → build workflow
│   ├── g3_adapter.py    # GermlineToG3Adapter - format conversion
│   ├── cli.py           # CLI commands
│   ├── renumbering_integration.py  # LocalHMMBuilder for renumbering
│   ├── providers/       # Database providers
│   │   ├── __init__.py
│   │   ├── base.py      # GermlineProvider abstract class
│   │   ├── imgt.py      # IMGTProvider
│   │   ├── ogrdb.py     # OGRDBProvider
│   │   ├── vdjbase.py   # VDJbaseProvider
│   │   └── custom.py    # CustomProvider
│   ├── builders/        # IgBLAST database builders
│   ├── sources/         # Raw data by provider
│   │   ├── imgt/
│   │   ├── ogrdb/
│   │   ├── vdjbase/
│   │   └── custom/
│   ├── normalized/      # Merged/processed sequences
│   ├── igblast/         # IgBLAST-ready databases
│   │   ├── Ig/
│   │   │   └── internal_data/{species}/
│   │   └── aux_db/
│   ├── hmms/            # Generated HMM files
│   ├── stockholms/      # Stockholm alignment files
│   ├── scripts/         # Data processing scripts
│   └── utils/           # Germlines utilities
│
├── renumbering/         # Antibody numbering module
│   ├── __init__.py
│   ├── renumbering.py   # Main renumbering logic
│   ├── result.py        # Renumbering results
│   ├── constants.py     # Numbering scheme constants
│   ├── exception.py     # Custom exceptions
│   ├── numbering_translator.py  # Scheme translation
│   ├── aligners/        # Alignment implementations
│   │   └── hmmer.py     # HMMER-based alignment
│   ├── clients/         # External service clients
│   └── data/            # Reference data
│
├── numbering/           # Legacy numbering (being deprecated)
│
├── receptor/            # Receptor utilities
│
├── cluster/             # Sequence clustering
│
├── typing/              # Type definitions
│   └── __init__.py      # Species, Chain, Source types
│
└── utility/             # Shared utilities
```

## Tests Structure

```
tests/
├── __init__.py
├── conftest.py          # pytest fixtures
├── data/                # Test data files
├── integration/         # Integration tests
└── unit/                # Unit tests
    ├── airr/            # AIRR module tests
    │   ├── test_airr.py
    │   └── test_methods.py
    ├── germlines/       # Germlines module tests
    │   ├── test_reference_integration.py
    │   └── ... (20+ test files)
    ├── reference/       # Reference module tests
    ├── renumbering/     # Renumbering tests
    ├── aligners/        # Aligner tests
    ├── cluster/         # Cluster tests
    ├── receptor/        # Receptor tests
    ├── typing/          # Typing tests
    └── utility/         # Utility tests
```

## Key Entry Points

### For Annotation
```python
# Primary: Airr class
from sadie.airr import Airr
airr = Airr("human")
result = airr.run_fasta("sequences.fasta")

# With prebuilt database
airr = Airr("custom", database="/path/to/database")
```

### For Reference Building
```python
# Build custom reference database
from sadie.reference import References
refs = References.from_yaml("reference.yml", use_germlines=True)
refs.make_airr_database("/output/path")
```

### For Germlines Access
```python
# Direct gene access
from sadie.germlines import get_germline_genes, GermlineManager

# Simple API
genes = get_germline_genes("human", "V", "H")

# Advanced with custom providers
manager = GermlineManager(providers=["custom", "imgt"])
genes = manager.get_genes("human", "V", "H")
```

## Configuration Files

### reference.yml Format
```yaml
# Reference name (becomes database name)
human:
  # Source database
  imgt:
    # Species
    human:
      - IGHV1-69*01
      - IGHD3-3*01
      - IGHJ6*01
      # ... V, D, J genes required

# Chimeric example (multiple species)
chimeric_ref:
  imgt:
    human:
      - IGHV1-69*01
    mouse:
      - IGHV1-18*01
```

### IgBLAST Database Structure (Generated)
```
{output_path}/
├── Ig/
│   ├── blastdb/
│   │   └── {name}/
│   │       ├── {name}_V.nhr
│   │       ├── {name}_V.nin
│   │       ├── {name}_V.nsq
│   │       ├── {name}_D.*
│   │       └── {name}_J.*
│   └── internal_data/
│       └── {name}/
│           ├── {name}.ndm.imgt    # FWR/CDR boundaries
│           └── {name}_V.*         # V gene database
├── aux_db/
│   └── imgt/
│       └── {name}_gl.aux          # J gene reading frames
└── .references_dataframe.csv.gz   # Full reference data
```

## Where to Put New Code

### New Provider
```
src/sadie/germlines/providers/{provider_name}.py
- Extend GermlineProvider base class
- Implement: fetch_genes, fetch_gene_by_name, get_metadata, is_available, download
- Register in manager._create_provider()
```

### New Data Source
```
src/sadie/germlines/sources/{source_name}/
- Raw data files organized by species
- Must follow FASTA naming: IG{chain}{segment}.fasta
```

### New AIRR Feature
```
src/sadie/airr/airr.py           # Main Airr class
src/sadie/airr/airrtable/        # AirrTable extensions
src/sadie/airr/methods.py        # Helper methods
```

### New Reference Configuration
```
src/sadie/reference/data/        # Default configs
reference.yml (project root)     # Project-specific config
```

### New Tests
```
tests/unit/{module}/test_{feature}.py    # Unit tests
tests/integration/                        # Integration tests
tests/data/                              # Test fixtures
```

## Critical Files for Reference Module (v1.2)

| File | Purpose |
|------|---------|
| `reference/reference.py` | Reference and References classes - YAML → database |
| `reference/models.py` | Pydantic validation (GeneEntry, GeneEntries) |
| `germlines/manager.py` | GermlineManager - multi-source priority lookup |
| `germlines/g3_adapter.py` | GermlineGene → G3 format adapter |
| `germlines/providers/base.py` | GermlineProvider abstract interface |
| `airr/igblast/germline.py` | GermlineData - database path resolution |
| `airr/airr.py` | Airr class - uses database or germlines |
| `reference.yml` | Default reference configuration |
