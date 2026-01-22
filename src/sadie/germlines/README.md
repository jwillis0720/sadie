# Sadie Germlines Module

**Self-contained germline database management for immunoglobulin repertoires.**

This module replaces G3 API dependency with local, multi-source germline database management.

## Overview

The germlines module provides:

✅ **Multiple Data Sources**: IMGT, OGRDB, VDJbase, and custom user sequences
✅ **Priority-Based Merging**: Custom overrides OGRDB overrides VDJbase overrides IMGT
✅ **Automatic Processing**: Auto-gaps ungapped sequences
✅ **IgBLAST Integration**: Builds databases and auxiliary files
✅ **Offline Operation**: Works without internet after initial setup
✅ **Self-Contained**: Can be extracted as standalone package

## Architecture

### Staged Pipeline

```
sources/          normalized/       igblast/
├── custom/  →   ├── gapped/   →   ├── database/
├── imgt/         └── ungapped/     ├── aux_db/
└── ogrdb/                          └── internal_data/

Stage 1:          Stage 2:          Stage 3:
Fetch & Parse     Merge & Gap       Build DBs
```

**Stage 1: Sources**
- Raw FASTA files from each source
- User adds custom sequences here
- IMGT and OGRDB data downloaded once

**Stage 2: Normalized**
- Merged sequences (priority-based deduplication)
- Both gapped and ungapped versions
- Change detection triggers rebuild

**Stage 3: IgBLAST**
- BLAST databases for V, D, J segments
- Auxiliary files with CDR/FWR annotations
- Ready for IgBLAST execution

### Priority System

Default priority order: **custom > ogrdb > vdjbase > imgt**

Deduplication rules:
1. Same gene name → first provider wins
2. Same exact sequence → first provider wins
3. Novel gene → include from any source

```python
# Default priority
from sadie.germlines import get_germline_genes
genes = get_germline_genes("human", "V", "H")
# Uses: custom, then IMGT, then OGRDB

# Custom priority
from sadie.germlines import GermlineManager
manager = GermlineManager(providers=["ogrdb", "imgt"])
genes = manager.get_genes("human", "V", "H")
# Uses only OGRDB and IMGT (no custom)
```

## Directory Structure

```
germlines/
├── __init__.py              # Public API
├── README.md                # This file
├── models.py                # Data models
├── manager.py               # Priority-based manager
├── pipeline.py              # Staged processing
│
├── providers/               # Data source abstractions
│   ├── base.py             # Provider interface
│   ├── custom.py           # Custom sequences
│   ├── imgt.py             # IMGT provider
│   ├── ogrdb.py            # OGRDB provider
│   └── vdjbase.py          # VDJbase provider
│
├── builders/                # Database builders
│   ├── blast.py            # BLAST database builder
│   └── aux.py              # Auxiliary file builder
│
├── sources/                 # Raw data (Stage 1)
│   ├── custom/             # USER EDITABLE
│   │   ├── README.md
│   │   ├── _template/
│   │   └── human/
│   │       ├── IGHV.fasta  ← Add your sequences here
│   │       └── ...
│   ├── imgt/               # IMGT data
│   │   ├── IMGT_DATA.md   # How to populate
│   │   └── human/
│   │       ├── IGHV.fasta
│   │       └── ...
│   └── ogrdb/              # OGRDB data
│       ├── OGRDB_DATA.md  # How to populate
│       └── human/
│           └── ...
│
├── normalized/              # Processed data (Stage 2)
│   └── human/
│       ├── gapped/         # IMGT-gapped sequences
│       │   ├── IGHV.fasta
│       │   └── ...
│       └── ungapped/       # Ungapped sequences
│           ├── IGHV.fasta
│           └── ...
│
├── igblast/                 # IgBLAST format (Stage 3)
│   ├── database/
│   │   └── human/
│   │       ├── human_V.fasta
│   │       ├── human_V.n{hr,in,sq}
│   │       ├── human_D.*
│   │       └── human_J.*
│   ├── aux_db/
│   │   └── human_gl.aux
│   └── internal_data/      # IgBLAST internal data
│
├── scripts/                 # Maintenance scripts
│   ├── download_imgt.py
│   ├── download_ogrdb.py
│   └── validate.py
│
└── tests/                   # Unit tests
    ├── test_manager.py
    ├── test_pipeline.py
    └── test_providers.py
```

## Quick Start

### 1. Add Custom Sequences

```bash
# Create file
$ cat > sources/custom/human/IGHV.fasta << 'EOF'
>IGHV1-NOVEL*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCT...
EOF
```

### 2. Add IMGT Data

See `sources/imgt/IMGT_DATA.md` for instructions on downloading IMGT data.

### 3. Add OGRDB Data (Optional)

See `sources/ogrdb/OGRDB_DATA.md` for instructions.

### 4. Run Sadie

```python
from sadie.germlines import get_germline_genes

# Get genes (automatically processes if changed)
genes = get_germline_genes("human", "V", "H")

print(f"Found {len(genes)} genes:")
for gene in genes[:5]:
    print(f"  {gene.name} ({gene.source})")
```

### 5. Use with IgBLAST

The module automatically builds IgBLAST databases. Update your Sadie config:

```python
# See INTEGRATION_GUIDE.md for details
from pathlib import Path

germline_base = Path("src/sadie/germlines/igblast")
v_db = germline_base / "database" / "human" / "human_V"
aux_file = germline_base / "aux_db" / "human_gl.aux"
```

## API Documentation

### Simple API

```python
from sadie.germlines import get_germline_genes, get_gene_by_name

# Get all genes for segment/chain
genes = get_germline_genes("human", "V", "H")

# Get specific gene
gene = get_gene_by_name("IGHV1-69*01", "human")

# Custom providers
genes = get_germline_genes(
    "human", "V", "H",
    providers=["custom", "imgt"]  # Skip OGRDB
)
```

### Advanced API

```python
from sadie.germlines import GermlineManager, GermlinePipeline
from pathlib import Path

# Custom priority order
manager = GermlineManager(providers=["ogrdb", "imgt", "custom"])
genes = manager.get_genes("human", "V", "H")

# Manual pipeline control
pipeline = GermlinePipeline(Path("src/sadie/germlines"))
pipeline.force_rebuild("human")  # Force rebuild

# Check available species
species = manager.get_available_species()
print(f"Available species: {species}")
```

## Change Detection

The pipeline automatically detects when files change:

```bash
# Add new custom sequence
$ echo ">IGHV1-NEW*01" >> sources/custom/human/IGHV.fasta
$ echo "CAGGTGCAG..." >> sources/custom/human/IGHV.fasta

# Next run automatically rebuilds
$ sadie airr -i sequences.fasta -o results.tsv
# INFO: Sources changed, rebuilding normalized...
# INFO: Normalized changed, rebuilding IgBLAST...
```

No manual rebuild needed!

## Data Model

### GermlineGene

Unified gene model across all providers:

```python
@dataclass
class GermlineGene:
    name: str                    # IGHV1-69*01
    species: str                 # human
    segment: str                 # V, D, or J
    chain: str                   # H, K, or L
    sequence: str                # Ungapped
    sequence_gapped: str         # IMGT-gapped (optional)
    is_functional: bool          # Functional gene?
    functionality: str           # F, ORF, P
    source: str                  # custom, imgt, ogrdb
    regions: dict                # CDR/FWR regions (optional)
```

## Integration with Sadie

See `INTEGRATION_GUIDE.md` for complete integration instructions.

Summary:

1. Update `src/sadie/airr/igblast/germline.py` to use new paths
2. Update `src/sadie/reference/reference.py` to use germlines module
3. Update `src/sadie/renumbering/aligners/hmmer.py` to use new HMM builder

## Standalone Use

This module can be extracted as standalone package:

```bash
# Future: extract as separate repo
git subtree split --prefix=src/sadie/germlines -b germlines-standalone
```

Use in other projects:

```python
from germlines import get_germline_genes

genes = get_germline_genes("mouse", "V", "H")
```

## Testing

```bash
# Run tests
pytest tests/

# Validate installation
python -c "from sadie.germlines import get_germline_genes; print('OK')"

# Check data
python scripts/validate.py
```

## Performance

- **First run**: Downloads/processes data (~1-2 minutes for human)
- **Subsequent runs**: Uses cached data (instant)
- **Change detection**: Only rebuilds what changed
- **Offline**: Works completely offline after initial setup

## Troubleshooting

### No genes found

**Check**:
1. Are source files present? `ls sources/*/human/`
2. Are they named correctly? `IG{H|K|L}{V|D|J}.fasta`
3. Try force rebuild: `pipeline.force_rebuild("human")`

### Custom sequences not loading

**Check**:
1. File in correct directory? `sources/custom/human/`
2. Valid FASTA format? Headers start with `>`
3. Valid nucleotides? Only ACGT
4. Check logs: `import logging; logging.basicConfig(level=logging.DEBUG)`

### IgBLAST errors

**Check**:
1. BLAST databases built? `ls igblast/database/human/*.nhr`
2. Auxiliary file exists? `ls igblast/aux_db/human_gl.aux`
3. Run pipeline: `pipeline.update("human")`

## Contributing

This module follows:
- **Zen of Python**: Simple, explicit, readable
- **PEP 8**: Code style
- **Type hints**: All public functions
- **Docstrings**: Google style

## References

- **IMGT**: https://www.imgt.org/
- **OGRDB**: https://ogrdb.airr-community.org/
- **IgBLAST**: https://www.ncbi.nlm.nih.gov/igblast/
- **AIRR Standards**: https://docs.airr-community.org/

## Support

- **Issues**: https://github.com/jwillis0720/sadie/issues
- **Documentation**: https://sadie.jordanrwillis.com/
- **Email**: jwillis0720@gmail.com

## License

Same as Sadie (MIT License)
