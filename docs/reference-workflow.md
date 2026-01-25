# Custom Reference Databases

Learn how to create custom germline reference databases for AIRR annotation using YAML configuration files.

## Overview

SADIE supports creating custom reference databases that combine germline sequences from multiple sources:

| Source | Description | Website | Best For |
|--------|-------------|---------|----------|
| **imgt** | International ImMunoGeneTics | [imgt.org](https://www.imgt.org) | Comprehensive baseline |
| **ogrdb** | Open Germline Receptor Database | [ogrdb.airr-community.org](https://ogrdb.airr-community.org/) | Novel/validated alleles |
| **vdjbase** | VDJbase curated germlines | [vdjbase.org](https://www.vdjbase.org/) | Population-specific alleles |
| **custom** | Your own sequences | - | Novel discoveries |

!!! tip "New in v1.2"
    SADIE now supports **OGRDB** and **VDJbase** as first-class data sources alongside IMGT. You can mix sources within a single reference database.

This enables:

- Focusing on specific genes relevant to your research
- Combining curated alleles from different databases
- Including custom/novel germline sequences
- Creating species-specific or project-specific databases

## Workflow

```
Write YAML → Build Database → Use in Airr → Annotate Sequences
```

### Step 1: Write YAML Configuration

Create a YAML file defining which genes to include from which sources.

**Structure:**
```yaml
reference_name:        # Name for this database (used in code)
  source:              # Provider: imgt, ogrdb, vdjbase, custom
    species:           # Species name
      - GENE*ALLELE    # List of genes with allele numbers
```

**Example** (`my-reference.yml`):
```yaml
my-project:
  imgt:
    human:
      - IGHV1-2*02
      - IGHV3-23*01
      - IGHD3-10*01
      - IGHJ4*02
```

**Tip:** Always include the allele number (e.g., `*01`, `*02`). Find gene names at:

- IMGT genes: [IMGT Gene Database](https://www.imgt.org/genedb/)
- OGRDB genes: [OGRDB](https://ogrdb.airr-community.org/)

### Step 2: Build the Database

#### Option A: CLI (Recommended)

```bash
# New build command (v1.2+)
sadie reference build my-reference.yml --output ./my-databases

# Or legacy command
sadie reference make --reference my-reference.yml --outpath ./my-databases
```

#### Option B: Python API

```python
from sadie.reference import References

# Load and build from YAML
references = References.from_yaml("my-reference.yml")

# Build the database
output_path = references.make_airr_database("./my-databases")
print(f"Database built at: {output_path}")
```

**Output structure:**
```
my-databases/
├── Ig/
│   ├── blastdb/
│   │   └── my-project/
│   │       ├── my-project_V.*
│   │       ├── my-project_D.*
│   │       └── my-project_J.*
│   └── internal_data/
│       └── my-project/
│           └── my-project.ndm.imgt
├── aux_db/
│   └── imgt/
│       └── my-project_gl.aux
└── .references_dataframe.csv.gz
```

### Step 3: Use in AIRR Annotation

#### Option A: Using Prebuilt Database (v1.2+)

```python
from sadie.airr import Airr

# Use prebuilt database directly
airr = Airr("my-project", database="./my-databases")

# Annotate sequences
result = airr.run_file("sequences.fasta")
print(result)
```

#### Option B: Using References Object

```python
from sadie.airr import Airr
from sadie.reference import References

# Load custom reference
references = References(default_output_path="./my-databases")

# Create Airr annotator
airr = Airr(reference_name="my-project", references=references)

# Annotate sequences
result = airr.run_file("sequences.fasta")
```

## Complete Example

### 1. Create reference.yml

```yaml
# research-project.yml
therapeutic-antibodies:
  imgt:
    human:
      # Common therapeutic V genes
      - IGHV1-69*01
      - IGHV3-23*01
      - IGHV4-34*01
      # D genes
      - IGHD3-10*01
      - IGHD3-22*01
      # J genes
      - IGHJ4*02
      - IGHJ6*02
      # Kappa light chain
      - IGKV1-39*01
      - IGKV3-20*01
      - IGKJ1*01
```

### 2. Build Database

```bash
sadie reference build research-project.yml --output ./therapeutic-db
```

### 3. Run Annotation

```python
from sadie.airr import Airr

# Use prebuilt database
airr = Airr("therapeutic-antibodies", database="./therapeutic-db")

# Annotate your sequences
results = airr.run_file("my-sequences.fasta")

# View results
print(results[['sequence_id', 'v_call', 'd_call', 'j_call', 'productive']])
```

**Example Output:**
```
   sequence_id      v_call     d_call    j_call  productive
0  seq_001     IGHV1-69*01  IGHD3-10*01  IGHJ4*02       True
1  seq_002     IGHV3-23*01  IGHD3-22*01  IGHJ6*02       True
```

## Multi-Source Configuration

Combine genes from different providers in one reference:

```yaml
comprehensive-human:
  # Baseline from IMGT
  imgt:
    human:
      - IGHV1-2*02
      - IGHV1-18*01
      - IGHD3-10*01
      - IGHJ4*02
  
  # Add curated alleles from OGRDB
  ogrdb:
    human:
      - IGHV1-69*13  # Novel allele from OGRDB
      - IGHV4-34*09  # Curated variant
```

## Using OGRDB Alleles

[OGRDB](https://ogrdb.airr-community.org/) provides curated germline alleles that have been validated through the AIRR Community review process. Use OGRDB when you need:

- Novel alleles not yet in IMGT
- Alleles with experimental validation
- Population-specific variants

```yaml
# Example: Human OGRDB-only reference
human-ogrdb-curated:
  ogrdb:
    human:
      # OGRDB-validated V genes
      - IGHV1-69*01
      - IGHV1-69*13    # Novel allele
      - IGHV3-30*18    # Validated variant
      - IGHV4-34*09
      # D genes
      - IGHD3-10*01
      - IGHD3-22*01
      # J genes  
      - IGHJ4*02
      - IGHJ6*02
```

!!! note "OGRDB Gene Names"
    OGRDB uses the same naming convention as IMGT. Browse available alleles at [ogrdb.airr-community.org](https://ogrdb.airr-community.org/).

## Using VDJbase Alleles

[VDJbase](https://www.vdjbase.org/) provides population-level germline inference data. Use VDJbase when you need:

- Population-specific allele frequencies
- Haplotype-aware references
- Alleles inferred from repertoire data

```yaml
# Example: Human VDJbase reference
human-vdjbase:
  vdjbase:
    human:
      # VDJbase-inferred alleles
      - IGHV1-2*02
      - IGHV1-69*01
      - IGHV3-23*01
      - IGHD2-15*01
      - IGHJ4*02
```

## Combining All Sources

For maximum coverage, combine IMGT baseline with curated alleles:

```yaml
# Complete multi-source reference
human-complete:
  # IMGT baseline - comprehensive coverage
  imgt:
    human:
      - IGHV1-2*02
      - IGHV1-18*01
      - IGHV3-23*01
      - IGHD3-10*01
      - IGHD3-22*01
      - IGHJ4*02
      - IGHJ6*02
      - IGKV1-39*01
      - IGKJ1*01
  
  # OGRDB - add validated novel alleles
  ogrdb:
    human:
      - IGHV1-69*13    # Novel validated allele
      - IGHV4-34*09    # Population variant
  
  # VDJbase - add population-inferred alleles
  vdjbase:
    human:
      - IGHV3-30*18    # Inferred from repertoires
```

!!! warning "Source Priority"
    When the same allele appears in multiple sources, SADIE uses explicit source selection - each gene is fetched from its specified source. There is no automatic fallback between sources.

## Multi-Species Configuration

Create databases spanning multiple species:

```yaml
cross-species:
  imgt:
    human:
      - IGHV1-69*01
      - IGHD3-10*01
      - IGHJ4*02
    mouse:
      - IGHV1-4*01
      - IGHD1-1*01
      - IGHJ1*01
```

## Validation

### Check YAML Syntax

```bash
python -c "import yaml; yaml.safe_load(open('my-reference.yml'))"
```

### Verify Built Database

```bash
# Check database files exist
ls -la ./my-databases/Ig/blastdb/*/

# View reference dataframe
python -c "import pandas as pd; print(pd.read_csv('./my-databases/.references_dataframe.csv.gz').head())"
```

## Troubleshooting

### Gene Not Found

```
G3Error: Gene IGHV1-2*99 not found in imgt database for species human
```

**Solution:** Verify the gene name and allele exist in the source database.

### Missing D-REGION Error

```
ValueError: No D-REGION found in reference object
```

**Solution:** Include D genes in your YAML. Heavy chain annotation requires V, D, and J genes.

### Species Not Found

**Solution:** Use standard species names: `human`, `mouse`, `rat`, `rabbit`, `rhesus_macaque`, etc.

## Sample Files

Download sample configurations:

- [`reference-sample.yml`](../examples/reference-sample.yml) - Multi-source example

## See Also

- [Germline Management](germlines/index.md) - Managing germline databases
- [Provider Guide](germlines/provider-guide.md) - Understanding data sources
