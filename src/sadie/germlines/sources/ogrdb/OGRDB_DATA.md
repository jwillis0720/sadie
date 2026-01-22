# OGRDB Germline Data

This directory contains OGRDB germline sequences from the Open Germline Receptor Database.

**Data Source**: https://ogrdb.airr-community.org/
**Zenodo Archive**: https://zenodo.org/records/18145568

## Directory Structure

```
ogrdb/
├── OGRDB_DATA.md          # This file
├── human/
│   ├── IGHV.fasta          # Ungapped V genes (heavy chain)
│   ├── IGHV_gapped.fasta   # IMGT-gapped V genes (heavy chain)
│   ├── IGHD.fasta          # D genes (heavy chain)
│   ├── IGHJ.fasta          # J genes (heavy chain)
│   ├── IGHJ_gapped.fasta   # IMGT-gapped J genes
│   ├── IGKV.fasta          # Kappa V genes
│   ├── IGKJ.fasta          # Kappa J genes
│   ├── IGLV.fasta          # Lambda V genes
│   └── IGLJ.fasta          # Lambda J genes
├── mouse/
│   └── ...
└── ... (other species)
```

## How to Obtain OGRDB Data

### Method 1: Automated Download (Recommended)

Download OGRDB data from the Zenodo archive:

```bash
# Download human data
python -m sadie.germlines.scripts.download_ogrdb --species human

# Download multiple species
python -m sadie.germlines.scripts.download_ogrdb --species human mouse

# Force re-download (bypass cache)
python -m sadie.germlines.scripts.download_ogrdb --species human --force
```

The download script:
1. Downloads the OGRDB archive from Zenodo (~200MB)
2. Extracts the SQL dump
3. Parses gapped (`coding_seq_imgt`) and ungapped (`sequence`) sequences
4. Writes organized FASTA files per species/segment/chain

### Method 2: Using Python API

```python
from sadie.germlines.providers.ogrdb import OGRDBProvider

provider = OGRDBProvider()
provider.download(["human", "mouse"])
```

### Method 3: Manual Download from OGRDB Website

1. Visit https://ogrdb.airr-community.org/
2. Browse germline sets
3. Select species and locus
4. Download sequences in FASTA format
5. Save with naming convention: `IG{chain}{segment}.fasta`

## FASTA Format

### Ungapped FASTA (e.g., IGHV.fasta)

```
>IGHV1-69*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCT...
>IGHV3-23*01
GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCT...
```

### Gapped FASTA (e.g., IGHV_gapped.fasta)

```
>IGHV1-69*01
CAGGTGCAGCTGGTGCAG......TCTGGGGCT...GAGGTGAAGAAGCCT...
>IGHV3-23*01
GAGGTGCAGCTGTTGGAG......TCTGGGGGAGGC...TTGGTACAGCCT...
```

**Key Features:**
- Headers contain gene name (allele designation)
- Gapped sequences use `.` (period) as gap character per IMGT convention
- Novel alleles discovered through repertoire sequencing

## Data from SQL Dump

The Zenodo archive contains an SQL dump with the `gene_description` table:

| Column | Description |
|--------|-------------|
| `sequence_name` | Gene allele name (e.g., IGHV1-69*01) |
| `species` | Species name (e.g., Homo sapiens) |
| `sequence` | Ungapped nucleotide sequence |
| `coding_seq_imgt` | IMGT-gapped nucleotide sequence |

## Priority with IMGT

Default priority order: `custom > ogrdb > vdjbase > imgt`

- If a gene exists in both OGRDB and IMGT → OGRDB version is used (higher priority)
- Novel genes only in OGRDB → included
- Customize priority:

```python
from sadie.germlines import GermlineManager

# Use IMGT first, then OGRDB for novel alleles
manager = GermlineManager(providers=["custom", "imgt", "ogrdb"])
genes = manager.get_genes("human", "V", "H")
```

## Validation

After downloading:

```python
from sadie.germlines import get_manager

manager = get_manager()
genes = manager.get_genes("human", "V", "H")

# Check how many came from OGRDB
ogrdb_genes = [g for g in genes if g.source == "ogrdb"]
print(f"Loaded {len(ogrdb_genes)} OGRDB genes")

# Check for gapped sequences
gapped = [g for g in ogrdb_genes if g.sequence_gapped]
print(f"  {len(gapped)} with IMGT gapping")
```

## Important Notes

- **Community-curated**: OGRDB sequences are reviewed by the AIRR Community
- **Novel alleles**: Contains alleles discovered via repertoire sequencing not in IMGT
- **Evidence-based**: Each sequence has supporting evidence from repertoire data
- **Gapped sequences**: Both gapped and ungapped versions extracted from SQL dump
- **D segments**: No gapped versions (D segments are not gapped in IMGT)

## References

- OGRDB Homepage: https://ogrdb.airr-community.org/
- OGRDB API Docs: https://ogrdb.airr-community.org/api/docs
- Zenodo Archive: https://zenodo.org/records/18145568
- AIRR Community: https://www.antibodysociety.org/the-airr-community/
- Publication: https://academic.oup.com/nar/article/48/D1/D964/5576123
