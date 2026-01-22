# IMGT Germline Data

This directory contains IMGT germline sequences downloaded from the V-QUEST reference directory.

**Data Source**: https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/

## Directory Structure

```
imgt/
├── human/
│   ├── IGHV.fasta          # Ungapped (dots removed)
│   ├── IGHV_gapped.fasta   # IMGT-gapped (original)
│   ├── IGHD.fasta          # D segments (no gaps)
│   ├── IGHD_gapped.fasta
│   ├── IGHJ.fasta          # J segments (no gaps)
│   ├── IGHJ_gapped.fasta
│   ├── IGKV.fasta
│   ├── IGKV_gapped.fasta
│   ├── IGKJ.fasta
│   ├── IGKJ_gapped.fasta
│   ├── IGLV.fasta
│   ├── IGLV_gapped.fasta
│   ├── IGLJ.fasta
│   └── IGLJ_gapped.fasta
├── mouse/
│   └── ... (same structure)
└── ... (other species)
```

## How to Obtain IMGT Data

### Method 1: Automated Download (Recommended)

```bash
# Download human and mouse data
python -m sadie.germlines.scripts.download_imgt --species human mouse

# List all available species
python -m sadie.germlines.scripts.download_imgt --list-species

# Include T-cell receptor sequences
python -m sadie.germlines.scripts.download_imgt --species human --include-tr

# Force re-download
python -m sadie.germlines.scripts.download_imgt --species human --force
```

### Method 2: Manual Download from IMGT

1. Visit https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/
2. Navigate to species folder (e.g., `Homo_sapiens/IG/`)
3. Download FASTA files (IGHV.fasta, IGHD.fasta, etc.)
4. Place in appropriate species directory
5. Create ungapped versions by removing dots from sequences

## FASTA Format

### Header Format

```
>accession|gene_name|species|functionality|region|positions|length|codon_start|...
```

**Example**:
```
>M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |
caggttcagctggtgcagtctggagct...gaggtgaagaagcctggggcctcagtgaag...
```

**Key Header Fields**:
- Position 1: Accession number (GenBank/EMBL)
- Position 2: Gene name (e.g., `IGHV1-18*01`)
- Position 3: Species (e.g., `Homo sapiens`)
- Position 4: Functionality (`F` = functional, `ORF` = open reading frame, `P` = pseudogene)
- Position 5: Region type (`V-REGION`, `D-REGION`, `J-REGION`)

### Sequence Characteristics

- **V-regions**: IMGT-gapped with dots (`.`) per IMGT unique numbering (positions 1-104)
- **D/J-regions**: NOT gapped (IMGT numbering only applies to V regions)
- **Case**: Lowercase nucleotides
- **Gap character**: Period (`.`)

## Data Volume

| Species | Total | IGHV | IGHD | IGHJ | IGKV | IGKJ | IGLV | IGLJ |
|---------|-------|------|------|------|------|------|------|------|
| Human | 794 | 460 | 47 | 15 | 132 | 10 | 119 | 11 |
| Mouse | 953 | 678 | 61 | 9 | 168 | 10 | 19 | 8 |

## Important Notes

- **Gapped vs Ungapped**: `*_gapped.fasta` contains original IMGT format with dots; `*.fasta` has dots removed
- **Functionality**: Filter by functionality code if you only want functional genes (`F`)
- **Multiple alleles**: Each gene may have multiple alleles (e.g., `*01`, `*02`, `*03`)
- **Partial sequences**: Some sequences are partial (noted in header as "partial in 3'" or similar)

## Validation

After downloading, validate with:

```python
from sadie.germlines import get_manager

manager = get_manager()
genes = manager.get_genes("human", "V", "H")
print(f"Loaded {len(genes)} IGHV genes from IMGT")
```

## References

- IMGT Homepage: https://www.imgt.org/
- V-QUEST Reference Directory: https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/
- IMGT/GENE-DB: https://www.imgt.org/genedb/
- IMGT Scientific Chart: https://www.imgt.org/IMGTScientificChart/
