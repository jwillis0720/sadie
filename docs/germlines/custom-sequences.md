# Custom Sequences: Adding Your Own Germlines

## Overview

Custom sequences allow you to add your own novel germline alleles to SADIE's annotation databases. This is essential when you've discovered new alleles not yet in public databases (IMGT, OGRDB) or when working with proprietary germline data.

**Key Features:**

- **Highest Priority**: Custom sequences override all provider sequences
- **Simple Format**: Standard FASTA files
- **Species-Specific**: Organized by species and locus
- **Automatic Integration**: Merged during population/rebuild

---

## When to Add Custom Sequences

### Use Cases

**Novel Allele Discovery:**

You've identified a new germline allele through repertoire sequencing analysis that isn't in IMGT or OGRDB.

```
# Example: New IGHV3-30 variant
>IGHV3-30*20_custom
CAGGTGCAGCTGGTGCAGTCTGGGGCT...
```

**Model Organisms:**

Working with a species or strain not covered by public databases.

```
# Example: Custom mouse strain
>IGHV1-1*01_custom_strain
CAGGTCCAGCTGGTGCAGTCTGGAGCT...
```

**Validation Studies:**

Testing specific sequence variants or mutations for research purposes.

**Proprietary Data:**

Commercial or unpublished germline sequences that you have rights to use.

---

## File Format Requirements

### FASTA Format

Custom sequences must be in standard FASTA format:

```
>SEQUENCE_NAME
ATCGATCGATCGATCG...
```

**Requirements:**

1. **Header line** starts with `>`
2. **Sequence name** follows IMGT nomenclature (recommended)
3. **Sequence** on subsequent lines
4. **Nucleotide sequences** (A, T, C, G, N allowed)
5. **No amino acid sequences** (protein sequences not supported)

### Naming Conventions

**Recommended Format:**

```
>LOCUS_GENE*ALLELE[_suffix]
```

**Examples:**

```
>IGHV3-30*20
>IGHV1-2*05_novel
>IGHV4-59*01_custom
>IGHD3-10*02
>IGHJ4*03
```

**Guidelines:**

- Follow IMGT nomenclature when possible
- Use `*` to separate gene and allele numbers
- Optional `_suffix` for tracking (e.g., `_custom`, `_novel`, `_lab`)
- Avoid special characters except `*`, `-`, `_`

### Example Custom FASTA File

```fasta
>IGHV3-30*20_novel
CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTT
ACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACA
AACTATGCACAGAAGTTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGA
TCTGACGACACGGCCGTGTATTACTGTGCGAGAGA

>IGHV1-2*05_lab
CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTTTCCTGCAAGGCATCTGGATACACCTTC
ACCGGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGACGGATCAACCCTAACAGTGGTGGCACA
AACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGA
TCTGACGACACGGCCGTGTATTACTGTGCGAGAGA

>IGHD3-10*02_variant
GTATTACTATGATTACGACTGGAAC
```

---

## Directory Structure

### Location

Custom sequences are stored in:

```
~/.sadie/germlines/sources/custom/<species>/<locus>.fasta
```

**Example Structure:**

```
~/.sadie/germlines/sources/
├── imgt/
│   └── human/
│       ├── IGHV.fasta
│       ├── IGHD.fasta
│       └── IGHJ.fasta
├── ogrdb/
│   └── human/
│       └── IGHV.fasta
└── custom/
    ├── human/
    │   ├── IGHV.fasta    ← Your custom heavy chain V genes
    │   ├── IGHD.fasta    ← Your custom heavy chain D genes
    │   └── IGHJ.fasta    ← Your custom heavy chain J genes
    └── mouse/
        └── IGHV.fasta    ← Custom mouse sequences
```

### Locus Files

**Heavy Chain:**
- `IGHV.fasta` - Heavy chain V genes
- `IGHD.fasta` - Heavy chain D genes
- `IGHJ.fasta` - Heavy chain J genes

**Light Chain:**
- `IGKV.fasta` - Kappa V genes
- `IGKJ.fasta` - Kappa J genes
- `IGLV.fasta` - Lambda V genes
- `IGLJ.fasta` - Lambda J genes

---

## Step-by-Step Instructions

### Step 1: Prepare Your FASTA File

Create a FASTA file with your custom sequences:

```bash
# Create file on your computer
nano my_custom_germlines.fasta
```

Add your sequences:

```fasta
>IGHV3-30*20_novel
CAGGTGCAGCTGGTGCAGTCTGGGGCT...

>IGHV1-2*05_lab
CAGGTTCAGCTGGTGCAGTCTGGAGCT...
```

Save the file.

### Step 2: Create Custom Directory

```bash
# Create directory for your species
mkdir -p ~/.sadie/germlines/sources/custom/human
```

**For other species:**

```bash
# Mouse
mkdir -p ~/.sadie/germlines/sources/custom/mouse

# Rabbit
mkdir -p ~/.sadie/germlines/sources/custom/rabbit
```

### Step 3: Copy FASTA Files

Copy your FASTA file to the appropriate locus file:

```bash
# For heavy chain V genes
cp my_custom_germlines.fasta ~/.sadie/germlines/sources/custom/human/IGHV.fasta

# For heavy chain D genes
cp my_custom_d_genes.fasta ~/.sadie/germlines/sources/custom/human/IGHD.fasta

# For heavy chain J genes
cp my_custom_j_genes.fasta ~/.sadie/germlines/sources/custom/human/IGHJ.fasta
```

**Or combine multiple files:**

```bash
# Append to existing custom file
cat my_additional_sequences.fasta >> ~/.sadie/germlines/sources/custom/human/IGHV.fasta
```

### Step 4: Rebuild Databases

After adding custom sequences, rebuild the databases:

```bash
# Rebuild for the specific species
sadie germlines populate -p imgt -s human --force

# This will:
# 1. Merge IMGT + OGRDB + custom sequences
# 2. Deduplicate with custom taking priority
# 3. Build BLAST databases
# 4. Create auxiliary files
```

### Step 5: Verify Integration

Check that your custom sequences were integrated:

```bash
# Count sequences in custom source
grep -c "^>" ~/.sadie/germlines/sources/custom/human/IGHV.fasta

# Count sequences in merged database
grep -c "^>" ~/.sadie/germlines/databases/human/IGHV.fasta

# The merged count should be >= custom count
```

**Inspect specific sequences:**

```bash
# Check if your custom sequence is in the final database
grep "IGHV3-30\*20" ~/.sadie/germlines/databases/human/IGHV.fasta
```

### Step 6: Test Annotation

Test that your custom sequences are used in annotation:

```bash
# Create test sequence matching your custom germline
echo ">test_seq" > test.fasta
echo "CAGGTGCAGCTGGTGCAGTCTGGGGCT..." >> test.fasta

# Annotate
sadie airr test.fasta -o test_results.tsv

# Check v_call field for your custom allele
grep "IGHV3-30\*20" test_results.tsv
```

---

## Validation

### FASTA Format Validation

**Common Issues:**

```bash
# Check for invalid characters
grep -n "[^ATCGN>]" my_custom_germlines.fasta

# Check all headers start with >
grep -n "^[^>ATCGN]" my_custom_germlines.fasta

# Check for empty lines (should be minimal)
grep -n "^$" my_custom_germlines.fasta
```

### Sequence Length

**Typical Lengths:**

- **V genes**: 250-350 bp (coding region)
- **D genes**: 10-35 bp
- **J genes**: 40-60 bp

**Check lengths:**

```bash
# Count sequence lengths
awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' IGHV.fasta
```

### Duplicate Detection

Check for duplicate sequence names:

```bash
# Find duplicate headers
grep "^>" ~/.sadie/germlines/sources/custom/human/IGHV.fasta | sort | uniq -d
```

If duplicates exist, SADIE uses the first occurrence.

---

## Advanced Usage

### Multiple Custom Sequences

**Organize by project:**

```bash
# Keep source files organized
~/my_germlines/
├── project_A/
│   └── novel_alleles.fasta
├── project_B/
│   └── strain_specific.fasta
└── merged/
    └── IGHV.fasta  # Combine all for SADIE
```

**Merge before copying:**

```bash
cat ~/my_germlines/project_A/novel_alleles.fasta \
    ~/my_germlines/project_B/strain_specific.fasta \
    > ~/my_germlines/merged/IGHV.fasta

cp ~/my_germlines/merged/IGHV.fasta ~/.sadie/germlines/sources/custom/human/
```

### Version Control

Track your custom sequences with git:

```bash
cd ~/.sadie/germlines/sources/custom
git init
git add .
git commit -m "Initial custom germlines"

# After updates
git add human/IGHV.fasta
git commit -m "Added IGHV3-30*20 novel allele"
```

### Backing Up Custom Sequences

```bash
# Backup custom sequences
tar -czf custom_germlines_backup_$(date +%Y%m%d).tar.gz \
  ~/.sadie/germlines/sources/custom/

# Restore from backup
tar -xzf custom_germlines_backup_20260122.tar.gz -C ~/
```

---

## Troubleshooting

### Sequences Not Appearing in Annotation

**Causes:**

1. Database not rebuilt after adding custom sequences
2. Sequence name doesn't match IMGT nomenclature
3. File placed in wrong directory

**Solutions:**

```bash
# Rebuild databases
sadie germlines populate -s human --force

# Verify file location
ls -lh ~/.sadie/germlines/sources/custom/human/

# Check sequence names
grep "^>" ~/.sadie/germlines/sources/custom/human/IGHV.fasta
```

### Custom Sequence Overridden

If your custom sequence has the same name and sequence as an IMGT/OGRDB sequence, custom should win. Verify:

```bash
# Compare sequences
diff <(grep -A1 "IGHV3-30\*01" ~/.sadie/germlines/sources/custom/human/IGHV.fasta) \
     <(grep -A1 "IGHV3-30\*01" ~/.sadie/germlines/databases/human/IGHV.fasta)
```

### Invalid FASTA Format

```bash
# Validate FASTA format
python -c "
from Bio import SeqIO
try:
    list(SeqIO.parse('IGHV.fasta', 'fasta'))
    print('Valid FASTA')
except:
    print('Invalid FASTA format')
"
```

---

## Best Practices

1. **Use Clear Names**: Include suffixes like `_custom`, `_novel`, `_lab` to track origin
2. **Document Source**: Keep a README with provenance information
3. **Version Control**: Use git to track changes
4. **Backup Regularly**: Custom sequences are valuable research data
5. **Test Thoroughly**: Verify custom sequences are used in annotation
6. **Update Periodically**: Check if custom alleles are added to IMGT/OGRDB

---

## See Also

- [Provider Guide](provider-guide.md) - Understanding germline providers
- [CLI Reference](cli-reference.md) - Command documentation
- [Troubleshooting](troubleshooting.md) - Common issues
