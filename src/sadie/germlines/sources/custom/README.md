# Custom Germline Sequences

**USER EDITABLE AREA** - Add your own germline sequences here!

This is where you add novel germline sequences to use alongside IMGT and OGRDB databases.

## Quick Start

1. **Create species directory** (e.g., `human/`, `mouse/`)
2. **Add FASTA files** with naming: `IG{chain}{segment}.fasta`
3. **Run Sadie** - custom sequences are automatically processed!

## Directory Structure

```
custom/
├── README.md (this file)
├── _template/
│   └── example_IGHV.fasta  (example file)
├── human/
│   ├── IGHV.fasta         ← Add your sequences here
│   ├── IGHD.fasta
│   ├── IGHJ.fasta
│   ├── IGKV.fasta
│   ├── IGKJ.fasta
│   ├── IGLV.fasta
│   ├── IGLJ.fasta
│   └── .processed/        (auto-generated - don't edit)
│       ├── metadata.json
│       └── last_modified.txt
└── mouse/
    └── ... (same structure)
```

## File Naming Convention

**Required Pattern:** `IG{chain}{segment}.fasta`

Examples:
- Heavy chain V genes: `IGHV.fasta`
- Heavy chain D genes: `IGHD.fasta`
- Heavy chain J genes: `IGHJ.fasta`
- Kappa V genes: `IGKV.fasta`
- Kappa J genes: `IGKJ.fasta`
- Lambda V genes: `IGLV.fasta`
- Lambda J genes: `IGLJ.fasta`

## FASTA Format

### Option 1: Ungapped Sequences (Recommended)

**Easiest for users** - just provide raw nucleotide sequences:

```fasta
>IGHV1-NL1*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAG
GTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGAC

>IGHV3-NOVEL*01 My novel sequence from patient cohort
GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGAC
TCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCCATGAGCTGGGTCCGCCA
```

**System automatically:**
- ✅ Aligns to IMGT reference
- ✅ Generates gapped version
- ✅ Builds IgBLAST databases
- ✅ Creates auxiliary files

### Option 2: Pre-Gapped Sequences (Advanced)

If you already have IMGT-numbered sequences:

```fasta
>IGHV1-NL1*01
cag.gtgcagctggtgcag...tctggggctgag...gtgaag...aagcctggggcc
```

Use dots (`.`) for IMGT gaps. System will detect and use directly.

## Gene Naming

**Flexible naming supported:**

### Standard IMGT Format (Recommended)
```
>IGHV1-69*01
>IGKV3-20*01
```

### Custom Names
```
>my_novel_sequence
→ Auto-formatted to: IGHV-CUSTOM-my_novel_sequence*01
```

### With Descriptions
```
>IGHV1-69*01|found in patient XYZ cohort
→ Description ignored, gene name extracted
```

### Missing Chain Designation
```
>IGV1-2*01
→ Auto-corrected based on file (IGHV.fasta → IGHV1-2*01)
```

## Priority System

**Custom sequences fill gaps from standard databases.**

Default priority order: `vdjbase > ogrdb > imgt > custom`

Custom sequences are applied LAST, after validated databases. This ensures:
- **Validated data first**: VDJbase and OGRDB provide curated, community-validated alleles
- **Gap filling**: Custom sequences add lab-specific alleles not in public databases
- **Override capability**: To override a standard allele, use the same gene name in custom

To explicitly override a VDJbase/OGRDB gene, add a custom sequence with the identical gene name.

### Examples:

**Override IMGT gene:**
```bash
# Your custom/human/IGHV.fasta contains:
>IGHV1-69*01 Corrected version for our population
GAGGTGCAGCTGGTGCAGTCT...

# Result: Your version used instead of IMGT
```

**Add novel gene:**
```bash
# Your custom/human/IGHV.fasta contains:
>IGHV1-NEW*01 Novel allele discovered in our cohort
CAGGTGCAGCTGGTGCAGTCT...

# Result: Added to database alongside IMGT genes
```

**Mix all sources:**
```bash
# custom: 3 genes
# IMGT: 400 genes
# OGRDB: 50 genes

# Result: ~450 total genes (deduplicated by priority)
```

## Workflow Examples

### Example 1: Add Single Novel Sequence

```bash
# 1. Create file
$ cat > custom/human/IGHV.fasta << 'EOF'
>IGHV1-NOVEL*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCT
EOF

# 2. Run Sadie
$ sadie airr -i sequences.fasta -o results.tsv

# Output:
# INFO: Detected new custom file: IGHV.fasta
# INFO: Processing 1 sequence
# INFO: Aligning IGHV1-NOVEL*01...
# INFO: IgBLAST found IGHV1-NOVEL*01 in 15 sequences!
```

### Example 2: Override IMGT Gene

```bash
# Use corrected version of IGHV1-69*01
$ cat > custom/human/IGHV.fasta << 'EOF'
>IGHV1-69*01 Corrected for our population
GAGGTGCAGCTGGTGCAGTCT...
EOF

# Your version now used instead of IMGT version
```

### Example 3: Add Multiple Custom Genes

```bash
# Add several novel genes to one file
$ cat > custom/human/IGHV.fasta << 'EOF'
>IGHV1-NL1*01 Novel allele 1
CAGGTGCAGCTGGTGCAGTCTGGG...

>IGHV1-NL2*01 Novel allele 2
CAGGTGCAGCTGGTGGAGTCTGGG...

>IGHV3-COHORT*01 From patient cohort
GAGGTGCAGCTGGTGGAGTCTGGG...
EOF
```

## Automatic Processing

System automatically handles:

✅ **Change Detection**: Modified files trigger rebuild
✅ **Sequence Validation**: Checks for valid nucleotides
✅ **Auto-Gapping**: Aligns ungapped to IMGT reference
✅ **Database Building**: Creates BLAST databases
✅ **Auxiliary Files**: Generates IgBLAST aux files
✅ **Metadata Tracking**: Records processing info

## Troubleshooting

### "Could not gap sequence"

**Cause**: Sequence too divergent from IMGT reference

**Solutions**:
1. Provide pre-gapped sequence with IMGT numbering
2. Check sequence is in correct reading frame
3. Verify sequence is from correct species/chain

### "No sequences found"

**Cause**: File naming or format issue

**Check**:
- File name matches pattern: `IG{H|K|L}{V|D|J}.fasta`
- FASTA headers start with `>`
- File is in correct species directory
- No empty lines before first sequence

### "Invalid nucleotide"

**Cause**: Non-ACGT characters in sequence

**Solution**:
- Use only A, C, G, T, N
- Remove spaces, numbers, or other characters
- For gapped sequences, use only dots (`.`)

## Validation

Check your custom sequences loaded correctly:

```python
from sadie.germlines import get_manager

manager = get_manager()

# Get all genes (including custom)
genes = manager.get_genes("human", "V", "H")

# Filter to custom only
custom_genes = [g for g in genes if g.source == "custom"]

print(f"Custom genes loaded: {len(custom_genes)}")
for gene in custom_genes:
    print(f"  - {gene.name}")
```

## Advanced: Pre-Gapped Sequences

If you have sequences with IMGT numbering:

### IMGT Numbering Scheme

```
Position:  1-26   27-38  39-55  56-65  66-104 105-117 118-128
Region:    FWR1   CDR1   FWR2   CDR2   FWR3   CDR3    FWR4
```

### Example Gapped Sequence

```fasta
>IGHV1-69*01
cag.gtgcagctggtgcag...tctggggctgag...gtgaag...aagcctggggcc
tcagtgaag...gtctcctgcaaggcttctggt...tacacctttacc...agctatggtatc
agc...tgggtgcgac...aggcccctgga...caagggcttgagtgg...atgggat
ggat...cagcgcttacaatggt...aacacaaactatgcacag...aagctc
```

Dots (`.`) indicate IMGT gaps at standardized positions.

## More Information

- **Documentation**: https://sadie.jordanrwillis.com/germlines/custom
- **IMGT Numbering**: https://www.imgt.org/IMGTScientificChart/Numbering/
- **Examples**: See `_template/example_IGHV.fasta`

## Support

Questions? Issues?
- GitHub Issues: https://github.com/jwillis0720/sadie/issues
- Email: jwillis0720@gmail.com
