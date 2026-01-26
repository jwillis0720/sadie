# Provider Guide: Understanding Germline Data Sources

## Overview

SADIE's germlines module supports multiple germline sequence providers, each with unique strengths and coverage. Understanding the differences between providers helps you choose the right data sources for your research.

**Available Providers:**

- **IMGT** - International ImMunoGeneTics information system
- **OGRDB** - Open Germline Receptor Database
- **VDJbase** - Curated population-specific germline database
- **Custom** - Your own novel germline sequences

All providers can be used simultaneously. SADIE automatically merges sequences with priority-based deduplication to give you the most comprehensive germline coverage.

---

## Provider Comparison

### IMGT (International ImMunoGeneTics)

**What is IMGT?**

IMGT is the international reference in immunogenetics and immunoinformatics, established in 1989. It provides the most comprehensive collection of germline sequences with standardized nomenclature.

**Key Features:**

- **Species Coverage**: 33+ species including human, mouse, rat, rabbit, rhesus macaque, pig, dog, cat, cattle, sheep, horse, and more
- **Gene Types**: V, D, J segments for all IG and TR loci
- **Update Frequency**: Monthly releases (release-YYYYMM format)
- **Nomenclature**: Official IMGT nomenclature standard
- **Validation**: Expert-curated with peer review

**When to Use IMGT:**

- ✅ Working with any of the 33+ supported species
- ✅ Need official nomenclature for publications
- ✅ Require comprehensive V/D/J coverage
- ✅ Standard research workflows

**Limitations:**

- Contains primarily common/consensus alleles
- May lack recently discovered novel alleles
- Conservative curation (new alleles take time to include)

**Example:**

```bash
# Download IMGT data only
sadie germlines populate -p imgt

# Download specific species from IMGT
sadie germlines populate -p imgt -s human -s mouse
```

---

### OGRDB (Open Germline Receptor Database)

**What is OGRDB?**

OGRDB is a community-driven database for novel germline alleles, launched in 2020. It provides a platform for researchers to submit and validate newly discovered germline sequences.

**Key Features:**

- **Species Coverage**: Currently human and mouse
- **Novel Alleles**: Contains validated novel alleles not yet in IMGT
- **Validation Levels**: Multiple validation tiers (inferred, affirmed, etc.)
- **Community-Driven**: Researchers can submit discoveries
- **Update Frequency**: Regular updates as new sequences are validated

**When to Use OGRDB:**

- ✅ Working with human or mouse samples
- ✅ Analyzing diverse populations with rare alleles
- ✅ Need the most up-to-date novel allele discoveries
- ✅ Working with non-European populations (better coverage)

**Limitations:**

- Limited to human and mouse currently
- Smaller total number of sequences compared to IMGT
- Some sequences may have different nomenclature

**Example:**

```bash
# Download OGRDB data only
sadie germlines populate -p ogrdb

# OGRDB only has human and mouse
# Species filter not needed
```

**Novel Allele Example:**

OGRDB might contain `IGHV3-30*19` which represents a validated novel allele discovered through high-throughput sequencing studies but not yet included in IMGT's official database.

---

### VDJbase (Population Germline Database)

**What is VDJbase?**

VDJbase is a curated database of population-specific germline sequences derived from AIRR-seq datasets. It provides haplotype information and population-specific allele frequencies.

**Key Features:**

- **Population-Specific**: Sequences organized by population/ethnicity
- **Haplotypes**: Contains haplotype information
- **Frequency Data**: Allele frequency information
- **Derived from AIRR-seq**: Based on high-throughput repertoire data

**When to Use VDJbase:**

- ✅ Population genetics studies
- ✅ Rare allele discovery in specific populations
- ✅ Need allele frequency information
- ✅ Haplotype-aware analysis

**Limitations:**

- Limited species coverage
- Availability depends on population studies
- Not all populations equally represented

**Status in SADIE:**

```bash
# VDJbase support is configured but requires data availability
sadie germlines status
# Shows: vdjbase │ - │ - │ 0 │ Not downloaded
```

VDJbase integration is planned for future releases as the database expands.

---

### Custom Sequences

**What are Custom Sequences?**

Custom sequences are your own novel germline alleles that you've discovered and want to include in annotation. These might be unpublished novel alleles specific to your research.

**Key Features:**

- **Highest Priority**: Custom sequences override provider sequences
- **Flexible Format**: Simple FASTA format
- **Species-Specific**: Organized by species
- **Complete Control**: You manage the sequences

**When to Use Custom Sequences:**

- ✅ Discovered novel alleles not in public databases
- ✅ Working with model organisms or rare species
- ✅ Need to test specific sequence variants
- ✅ Proprietary germline data

See [Custom Sequences Guide](custom-sequences.md) for detailed instructions.

---

## Priority-Based Merging

### How It Works

When multiple providers have sequences with the same name, SADIE uses priority-based merging to resolve conflicts:

**Priority Order (Highest to Lowest):**

1. **VDJbase** sequences (curated, validated alleles from population studies)
2. **OGRDB** sequences (community-curated novel alleles)
3. **IMGT** sequences (comprehensive reference database)
4. **Custom** sequences (internal lab sequences for edge cases)

**Example Scenario:**

Suppose VDJbase, OGRDB, IMGT, and your custom sequences all have `IGHV3-30*01`:

```
VDJbase: IGHV3-30*01  → CAGGTGCAGCT... (500bp)
OGRDB:   IGHV3-30*01  → CAGGTGCAGCT... (500bp)
IMGT:    IGHV3-30*01  → CAGGTGCAGCT... (495bp)
Custom:  IGHV3-30*01  → CAGGTGCAGCT... (500bp)
```

**Result:** SADIE uses the **VDJbase** sequence because it has the highest priority.

### Deduplication Strategy

**Exact Match Deduplication:**

```python
# Sequences are considered duplicates if:
# 1. Same gene name (e.g., "IGHV3-30*01")
# 2. Same sequence (after normalization)

# Lower priority sequences are dropped
```

**No Partial Matching:**

SADIE does **not** merge partial sequences or attempt sequence alignment. Each sequence is treated as an independent entry.

### Benefits

**Comprehensive Coverage:**

- Get curated alleles from VDJbase (validated population data)
- Plus novel alleles from OGRDB (recent discoveries)
- Plus comprehensive reference from IMGT (species diversity)
- Plus your custom sequences (unpublished findings)

**Automatic Conflict Resolution:**

- No manual intervention needed
- Deterministic priority rules
- VDJbase curated data takes precedence

**Reproducibility:**

- Same input providers → same output database
- Version-locked for reproducible analysis

---

## Configuring Providers

### Download Specific Providers

```bash
# IMGT only (default for most users)
sadie germlines populate -p imgt

# OGRDB only (human and mouse novel alleles)
sadie germlines populate -p ogrdb

# All available providers (recommended)
sadie germlines populate -p all
```

### Combining Providers

**Recommended Setup for Human/Mouse:**

```bash
# Download both IMGT and OGRDB for maximum coverage
sadie germlines populate -p imgt -s human -s mouse
sadie germlines populate -p ogrdb

# Verify both are installed
sadie germlines status
```

Expected output:
```
┏━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━┓
┃ Provider ┃ Version       ┃ Downloaded   ┃ Species ┃ Status         ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━┩
│ imgt     │ release-202601│ 2026-01-22   │ 2       │ Up-to-date     │
│ ogrdb    │ release-202601│ 2026-01-22   │ 2       │ Up-to-date     │
└──────────┴───────────────┴──────────────┴─────────┴────────────────┘
```

**For Other Species:**

```bash
# IMGT is the only option for non-human/mouse species
sadie germlines populate -p imgt -s rabbit -s rhesus_macaque
```

---

## Verification and Inspection

### Check Downloaded Providers

```bash
sadie germlines status
```

### Inspect Merged Database

The merged database is stored in:

```
~/.sadie/germlines/
├── sources/              # Raw provider data
│   ├── imgt/
│   │   └── human/
│   │       ├── IGHV.fasta
│   │       ├── IGHD.fasta
│   │       └── IGHJ.fasta
│   ├── ogrdb/
│   │   └── human/
│   │       └── IGHV.fasta
│   └── custom/
│       └── human/
│           └── IGHV.fasta
└── databases/            # Merged BLAST databases
    └── human/
        ├── IGHV.ndb
        ├── IGHV.aux
        └── IGHV.ndm.imgt
```

### Count Sequences by Provider

```bash
# Count IMGT sequences
grep -c "^>" ~/.sadie/germlines/sources/imgt/human/IGHV.fasta

# Count OGRDB sequences
grep -c "^>" ~/.sadie/germlines/sources/ogrdb/human/IGHV.fasta

# Count merged sequences (after deduplication)
grep -c "^>" ~/.sadie/germlines/databases/human/IGHV.fasta
```

---

## Best Practices

### For Most Users

**Recommended:**
```bash
# Download all providers for comprehensive coverage
sadie germlines populate
```

This gives you:
- IMGT for comprehensive baseline
- OGRDB for novel human/mouse alleles
- Option to add custom sequences later

### For Human/Mouse Research

**Recommended:**
```bash
# Maximize coverage with IMGT + OGRDB
sadie germlines populate -p imgt -s human -s mouse
sadie germlines populate -p ogrdb
```

### For Other Species

**Required:**
```bash
# IMGT is the only option
sadie germlines populate -p imgt -s <species>
```

### For Custom Species/Alleles

**Required:**
```bash
# Start with IMGT baseline
sadie germlines populate -p imgt -s human

# Add your custom sequences
# See Custom Sequences Guide
```

---

## Updating Provider Data

Providers release updates monthly (IMGT, OGRDB) or as needed (VDJbase, Custom).

**Check for Updates:**

```bash
sadie germlines status
# Look for "Update available" status
```

**Update All Providers:**

```bash
sadie germlines populate --force
```

**Update Specific Provider:**

```bash
sadie germlines populate -p imgt --force
```

!!! note "Update Frequency"
    - **IMGT**: Monthly releases (around the 1st of each month)
    - **OGRDB**: Irregular updates as sequences are validated
    - **Custom**: You control when to update

    We recommend checking monthly and updating quarterly unless working with rapidly evolving allele databases.

---

## Troubleshooting

**Provider download fails:**
```bash
# Check internet connection
# Try downloading one provider at a time
sadie germlines populate -p imgt
sadie germlines populate -p ogrdb
```

**Sequences missing after merge:**
- Expected behavior if sequences are duplicates
- Higher priority provider sequences take precedence
- Use inspection commands above to verify

**OGRDB data not found:**
- OGRDB only has human and mouse
- Ensure you're requesting available species

**Need more help?**
- [Troubleshooting Guide](troubleshooting.md)
- [GitHub Issues](https://github.com/jwillis0720/sadie/issues)

---

## See Also

- [CLI Reference](cli-reference.md) - Command documentation
- [Custom Sequences](custom-sequences.md) - Adding your own sequences
- [Migration Guide](migration-guide.md) - Moving from G3 API
