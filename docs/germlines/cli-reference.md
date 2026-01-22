# CLI Reference: Germline Commands

Complete command-line reference for SADIE germline database management.

## sadie germlines populate

Download germline data from providers and build databases.

### Synopsis

```bash
sadie germlines populate [OPTIONS]
```

### Description

Downloads germline sequences from IMGT, OGRDB, and/or VDJbase providers and builds the necessary databases for AIRR annotation and renumbering.

The command performs:
1. Downloads germline FASTA files from selected providers
2. Normalizes and merges sequences with priority-based deduplication
3. Builds BLAST databases for IgBLAST
4. Creates auxiliary files (*.aux) for CDR3 reconstruction
5. Generates internal_data files (*.ndm.imgt)
6. Validates downloaded data

### Options

#### `-p, --provider [imgt|ogrdb|vdjbase|all]`

Provider to download from. Default: `all`

**Examples:**

```bash
# Download from all providers (default)
sadie germlines populate

# Download from IMGT only
sadie germlines populate -p imgt

# Download from OGRDB only
sadie germlines populate -p ogrdb
```

#### `-s, --species SPECIES`

Specific species to download. Can be repeated for multiple species.

**Examples:**

```bash
# Download human only
sadie germlines populate -s human

# Download human and mouse
sadie germlines populate -s human -s mouse

# Download all species (default - don't specify -s)
sadie germlines populate
```

**Available Species** (IMGT):
`human`, `mouse`, `rat`, `rabbit`, `rhesus_macaque`, `pig`, `dog`, `cat`, `cattle`, `sheep`, `horse`, `chicken`, `duck`, `turkey`, `zebrafish`, and 18 more.

#### `-f, --force`

Force re-download even if data is up-to-date.

**Example:**

```bash
# Re-download all data
sadie germlines populate --force

# Re-download specific species
sadie germlines populate -p imgt -s human --force
```

!!! warning "Force Re-download"
    The `--force` flag will delete existing data and re-download everything for the selected provider/species. Use with caution.

#### `--dry-run`

Show what would be downloaded without actually downloading.

**Example:**

```bash
# Preview what will be downloaded
sadie germlines populate --dry-run

# Preview specific provider
sadie germlines populate -p imgt -s human --dry-run
```

**Output:**
```
SADIE Germline Database Population
==================================================
DRY RUN - no changes will be made

imgt: Would download 1 species
  - human
```

### Examples

**Download all providers and species:**

```bash
sadie germlines populate
```

Expected time: 5-10 minutes
Downloaded: ~500MB
Species: 33

---

**Download only human from IMGT:**

```bash
sadie germlines populate -p imgt -s human
```

Expected time: 2-3 minutes
Downloaded: ~30MB
Species: 1

---

**Download human and mouse from all providers:**

```bash
sadie germlines populate -s human -s mouse
```

Expected time: 3-5 minutes
Downloaded: ~60MB
Species: 2

---

**Check what would be downloaded:**

```bash
sadie germlines populate --dry-run
```

Shows provider status and species counts without downloading.

---

**Force re-download of IMGT data:**

```bash
sadie germlines populate -p imgt --force
```

Deletes existing IMGT data and re-downloads.

### Output Explanation

The command displays:

1. **Provider Status Table** - Shows current versions and update status
2. **Download Progress** - Real-time progress with species counts
3. **Post-Download Build** - BLAST database and auxiliary file generation
4. **Validation** - Data integrity checks
5. **Summary Table** - Success/failure counts by provider

**Example Output:**

```
SADIE Germline Database Population
==================================================

Provider Status
┏━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ Provider ┃ Status        ┃ Species ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ imgt     │ Needs update  │ 0       │
│ ogrdb    │ Needs update  │ 0       │
│ vdjbase  │ Needs update  │ 0       │
└──────────┴───────────────┴─────────┘

Processing imgt...
imgt: human ████████████████████ 100%

Running post-download build pipeline...
Building auxiliary files... ✓
Building internal_data... ✓
Build complete

Validating imgt data...
imgt data validation passed

Summary
==================================================
┏━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━┓
┃ Provider ┃ Status  ┃ Downloaded ┃ Failed ┃
┡━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━┩
│ imgt     │ success │ 33         │ 0      │
│ ogrdb    │ success │ 2          │ 0      │
│ vdjbase  │ skipped │ 0          │ 0      │
└──────────┴─────────┴────────────┴────────┘
```

### Post-Download Process

After downloading, the command automatically:

1. **Builds BLAST databases** for IgBLAST using `makeblastdb`
2. **Creates auxiliary files** (*.aux) containing CDR/FWR annotations
3. **Generates internal_data** files (*.ndm.imgt) for IgBLAST
4. **Validates data** integrity and completeness

These steps require no user intervention.

### Resume Failed Downloads

If a download fails partway through, the command automatically resumes from the last successful species:

```bash
# First attempt fails at species 10/33
sadie germlines populate

# Re-run resumes from species 11
sadie germlines populate
```

Checkpoint files (`.populate_checkpoint.json`) track progress.

---

## sadie germlines status

Show status of local germline databases.

### Synopsis

```bash
sadie germlines status
```

### Description

Displays the current state of germline databases including:
- Downloaded providers and versions
- Number of species per provider
- Update availability status
- Last download timestamps

### Examples

```bash
sadie germlines status
```

### Output

```
Germline Database Status
┏━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━┓
┃ Provider ┃ Version       ┃ Downloaded   ┃ Species ┃ Status         ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━┩
│ imgt     │ release-202601│ 2026-01-21   │ 33      │ Up-to-date     │
│ ogrdb    │ release-202601│ 2026-01-20   │ 2       │ Up-to-date     │
│ vdjbase  │ -             │ -            │ 0       │ Not downloaded │
└──────────┴───────────────┴──────────────┴─────────┴────────────────┘
```

### Status Interpretations

- **Up-to-date** - Data is current for this month's release
- **Update available** - Newer data is available (re-run populate)
- **Not downloaded** - Provider has never been downloaded

### Version Format

Versions follow the format `release-YYYYMM` (e.g., `release-202601` = January 2026).

---

## Environment Variables

### SADIE_USE_GERMLINES_MODULE

Controls whether to use the germlines module or legacy G3 paths.

**Type**: Boolean
**Default**: `true`
**Values**: `true`, `false`, `1`, `0`, `yes`, `no`

#### Usage

**Use germlines module (default):**

```bash
export SADIE_USE_GERMLINES_MODULE=true
sadie airr input.fasta -o output.tsv
```

**Use legacy G3 paths (deprecated):**

```bash
export SADIE_USE_GERMLINES_MODULE=false
sadie airr input.fasta -o output.tsv
```

!!! danger "G3 Deprecation"
    G3 API support will be **removed after 2026-06-01**. The `SADIE_USE_GERMLINES_MODULE=false` option will stop working. [Migrate to germlines](migration-guide.md) before this deadline.

#### When to Use

**Set to `false` only if:**
1. You need to temporarily use legacy G3 data
2. You're comparing germlines output with G3 output
3. You're troubleshooting migration issues

**Default behavior (no environment variable):**
- Uses germlines module
- No action needed

---

## Common Workflows

### First-Time Setup

```bash
# 1. Populate databases
sadie germlines populate

# 2. Verify installation
sadie germlines status

# 3. Use SADIE normally
sadie airr sequences.fasta -o results.tsv
```

### Update Germline Data Monthly

```bash
# Check for updates
sadie germlines status

# If updates available, re-populate
sadie germlines populate --force
```

### Species-Specific Setup

```bash
# Download only the species you need
sadie germlines populate -p imgt -s human -s mouse

# Verify
sadie germlines status

# Add more species later
sadie germlines populate -p imgt -s rabbit
```

### Test Drive Before Full Download

```bash
# Check what will be downloaded
sadie germlines populate --dry-run

# Download single species to test
sadie germlines populate -p imgt -s human

# If satisfied, download all
sadie germlines populate
```

---

## Troubleshooting

**Command hangs during download:**
- Check internet connection
- Some species take longer (normal)
- Use `--dry-run` first to see what's being downloaded

**"Species not found" error:**
- Run `sadie germlines populate -s <species>` first
- Check spelling of species name
- See available species list above

**Download fails:**
- Re-run command - it will resume from checkpoint
- Check available disk space (~500MB needed for all species)
- See [Troubleshooting Guide](troubleshooting.md) for detailed solutions

---

## See Also

- [Migration Guide](migration-guide.md) - Migrate from G3 API
- [Provider Guide](provider-guide.md) - Understanding providers
- [Troubleshooting](troubleshooting.md) - Common issues
- [Overview](index.md) - Getting started
