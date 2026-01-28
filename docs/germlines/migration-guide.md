# Migration Guide: G3 API to Local Germlines

## Overview

This guide helps you migrate from the legacy G3 API to SADIE's local germlines module. The migration is straightforward and brings significant benefits including offline operation, faster performance, and multi-source germline support.

!!! danger "G3 API Deprecation"
    **The G3 API will be removed after 2026-06-01.** After this date, `SADIE_USE_GERMLINES_MODULE=false` will no longer work. Please complete your migration before this deadline.

## Why Migrate?

### Benefits of Local Germlines

**Offline Operation**
Work without internet connectivity. Once populated, all germline data is available locally with no external API dependency.

**Multi-Source Support**
Combine sequences from IMGT, OGRDB, VDJbase, and your own custom germlines. Priority-based merging ensures your custom sequences take precedence.

**Faster Performance**
Eliminate network latency. Local database access is instant with no API round-trips.

**Reproducibility**
Lock germline versions for reproducible analysis. Your results won't change due to upstream database updates.

**Extended Species Coverage**
Support for 33+ species including human, mouse, rat, rabbit, rhesus macaque, pig, dog, cat, cattle, sheep, and more.

### Migration Timeline

| Date | Event |
|------|-------|
| **2026-01** | Local germlines released as default |
| **2026-01 to 2026-06** | Migration period (G3 still available) |
| **2026-06-01** | G3 API support removed |

You have approximately **5 months** to complete migration.

---

## Step-by-Step Migration

### Step 1: Verify Current Setup

Check if you're currently using G3 API:

```bash
# Check environment variable
echo $SADIE_USE_GERMLINES_MODULE

# If output is "false" or empty, you're using G3
```

### Step 2: Populate Local Germlines

Download germline databases for the species you need:

```bash
# Option A: Download all species (recommended)
sadie germlines populate

# Option B: Download specific species only
sadie germlines populate -s human -s mouse
```

**Expected time:**
- All species: 5-10 minutes
- Single species: 2-3 minutes

**Disk space required:**
- All species: ~500MB
- Human only: ~30MB

### Step 3: Verify Installation

Check that germlines were downloaded successfully:

```bash
sadie germlines status
```

Expected output:
```
┏━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━┓
┃ Provider ┃ Version       ┃ Downloaded   ┃ Species ┃ Status         ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━┩
│ imgt     │ release-202601│ 2026-01-22   │ 33      │ Up-to-date     │
│ ogrdb    │ release-202601│ 2026-01-22   │ 2       │ Up-to-date     │
└──────────┴───────────────┴──────────────┴─────────┴────────────────┘
```

### Step 4: Enable Germlines Module

Remove or update the environment variable:

```bash
# Option A: Remove the variable (germlines becomes default)
unset SADIE_USE_GERMLINES_MODULE

# Option B: Explicitly enable germlines
export SADIE_USE_GERMLINES_MODULE=true
```

To make this permanent, update your shell profile:

```bash
# For bash users (~/.bashrc or ~/.bash_profile)
echo 'export SADIE_USE_GERMLINES_MODULE=true' >> ~/.bashrc

# For zsh users (~/.zshrc)
echo 'export SADIE_USE_GERMLINES_MODULE=true' >> ~/.zshrc
```

### Step 5: Test Your Workflows

Run your existing workflows to verify everything works:

```bash
# Command-line usage
sadie airr sequences.fasta -o results.tsv

# Python API usage
python your_analysis_script.py
```

Your code requires **no changes**. SADIE automatically uses local germlines.

### Step 6: Compare Results (Optional)

If you want to verify that results match G3 output:

```bash
# Run with germlines
export SADIE_USE_GERMLINES_MODULE=true
sadie airr sequences.fasta -o results_germlines.tsv

# Run with G3 (for comparison)
export SADIE_USE_GERMLINES_MODULE=false
sadie airr sequences.fasta -o results_g3.tsv

# Compare outputs
diff results_germlines.tsv results_g3.tsv
```

!!! note "Expected Differences"
    Minor differences may occur due to updated germline sequences in newer IMGT releases or additional sequences from OGRDB. These are expected and indicate improved annotation accuracy.

---

## What Changes in Workflows

### No Code Changes Required

Your existing Python code continues to work without modification:

```python
from sadie.airr import Airr

# This works exactly the same
airr = Airr("human")
results = airr.run_file("sequences.fasta")
```

### Command-Line Usage Unchanged

CLI commands work identically:

```bash
# Before (G3)
sadie airr sequences.fasta -o output.tsv

# After (germlines) - same command!
sadie airr sequences.fasta -o output.tsv
```

### New Capabilities

**Multi-provider support:**
```python
# Germlines automatically merges IMGT + OGRDB + custom sequences
# No configuration needed - priority-based merging is automatic
```

**Offline operation:**
```bash
# Works without internet after initial setup
sadie airr sequences.fasta -o output.tsv
```

**Custom sequences:**
```bash
# Add your own novel germlines (not possible with G3)
# See Custom Sequences guide for details
```

---

## Rollback Instructions

If you need to temporarily revert to G3 API (before 2026-06-01):

### Temporary Rollback

```bash
# Set environment variable to use G3
export SADIE_USE_GERMLINES_MODULE=false

# Run SADIE commands
sadie airr sequences.fasta -o output.tsv
```

### Persistent Rollback

Update your shell profile:

```bash
# For bash
echo 'export SADIE_USE_GERMLINES_MODULE=false' >> ~/.bashrc

# For zsh
echo 'export SADIE_USE_GERMLINES_MODULE=false' >> ~/.zshrc
```

!!! warning "Rollback Only Until 2026-06-01"
    G3 API support will be completely removed after 2026-06-01. Rollback is only available during the migration period.

---

## Frequently Asked Questions

### Do I need to change my code?

No. SADIE automatically detects and uses local germlines. Your existing Python scripts and command-line workflows continue to work without modification.

### Can I use both G3 and germlines?

Yes, during the migration period (until 2026-06-01) you can switch between them using the `SADIE_USE_GERMLINES_MODULE` environment variable.

### What if I only need human data?

Download only the species you need:

```bash
sadie germlines populate -s human
```

This reduces download time and disk usage.

### Will my results change?

Results may differ slightly due to:
- **Updated sequences**: IMGT databases are updated monthly with new alleles
- **Additional sources**: OGRDB provides validated novel alleles not in IMGT
- **Better coverage**: More comprehensive germline coverage improves annotation accuracy

These differences indicate improved annotation quality.

### How do I keep germlines up-to-date?

Check for updates monthly:

```bash
# Check status
sadie germlines status

# If updates available, re-populate
sadie germlines populate --force
```

### What happens if download fails?

The populate command automatically resumes from the last successful species:

```bash
# If download fails partway through
sadie germlines populate  # Automatically resumes from checkpoint
```

### Can I use custom germline sequences?

Yes! This is a key advantage over G3. See the [Custom Sequences Guide](custom-sequences.md) for details.

### What if I need a species not available?

IMGT supports 33+ species. If your species is available in IMGT, it's automatically included. For custom species, use the custom sequences feature.

### Do I need internet access after setup?

No. After the initial `sadie germlines populate`, SADIE works completely offline.

### How much disk space is needed?

- **All species**: ~500MB
- **Human only**: ~30MB
- **Human + mouse**: ~60MB

---

## Troubleshooting

**"Species not found" error:**
```bash
# Solution: Populate germlines first
sadie germlines populate -s human
```

**Download hangs:**
- Check internet connection
- Some species take longer (normal)
- Use `--dry-run` to preview what will be downloaded

**Results differ from G3:**
- Expected due to updated germline sequences
- Differences indicate improved annotation accuracy
- Compare using `diff` to verify changes

**Need more help?**
- See [Troubleshooting Guide](troubleshooting.md)
- Report issues: [GitHub Issues](https://github.com/jwillis0720/sadie/issues)

---

## Next Steps

After migration:

1. ✅ Remove `SADIE_USE_GERMLINES_MODULE=false` from shell profiles
2. ✅ Set up monthly germline updates
3. ✅ Explore multi-source capabilities (IMGT + OGRDB)
4. ✅ Consider adding custom sequences if needed

**Need more information?**
- [CLI Reference](cli-reference.md) - Complete command documentation
- [Provider Guide](provider-guide.md) - Understanding IMGT vs OGRDB
- [Custom Sequences](custom-sequences.md) - Adding your own germlines
