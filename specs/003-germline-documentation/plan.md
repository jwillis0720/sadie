# Documentation Plan: Germline Database Integration

**Branch**: `sadie_docs` | **Date**: 2026-01-22 | **Based on**: [002-germline-integration](../002-germline-integration/plan.md)

## Summary

Create comprehensive user-facing documentation for the germline database integration feature, including CLI command reference, migration guides, API documentation, and troubleshooting resources. This documentation will help users transition from G3 API to the new local germlines module and understand the new `sadie germlines` commands.

## Documentation Scope

### What Changed (User-Facing)

Based on the 002-germline-integration implementation:

1. **New CLI Commands** (app.py:612-701)
   - `sadie germlines populate` - Download and build germline databases
   - `sadie germlines status` - Check database status and versions

2. **Feature Flag** (airr/igblast/germline.py:14-23)
   - `SADIE_USE_GERMLINES_MODULE` environment variable (default: true)
   - Controls whether to use new germlines module or legacy G3 paths

3. **Multi-Species Support**
   - Expanded from human-only to 33+ species
   - Support for IMGT, OGRDB, VDJbase providers

4. **Offline Operation**
   - No longer requires G3 API connectivity
   - All data stored locally after initial setup

5. **Custom Germline Support**
   - Users can add custom sequences in `sources/custom/`
   - Priority-based merging (custom > ogrdb > vdjbase > imgt)

### What Didn't Change (User Experience)

- Existing `sadie airr` and `sadie renumbering` commands work unchanged
- Output formats remain the same
- API interfaces remain backward compatible
- Existing workflows continue to work (G3 still available as fallback)

## Documentation Structure

```text
docs/
├── germlines/                          # NEW: Main germlines documentation
│   ├── index.md                       # Overview and quick start
│   ├── cli-reference.md               # Complete CLI command reference
│   ├── migration-guide.md             # G3 → Germlines migration
│   ├── provider-guide.md              # Working with IMGT/OGRDB/VDJbase
│   ├── custom-sequences.md            # Adding custom germlines
│   ├── troubleshooting.md             # Common issues and solutions
│   └── architecture.md                # Technical architecture (for developers)
│
├── index.md                           # UPDATED: Add germlines section
├── annotation.md                      # UPDATED: Add germlines backend info
└── reference.md                       # UPDATED: Mention germlines option
```

## Documentation Deliverables

### 1. Germlines Overview (docs/germlines/index.md)

**Priority**: P0 (Critical)
**Target Audience**: All users
**Estimated Length**: 500-700 words

**Content**:
- What is the germlines module?
- Why use it over G3? (offline, multi-source, custom sequences)
- Quick start example
- Link to detailed guides

**Example Structure**:
```markdown
# Germline Database Management

## Overview
The germlines module provides local germline database management...

## Why Use Local Germlines?
- ✅ **Offline operation**: No internet required after setup
- ✅ **Multi-source**: IMGT, OGRDB, VDJbase, custom
- ✅ **Faster**: No API latency
- ✅ **Reproducible**: Fixed versions

## Quick Start
\```bash
# Download germline data (one-time setup)
sadie germlines populate

# Check status
sadie germlines status

# Use with AIRR annotation (automatic)
sadie airr input.fasta -o output.tsv
\```

## Next Steps
- [CLI Reference](cli-reference.md) - Complete command documentation
- [Migration Guide](migration-guide.md) - Migrate from G3
- [Custom Sequences](custom-sequences.md) - Add your own germlines
```

### 2. CLI Reference (docs/germlines/cli-reference.md)

**Priority**: P0 (Critical)
**Target Audience**: All users
**Estimated Length**: 1000-1500 words

**Content**:
- Complete reference for `sadie germlines populate`
  - All options: --provider, --species, --force, --dry-run
  - Examples for common use cases
  - Output explanation
- Complete reference for `sadie germlines status`
  - Reading the status table
  - Understanding version information
- Integration with existing commands
  - Using SADIE_USE_GERMLINES_MODULE flag
  - Specifying database paths

**Example Structure**:
```markdown
# CLI Reference: Germline Commands

## sadie germlines populate

Download germline data from providers and build databases.

### Synopsis
\```bash
sadie germlines populate [OPTIONS]
\```

### Options

**-p, --provider [imgt|ogrdb|vdjbase|all]**
- Provider to download from (default: all)
- Examples:
  \```bash
  sadie germlines populate -p imgt      # IMGT only
  sadie germlines populate -p all       # All providers
  \```

**-s, --species SPECIES**
- Specific species to download (repeatable)
- Examples:
  \```bash
  sadie germlines populate -s human
  sadie germlines populate -s human -s mouse
  \```

**-f, --force**
- Force re-download even if up-to-date

**--dry-run**
- Show what would be downloaded without downloading

### Examples

**Download all providers and species:**
\```bash
sadie germlines populate
\```

**Download only human from IMGT:**
\```bash
sadie germlines populate -p imgt -s human
\```

**Check what would be downloaded:**
\```bash
sadie germlines populate --dry-run
\```

**Force re-download:**
\```bash
sadie germlines populate --force
\```

### Output Explanation

The command displays:
1. Provider status table (current vs. available versions)
2. Download progress with species counts
3. Post-download build pipeline progress
4. Summary table with success/failure counts

### Post-Download Process

After downloading, the command automatically:
1. Builds BLAST databases for IgBLAST
2. Creates auxiliary files (*.aux) for CDR3 reconstruction
3. Generates internal_data files (*.ndm.imgt)
4. Validates downloaded data

---

## sadie germlines status

Show status of local germline databases.

### Synopsis
\```bash
sadie germlines status
\```

### Output

Displays table with:
- Provider name
- Current version (release-YYYYMM)
- Download date
- Number of species
- Update status

### Example Output

\```
┏━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━┓
┃ Provider ┃ Version       ┃ Downloaded   ┃ Species ┃ Status         ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━┩
│ imgt     │ release-202601│ 2026-01-21   │ 33      │ Up-to-date     │
│ ogrdb    │ release-202601│ 2026-01-20   │ 2       │ Up-to-date     │
│ vdjbase  │ -             │ -            │ 0       │ Not downloaded │
└──────────┴───────────────┴──────────────┴─────────┴────────────────┘
\```

---

## Environment Variables

### SADIE_USE_GERMLINES_MODULE

Controls whether to use germlines module or legacy G3 paths.

**Values**: `true` (default), `false`, `1`, `0`, `yes`, `no`

**Usage**:
\```bash
# Use germlines module (default)
export SADIE_USE_GERMLINES_MODULE=true
sadie airr input.fasta

# Use legacy G3 paths
export SADIE_USE_GERMLINES_MODULE=false
sadie airr input.fasta
\```

**Deprecation Notice**: G3 API support will be removed after 2026-06-01.
```

### 3. Migration Guide (docs/germlines/migration-guide.md)

**Priority**: P0 (Critical)
**Target Audience**: Existing G3 users
**Estimated Length**: 800-1000 words

**Content**:
- Why migrate from G3?
- Step-by-step migration instructions
- What changes in your workflows
- Rollback instructions if needed
- FAQ about migration

**Example Structure**:
```markdown
# Migration Guide: G3 to Germlines Module

## Why Migrate?

The G3 API is being deprecated and will be removed after **2026-06-01**.

Benefits of migrating:
- ✅ **Offline operation**: No internet dependency
- ✅ **Faster**: No API latency
- ✅ **More sources**: IMGT + OGRDB + VDJbase + custom
- ✅ **Reproducible**: Fixed data versions

## Migration Timeline

- **2026-01**: Germlines module released (default)
- **2026-06-01**: G3 API support removed

## Step-by-Step Migration

### Step 1: Check Current Version

\```bash
sadie --version
# Ensure you have version >= 0.10.0
\```

### Step 2: Download Germline Data

\```bash
# Download all providers (recommended)
sadie germlines populate

# Or download specific providers
sadie germlines populate -p imgt -s human
\```

This takes 2-5 minutes depending on species count.

### Step 3: Verify Installation

\```bash
sadie germlines status
\```

Check that your species shows as "Up-to-date".

### Step 4: Test Your Workflow

\```bash
# Your existing command should work unchanged
sadie airr input.fasta -o output.tsv
\```

The germlines module is now active by default.

### Step 5: Update Scripts (Optional)

No code changes required! The germlines module is backward compatible.

Optional: Set explicit environment variable in scripts:
\```bash
#!/bin/bash
export SADIE_USE_GERMLINES_MODULE=true
sadie airr "$@"
\```

## What Changes?

### Commands That Work Unchanged
- ✅ `sadie airr`
- ✅ `sadie renumbering`
- ✅ `sadie reference make`
- ✅ All Python API calls

### New Commands
- 🆕 `sadie germlines populate`
- 🆕 `sadie germlines status`

### Data Location
- **Old**: Downloaded from G3 API on demand
- **New**: Stored locally in `~/.sadie/germlines/` (or package location)

## Rollback to G3 (Temporary)

If you need to use G3 temporarily:

\```bash
export SADIE_USE_GERMLINES_MODULE=false
sadie airr input.fasta
\```

**Warning**: G3 support will be removed after 2026-06-01.

## FAQ

**Q: Will my existing results change?**
A: No. For IMGT data, results are identical. Custom sources may add new annotations.

**Q: How much disk space do I need?**
A: Approximately 500MB for human, 2-3GB for all 33 species.

**Q: Can I use both G3 and germlines?**
A: Yes, toggle with `SADIE_USE_GERMLINES_MODULE` environment variable (until 2026-06-01).

**Q: Do I need to re-download monthly?**
A: No. Use `sadie germlines populate --force` when you want updates.

**Q: What if my species isn't available?**
A: Check available species with `sadie germlines status`. Request additions via GitHub issues.
```

### 4. Provider Guide (docs/germlines/provider-guide.md)

**Priority**: P1 (Important)
**Target Audience**: Advanced users
**Estimated Length**: 1200-1500 words

**Content**:
- Overview of each provider (IMGT, OGRDB, VDJbase)
- Differences between providers
- When to use which provider
- Priority-based merging explained
- Configuring provider priority

### 5. Custom Sequences Guide (docs/germlines/custom-sequences.md)

**Priority**: P1 (Important)
**Target Audience**: Researchers with novel sequences
**Estimated Length**: 800-1000 words

**Content**:
- What are custom sequences?
- When to add custom sequences
- File format requirements
- Adding custom sequences step-by-step
- Validation and troubleshooting

### 6. Troubleshooting Guide (docs/germlines/troubleshooting.md)

**Priority**: P1 (Important)
**Target Audience**: All users
**Estimated Length**: 1000-1200 words

**Content**:
- Common error messages and solutions
- Database not found errors
- Download failures
- Species not available
- Performance issues
- Validation errors

**Example Structure**:
```markdown
# Troubleshooting Guide

## Common Issues

### Error: "Species 'human' not found in germlines module"

**Symptom**:
\```
ValueError: Species 'human' not found in germlines module at /path/to/internal_data/human.
Build germlines databases with: update_databases('human').
\```

**Cause**: Germlines data not downloaded

**Solution**:
\```bash
sadie germlines populate -p imgt -s human
\```

---

### Error: "Download failed for species"

**Symptom**:
\```
[ERROR] Failed to download human from imgt: Connection timeout
\```

**Causes**:
1. Network connectivity issues
2. IMGT server unavailable
3. Firewall blocking connection

**Solutions**:
1. Check internet connection
2. Retry with `sadie germlines populate -p imgt -s human`
3. Use `--force` flag if partially downloaded
4. Check firewall settings

---

### Warning: "G3 API is deprecated"

**Symptom**:
\```
WARNING: G3 API is deprecated. Set SADIE_USE_GERMLINES_MODULE=true.
G3 will be removed after 2026-06-01.
\```

**Solution**: Migrate to germlines module (see [Migration Guide](migration-guide.md))

---

### Slow downloads

**Symptom**: Download taking > 10 minutes

**Solutions**:
1. Download single species first: `sadie germlines populate -s human`
2. Check network speed
3. Use `--dry-run` to estimate download size
4. Resume failed downloads (automatic from checkpoint)

---

## Debug Mode

Enable detailed logging:

\```bash
sadie -vvvvv germlines populate
\```

This shows:
- Download progress per file
- BLAST database build steps
- Auxiliary file generation
- Validation checks
```

### 7. Architecture Documentation (docs/germlines/architecture.md)

**Priority**: P2 (Nice to have)
**Target Audience**: Developers and contributors
**Estimated Length**: 1500-2000 words

**Content**:
- Technical architecture overview
- Staged pipeline explanation
- Provider interface design
- Database builder components
- Integration points with AIRR/renumbering
- Code structure

## Updates to Existing Documentation

### Update docs/index.md

Add new section:

```markdown
## Germline Database Management

SADIE now includes local germline database management, replacing the G3 API dependency.

- **[Germlines Overview](germlines/index.md)** - Get started with local germlines
- **[CLI Reference](germlines/cli-reference.md)** - Complete command reference
- **[Migration Guide](germlines/migration-guide.md)** - Migrate from G3 API

```

### Update docs/annotation.md

Add section about germlines backend:

```markdown
## Germline Data Sources

SADIE supports multiple germline data sources:

- **IMGT** (default): International ImMunoGeneTics information system
- **OGRDB**: Open Germline Receptor Database
- **VDJbase**: Curated germline database
- **Custom**: Your own sequences

See [Germline Provider Guide](germlines/provider-guide.md) for details.

### Using Local Germlines (Recommended)

\```python
from sadie.airr import Airr

# Uses local germlines automatically (default)
airr = Airr("human")
results = airr.run_file("sequences.fasta")
\```

### Using Legacy G3 API (Deprecated)

\```bash
export SADIE_USE_GERMLINES_MODULE=false
sadie airr sequences.fasta
\```

**Note**: G3 support will be removed after 2026-06-01.
```

## Implementation Plan

### Phase 1: Critical Documentation (Week 1)
- [ ] docs/germlines/index.md
- [ ] docs/germlines/cli-reference.md
- [ ] docs/germlines/migration-guide.md
- [ ] Update docs/index.md
- [ ] Update docs/annotation.md

### Phase 2: Detailed Guides (Week 2)
- [ ] docs/germlines/provider-guide.md
- [ ] docs/germlines/custom-sequences.md
- [ ] docs/germlines/troubleshooting.md

### Phase 3: Developer Documentation (Week 3)
- [ ] docs/germlines/architecture.md
- [ ] Code examples in docstrings
- [ ] API reference generation

### Phase 4: Review and Polish (Week 4)
- [ ] Technical review
- [ ] User testing with example workflows
- [ ] Screenshots and diagrams
- [ ] Link validation
- [ ] Proofreading

## Success Criteria

- [ ] All new CLI commands documented with examples
- [ ] Migration path from G3 is clear and tested
- [ ] Common error messages have solutions in troubleshooting
- [ ] Users can complete setup without external help
- [ ] Documentation passes technical review
- [ ] Links are valid and render correctly
- [ ] Code examples are tested and working

## Documentation Style Guide

### Voice and Tone
- Clear, concise, professional
- Action-oriented (use imperatives: "Download...", "Check...")
- Friendly but not overly casual

### Code Examples
- Always include complete, runnable examples
- Show both command-line and Python API usage
- Include expected output where helpful

### Formatting
- Use code blocks with syntax highlighting
- Include command prompts (`$` for bash)
- Use tables for comparing options
- Add warnings/notes with emoji indicators (⚠️, ℹ️, ✅)

### Structure
- Start each page with brief overview
- Include "Prerequisites" section when needed
- End with "Next Steps" or "See Also" links

## Review Process

1. **Self-review**: Author reviews for completeness
2. **Technical review**: Developer validates accuracy
3. **User testing**: Non-expert user follows instructions
4. **Final polish**: Fix issues, improve clarity

## Related Specifications

- [002-germline-integration/spec.md](../002-germline-integration/spec.md) - Feature specification
- [002-germline-integration/plan.md](../002-germline-integration/plan.md) - Implementation plan
- [src/sadie/germlines/INTEGRATION_GUIDE.md](../../src/sadie/germlines/INTEGRATION_GUIDE.md) - Developer integration guide

## Notes

- Documentation should be complete before the 2026-06-01 G3 deprecation deadline
- Include migration warnings in appropriate places
- Consider creating video tutorials for complex workflows (future work)
- Update README.md with germlines quick start section
