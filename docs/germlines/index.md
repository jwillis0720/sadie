# Germline Database Management

## Overview

The germlines module provides local germline database management for SADIE, replacing the dependency on the G3 API. This enables offline operation, multi-source germline data, and complete control over your germline reference databases.

### What is the Germlines Module?

The germlines module is a self-contained system for managing immunoglobulin germline gene databases. It downloads, processes, and maintains germline sequences from multiple sources:

- **IMGT** (International ImMunoGeneTics information system)
- **OGRDB** (Open Germline Receptor Database)
- **VDJbase** (Curated germline database)
- **Custom** sequences (your own novel germlines)

All data is stored locally, allowing SADIE to operate completely offline after the initial setup.

### Why Use Local Germlines?

The germlines module offers several advantages over the legacy G3 API:

**✅ Offline Operation**
Work without internet connectivity. Once populated, all germline data is available locally with no API dependency.

**✅ Multi-Source Support**
Combine sequences from IMGT, OGRDB, VDJbase, and your own custom germlines. Priority-based merging ensures your custom sequences take precedence.

**✅ Faster Performance**
Eliminate network latency. Local database access is instant with no API round-trips.

**✅ Reproducibility**
Lock germline versions for reproducible analysis. Your results won't change due to upstream database updates.

**✅ Extended Species Coverage**
Support for 33+ species including human, mouse, rat, rabbit, rhesus macaque, and more.

!!! warning "G3 API Deprecation"
    The G3 API is deprecated and will be removed after **2026-06-01**. Please migrate to the local germlines module. See the [Migration Guide](migration-guide.md) for step-by-step instructions.

## Quick Start

### Initial Setup

Download and build germline databases for all species and providers:

```bash
sadie germlines populate
```

This one-time setup takes 5-10 minutes and downloads:
- IMGT data for 33 species
- OGRDB data for human and mouse
- VDJbase data (if available)

The command builds BLAST databases, auxiliary files, and internal data structures automatically.

### Verify Installation

Check the status of your germline databases:

```bash
sadie germlines status
```

You should see output similar to:

```
┏━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━┓
┃ Provider ┃ Version       ┃ Downloaded   ┃ Species ┃ Status         ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━┩
│ imgt     │ release-202601│ 2026-01-22   │ 33      │ Up-to-date     │
│ ogrdb    │ release-202601│ 2026-01-22   │ 2       │ Up-to-date     │
│ vdjbase  │ -             │ -            │ 0       │ Not downloaded │
└──────────┴───────────────┴──────────────┴─────────┴────────────────┘
```

### Use with AIRR Annotation

The germlines module is now the default. No code changes needed:

```python
from sadie.airr import Airr

# Automatically uses local germlines
airr = Airr("human")
results = airr.run_file("sequences.fasta")
```

Command-line usage also works automatically:

```bash
sadie airr sequences.fasta -o results.tsv
```

## Selective Species Download

Download specific species only:

```bash
# Download human only from IMGT
sadie germlines populate -p imgt -s human

# Download multiple species
sadie germlines populate -p imgt -s human -s mouse
```

## Next Steps

### For New Users
1. ✅ Run `sadie germlines populate` to set up databases
2. ✅ Verify with `sadie germlines status`
3. ✅ Use SADIE normally - germlines work automatically

### For Existing G3 Users
- **[Migration Guide](migration-guide.md)** - Step-by-step migration from G3 API
- **[CLI Reference](cli-reference.md)** - Complete command documentation

### Advanced Topics
- **[Provider Guide](provider-guide.md)** - Understanding IMGT vs OGRDB vs VDJbase
- **[Custom Sequences](custom-sequences.md)** - Adding your own novel germlines
- **[Troubleshooting](troubleshooting.md)** - Common issues and solutions
- **[Architecture](architecture.md)** - Technical details for developers

## Support

- **Issues**: [GitHub Issues](https://github.com/jwillis0720/sadie/issues)
- **Documentation Feedback**: Use our [feedback template](https://github.com/jwillis0720/sadie/issues/new?template=documentation-feedback.md)

## Key Features

| Feature | G3 API (Deprecated) | Local Germlines |
|---------|---------------------|-----------------|
| Requires Internet | Always | Only for initial setup |
| Multi-source Support | No | Yes (IMGT, OGRDB, VDJbase, custom) |
| Custom Germlines | No | Yes |
| Species Coverage | Limited | 33+ species |
| Performance | Network dependent | Instant (local) |
| Reproducibility | Version can change | Fixed versions |
| Maintenance | Requires G3 server | Self-contained |

**Ready to get started?** Run `sadie germlines populate` and you're all set!
