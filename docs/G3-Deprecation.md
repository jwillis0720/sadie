# G3 API Deprecation Notice

## Overview

The G3 API (Germline Gene/Genome database API) is deprecated and scheduled for removal from SADIE. This document provides migration guidance for users and developers.

## Timeline

| Milestone | Date | Status |
|-----------|------|--------|
| G3 API deprecated | 2025-01-01 | ✅ Complete |
| Germlines module default | 2025-01-22 | ✅ Complete |
| G3 API removal | After 2026-06-01 | ⏳ Scheduled |

## Migration Path

### For Users

The germlines module is now the default data provider. No action is required unless you explicitly disabled it.

**If you were using `SADIE_USE_GERMLINES_MODULE=false`:**

1. Remove the environment variable to use the new default
2. Or update to `SADIE_USE_GERMLINES_MODULE=true` explicitly

```bash
# Remove legacy override (recommended)
unset SADIE_USE_GERMLINES_MODULE

# Or explicitly enable germlines module
export SADIE_USE_GERMLINES_MODULE=true
```

### For Developers

Replace direct G3 API calls with the germlines module:

```python
# Before (deprecated)
from sadie.airr.igblast.germline import GermlineData
gd = GermlineData("human")
# Accessing gd.v_gene_dir, gd.aux_path, etc.

# After (current)
from sadie.germlines.manager import GermlineManager
manager = GermlineManager()
genes = manager.get_genes(species="human", locus="IGH", segment="V")
```

## Environment Variable Behavior

### `SADIE_USE_GERMLINES_MODULE`

Controls which backend is used for germline data:

| Value | Behavior |
|-------|----------|
| Not set | Uses germlines module (default) |
| `true` | Uses germlines module |
| `false` | Uses deprecated G3 API (shows deprecation warning) |

### Deprecation Warning

When `SADIE_USE_GERMLINES_MODULE=false`, SADIE logs a deprecation warning:

```
DeprecationWarning: G3 API is deprecated and will be removed after 2026-06-01.
Set SADIE_USE_GERMLINES_MODULE=true or unset to use the germlines module.
```

## Benefits of Germlines Module

The germlines module provides several advantages over the legacy G3 API:

1. **Multi-source support**: Access alleles from IMGT, OGRDB, VDJbase, and custom sources
2. **Configurable priority**: Control which source takes precedence for overlapping alleles
3. **Up-to-date data**: Uses current IMGT GENE-DB (vs G3's frozen snapshot)
4. **Offline capability**: Cached data eliminates network dependency
5. **Reference.yml integration**: Build custom germline databases from YAML configuration

## What Will Be Removed

After 2026-06-01, the following will be removed:

- G3 API fallback code paths in `src/sadie/airr/igblast/germline.py`
- `SADIE_USE_GERMLINES_MODULE` feature flag (germlines module becomes only option)
- Legacy G3-related tests and fixtures

## Questions?

For migration assistance, please open an issue on the SADIE repository.

---

*Last updated: 2026-01-25*
