---
created: 2026-01-26T19:11
title: Add prefix matching for IMGT C region downloads
area: germlines
files:
  - src/sadie/germlines/scripts/download_imgt.py:420-460
  - src/sadie/germlines/sources/imgt/macaque/species_mapping.json
---

## Problem

IMGT C region downloads are missing alleles like `IGHA*24` for macaque because of exact species name matching.

**Root cause:**
- V-QUEST (V/D/J source) has 41 macaque strain variants
- GENE-DB (C region source) has 59 macaque strain variants
- Species mapping JSON only contains V-QUEST variants
- `IGHA*24` requires `Macaca mulatta_RUp15` which only exists in GENE-DB

**Missing strains in GENE-DB but not V-QUEST:**
- `Macaca mulatta_RUp15` (has IGHA*24)
- `Macaca mulatta_Chinese RGr/RJz7/RUa8`
- `Macaca mulatta_Indian RGr/RLi1/RWy2`
- `Macaca mulatta_PMac`, `_UTS378`, `_UTS380`, etc.

**Current behavior:**
```
WARNING:Reference:Gene IGHA*24 not found in imgt database for macaque
```

## Solution

Change `_parse_genedb_header()` in `download_imgt.py` to use **prefix matching** instead of exact matching:

1. Extract base species name from mapping (e.g., "Macaca mulatta" from "Macaca mulatta_17573")
2. In C region filtering, match any GENE-DB species that starts with the base name
3. This will capture all strains like `Macaca mulatta_RUp15` for macaque

Alternative: Collect species variants from GENE-DB during first download and merge with V-QUEST variants in the JSON mapping.
