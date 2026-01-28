# Features Research: Reference Module Unification

Research Date: 2026-01-23
Focus: Multi-source allele selection for reference.yml

## Executive Summary

The Reference module unification aims to enable `reference.yml` to select alleles from all 4 germline sources (imgt, ogrdb, vdjbase, custom) using the germlines module as the data provider instead of the G3 API.

**Current State:**
- reference.yml only supports `imgt` and `custom` sources
- G3 API is the default data provider
- `use_germlines=True` flag exists but only for IMGT data

**Target State:**
- reference.yml supports all 4 sources: imgt, ogrdb, vdjbase, custom
- Germlines module replaces G3 API as default (`use_g3=False`)
- Chimeric references can mix species from different sources

---

## Table Stakes

Must-have features for this milestone to work. These are non-negotiable requirements.

### TS-1: Expand Source Validation in YamlRef

**Description:** YamlRef currently only validates `imgt` and `custom` sources. Must accept `ogrdb` and `vdjbase`.

**Complexity:** Low
**Existing Dependencies:** `src/sadie/reference/yaml.py` (YamlRef class)
**Implementation Notes:**
- Update source validation logic to accept ["imgt", "ogrdb", "vdjbase", "custom"]
- Maintain backward compatibility with existing reference.yml files

---

### TS-2: Route Source Selection Through Germlines Module

**Description:** When `use_germlines=True` (soon default), Reference must route gene lookups through the appropriate germlines provider based on the `source` field in reference.yml.

**Complexity:** Medium
**Existing Dependencies:**
- `src/sadie/reference/reference.py` (Reference class, `_get_gene()`, `_get_genes()`)
- `src/sadie/germlines/manager.py` (GermlineManager)
- `src/sadie/germlines/g3_adapter.py` (GermlineToG3Adapter)

**Implementation Notes:**
- Currently `use_germlines` path ignores `source` field and uses manager priority
- Need to create source-specific lookup: `get_gene_by_name(name, species, source=source)`
- G3 adapter already handles transformation to G3 format

---

### TS-3: Add use_g3=False Parameter

**Description:** Add `use_g3` parameter to Reference class with `False` as default (soft deprecation of G3 API).

**Complexity:** Low
**Existing Dependencies:**
- `src/sadie/reference/reference.py` (Reference.__init__)

**Implementation Notes:**
- Rename `use_germlines` to internal, add public `use_g3` parameter
- `use_g3=False` → use germlines module
- `use_g3=True` → use G3 API (legacy path)
- Update docstrings to reflect deprecation timeline

---

### TS-4: Create reference-sample.yml

**Description:** Create sample reference.yml demonstrating multi-source selection with mouse=imgt, human=ogrdb, macaque=vdjbase.

**Complexity:** Low
**Existing Dependencies:**
- `src/sadie/reference/data/reference.yml` (existing format)
- Germlines module with populated data for all 3 species

**Implementation Notes:**
```yaml
sample:
  imgt:
    mouse:
      - IGHV1-12*01
      - IGHD1-1*01
      - IGHJ1*01
  ogrdb:
    human:
      - IGHV1-69*01
      - IGHD3-3*01
      - IGHJ6*01
  vdjbase:
    rhesus_macaque:
      - IGHV1-2*01
      - IGHD3-3*01
      - IGHJ4*01
```

---

### TS-5: Source-Specific Gene Validation

**Description:** Validate that requested genes exist in the specified source for the specified species.

**Complexity:** Medium
**Existing Dependencies:**
- Germlines providers (IMGT, OGRDB, VDJbase, Custom)
- `src/sadie/germlines/manager.py` (GermlineManager.get_gene_by_name)

**Implementation Notes:**
- Each provider has different species coverage:
  - IMGT: 29+ species
  - OGRDB: human, mouse (primarily)
  - VDJbase: human, rhesus_macaque
  - Custom: user-defined
- Clear error messages when source doesn't have requested species
- Clear error messages when gene not found in source

---

### TS-6: Error Handling for Missing Source Data

**Description:** Provide clear error messages when requested source lacks data for species or gene.

**Complexity:** Low
**Existing Dependencies:**
- ERR-01 requirement from PROJECT.md (already implemented for providers)

**Implementation Notes:**
```
G3Error: "Gene IGHV1-69*01 not found in source 'ogrdb' for species 'mouse'.
         OGRDB only supports: human. Available for mouse: imgt, custom."
```

---

## Differentiators

Nice-to-have features that improve user experience or future extensibility.

### D-1: Auto-Discovery of Available Sources per Species

**Description:** Utility to list which sources have data for a given species.

**Complexity:** Low
**Existing Dependencies:**
- `GermlineManager.get_available_species()` (per-manager, not per-provider)

**Implementation Notes:**
- Add method `get_available_sources(species: str) -> List[str]`
- Returns ["imgt", "ogrdb"] for human, ["imgt"] for dog, etc.
- Useful for CLI and error messages

---

### D-2: Validation Report for reference.yml

**Description:** Dry-run validation that checks all genes in reference.yml exist in their specified sources.

**Complexity:** Medium
**Existing Dependencies:**
- YamlRef class
- Germlines providers

**Implementation Notes:**
- CLI command: `sadie reference validate --yaml reference.yml`
- Reports: found/missing/deprecated genes per source
- Does not modify databases or perform actual lookups

---

### D-3: Source Metadata in Output DataFrame

**Description:** Include `source` column in Reference dataframe to track gene provenance.

**Complexity:** Low
**Existing Dependencies:**
- `Reference.get_dataframe()` method
- G3 adapter already preserves source field

**Implementation Notes:**
- Already partially implemented via G3 adapter
- Ensure source field flows through to final output

---

### D-4: Reference Template Generator

**Description:** Generate reference.yml templates based on available data.

**Complexity:** Medium
**Existing Dependencies:**
- Germlines manager
- Provider metadata

**Implementation Notes:**
- CLI command: `sadie reference generate-template --species human --source ogrdb`
- Outputs YAML with all available genes from that source
- Helps users discover what's available

---

### D-5: Multi-Source Priority Within Species

**Description:** Allow reference.yml to specify fallback order for genes within a species entry.

**Complexity:** High (Constitution change)
**Existing Dependencies:**
- Constitution principle: single provider per run

**Implementation Notes:**
- NOT recommended for this milestone
- Would violate "no per-segment mixing" principle
- Defer to future milestone if truly needed

---

## Anti-Features

Features deliberately NOT being built in this milestone to maintain scope and avoid complexity.

### AF-1: Multi-Provider Blending Per Analysis

**Status:** OUT OF SCOPE (Constitution principle)

**Rationale:** The germlines module constitution explicitly states "single provider selection per run (no per-segment mixing)". This milestone extends which sources reference.yml can specify, but does NOT change how the germlines module merges data.

**What we ARE doing:** reference.yml can specify different sources for different species entries (chimeric references).

**What we are NOT doing:** Automatically merging IMGT + OGRDB + VDJbase for a single species.

---

### AF-2: Automatic G3 API Fallback

**Status:** OUT OF SCOPE (ERR-04 requirement)

**Rationale:** "No silent G3 fallback" is a validated requirement. If germlines module fails, we raise an error rather than silently using G3.

---

### AF-3: Real-Time Provider Synchronization

**Status:** OUT OF SCOPE (existing constraint)

**Rationale:** Germlines data update is manual via `sadie germlines populate`. This milestone does not add automatic sync with upstream databases.

---

### AF-4: Per-Gene Source Override

**Status:** OUT OF SCOPE (complexity)

**Rationale:** A syntax like:
```yaml
human:
  ogrdb:
    human:
      - IGHV1-69*01
      - IGHD3-3*01@imgt  # override source for this gene
```
Would add significant complexity to YamlRef parsing. Defer to future if needed.

---

### AF-5: T-Cell Receptor (TR) Support

**Status:** OUT OF SCOPE (project constraint)

**Rationale:** Project scope is immunoglobulin (IG) only. TR germlines are a separate concern.

---

### AF-6: GUI for Reference Configuration

**Status:** OUT OF SCOPE (project constraint)

**Rationale:** CLI and API only per project requirements.

---

### AF-7: Database Auto-Building on Reference Change

**Status:** OUT OF SCOPE (user requirement)

**Rationale:** "DB building is separate step" is an explicit user requirement. Reference.yml changes do not trigger automatic IgBLAST database regeneration.

---

## Feature Dependencies

Visual representation of what depends on what:

```
┌─────────────────────────────────────────────────────────────────┐
│                    reference-sample.yml (TS-4)                  │
│                              │                                  │
│         Requires: TS-1, TS-2, TS-5 to function                 │
└─────────────────────────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│              Source Validation (TS-1)                           │
│                              │                                  │
│         Prerequisite for all other features                    │
└─────────────────────────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│         Route Through Germlines (TS-2)                          │
│                              │                                  │
│  Depends on: TS-1, Germlines Manager, G3 Adapter               │
│  Enables: TS-5, TS-6, D-1, D-2, D-3                            │
└─────────────────────────────────────────────────────────────────┘
                               │
              ┌────────────────┼────────────────┐
              ▼                ▼                ▼
┌───────────────────┐ ┌────────────────┐ ┌────────────────┐
│ use_g3 param      │ │ Gene Validation │ │ Error Handling │
│ (TS-3)            │ │ (TS-5)          │ │ (TS-6)         │
│                   │ │                 │ │                │
│ Independent       │ │ Depends: TS-2   │ │ Depends: TS-5  │
└───────────────────┘ └────────────────┘ └────────────────┘
```

### Implementation Order (Recommended)

1. **TS-1: Expand Source Validation** - Foundation for all source-related work
2. **TS-3: Add use_g3=False Parameter** - API surface change, can parallel TS-1
3. **TS-2: Route Through Germlines Module** - Core integration work
4. **TS-5: Source-Specific Validation** - Enables proper error handling
5. **TS-6: Error Messages** - Polish for user experience
6. **TS-4: reference-sample.yml** - Demonstration and documentation

### Differentiators Order (If Time Permits)

1. D-3: Source Metadata in DataFrame (low effort, high value)
2. D-1: Auto-Discovery of Sources (enables D-2)
3. D-2: Validation Report (CLI user experience)
4. D-4: Template Generator (CLI user experience)

---

## Complexity Summary

| ID | Feature | Complexity | Phase |
|----|---------|------------|-------|
| TS-1 | Source Validation | Low | Core |
| TS-2 | Germlines Routing | Medium | Core |
| TS-3 | use_g3 Parameter | Low | Core |
| TS-4 | Sample YAML | Low | Core |
| TS-5 | Gene Validation | Medium | Core |
| TS-6 | Error Handling | Low | Core |
| D-1 | Source Discovery | Low | Nice-to-have |
| D-2 | Validation Report | Medium | Nice-to-have |
| D-3 | Source Metadata | Low | Nice-to-have |
| D-4 | Template Generator | Medium | Nice-to-have |

**Estimated Core Effort:** 3-5 development units
**With Differentiators:** +2-3 development units

---

## Quality Gate Verification

- [x] Categories are clear (Table Stakes / Differentiators / Anti-Features)
- [x] Complexity noted for each feature (Low / Medium / High)
- [x] Dependencies on existing features identified (with paths)
- [x] Implementation order provided
- [x] Anti-features explicitly justified
