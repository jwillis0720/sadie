# Implementation Plan: Germline Database Integration

**Branch**: `002-germline-integration` | **Date**: 2026-01-21 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/002-germline-integration/spec.md`

## Summary

Connect SADIE's new germline database module to existing `sadie.airr.Airr` and `sadie.renumbering.Renumbering` modules, enabling selection of germline providers (IMGT, OGRDB, VDJbase, custom) via a `germline_backend` parameter. The integration uses a feature flag (`SADIE_USE_GERMLINES_MODULE`) for backwards compatibility, with G3 API remaining the default while local germlines is opt-in. Test suite in `tests/unit/germlines/` mirrors critical path tests from airr and renumbering modules.

## Technical Context

**Language/Version**: Python 3.10+
**Primary Dependencies**: biopython (≥1.80), pyhmmer (^0.11.1), pydantic (≥2.0.0), pandas (≥1.5), PyYAML, click (≥8.0)
**Storage**: Local filesystem (FASTA files, BLAST databases, HMM models, YAML configs)
**Testing**: pytest with coverage
**Target Platform**: Linux/macOS (cross-platform Python)
**Project Type**: Single Python package with CLI
**Performance Goals**: Equivalent to G3 backend (NFR-001)
**Constraints**: No silent fallback to G3 on failure (NFR-002), offline-capable after initial setup
**Scale/Scope**: Supports human, mouse, rat, rabbit, macaque, dog species; V/D/J/C segments; H/K/L chains

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design and before implementation.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Provider-Based Architecture | ✅ PASS | Providers are independent; no cross-provider dependencies |
| II. Priority-Based Merging (NON-NEGOTIABLE) | ✅ PASS | Default priority (custom > ogrdb > vdjbase > imgt); duplicates drop lower priority with warning; no per-segment mixing (FR-014) |
| III. Local-First Operation | ✅ PASS | Runtime uses local data only; offline coverage in US4 |
| IV. Staged Pipeline Architecture | ✅ PASS | sources → normalized → igblast pipeline respected; no bypass |
| V. Integration Compatibility | ✅ PASS | Backward compatibility maintained via feature flag; G3 default preserved |

## Project Structure

### Documentation (this feature)

```text
specs/002-germline-integration/
├── spec.md              # Feature specification (complete)
├── plan.md              # This file
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output
├── quickstart.md        # Phase 1 output
├── contracts/           # Phase 1 output (Python API contracts)
└── tasks.md             # Phase 2 output (/speckit.tasks command)
```

### Source Code (existing repository structure)

```text
src/sadie/
├── germlines/                    # NEW: Local germline database module
│   ├── __init__.py              # Public API (get_germline_genes, get_manager, etc.)
│   ├── models.py                # GermlineGene, ProviderMetadata (pydantic)
│   ├── manager.py               # GermlineManager (multi-source lookup)
│   ├── pipeline.py              # GermlinePipeline (3-stage processing)
│   ├── renumbering_integration.py  # LocalHMMBuilder
│   ├── g3_adapter.py            # Format conversion to G3 response format
│   ├── providers/               # Data source handlers
│   │   ├── base.py, imgt.py, ogrdb.py, vdjbase.py, custom.py
│   ├── builders/                # Data processing
│   │   ├── gapper.py, blast.py, aux.py, hmm.py
│   ├── utils/
│   │   └── feature_flags.py     # SADIE_USE_GERMLINES_MODULE flag
│   ├── sources/                 # Raw germline data (FASTA)
│   ├── normalized/              # Processed sequences
│   ├── igblast/                 # BLAST-ready databases
│   └── hmms/                    # Cached HMM models
│
├── airr/                         # EXISTING: AIRR annotation
│   ├── igblast/
│   │   └── germline.py          # ✅ INTEGRATED: GermlineData paths updated
│   └── ...
│
├── renumbering/                  # EXISTING: Antibody numbering
│   ├── aligners/
│   │   └── hmmer.py             # ✅ INTEGRATED: LocalHMMBuilder support
│   └── ...
│
└── reference/                    # EXISTING: Reference system
    └── reference.py             # ✅ INTEGRATED: GermlineManager + G3Adapter

tests/
├── unit/
│   ├── airr/                    # Existing AIRR tests (G3-based)
│   ├── renumbering/             # Existing renumbering tests (G3-based)
│   └── germlines/               # NEW: Mirrored tests using germlines backend
│       ├── test_airr_integration.py       # ✅ CREATED
│       ├── test_renumbering_integration.py # ✅ CREATED
│       └── test_reference_integration.py   # ✅ CREATED
└── integration/
    └── ...
```

**Structure Decision**: Existing single-package structure maintained. New `germlines/` module added under `src/sadie/`. Test suite mirrored in `tests/unit/germlines/`.

## Complexity Tracking

No constitution violations requiring justification. The integration follows established patterns:
- Feature flag for gradual migration (standard practice)
- Adapter pattern for format compatibility (G3Adapter)
- Provider pattern for multiple data sources (already implemented in germlines module)

## Integration Status (Current)

Based on current tasks, the integration is **81% complete** (60/74 tasks):

| Component | File | Status |
|-----------|------|--------|
| IgBLAST paths | `airr/igblast/germline.py` | ✅ Complete |
| Renumbering HMM | `renumbering/aligners/hmmer.py` | ✅ Complete |
| Reference system | `reference/reference.py` | ✅ Complete |
| Feature flag | `germlines/utils/feature_flags.py` | ✅ Complete |
| Test suite | `tests/unit/germlines/` | ✅ 91% Complete |
| BLAST databases | `germlines/igblast/` | ⏳ Human only (Phase 10 pending) |
| Documentation | INTEGRATION_GUIDE.md | ✅ Complete |
| Multi-species | 33 species support | 🚧 Phase 10 in progress |

## Remaining Work

### Phase 10: Species Expansion (Current Focus)
1. **Download IMGT data** for all 33 species in SPECIES_MAP
2. **Build IgBLAST BLAST databases** for each downloaded species
3. **Create auxiliary files** (*.aux) for J gene CDR3 start positions
4. **Update organism.yaml** with all species configurations
5. **Multi-species testing** - AIRR and renumbering tests for mouse, rabbit, rhesus_macaque, chicken/zebrafish

### Outstanding Gaps
- T004a: Verify gapped AA/NT sequences for all V/J genes
- T035a: Test gapped AA fallback translation
