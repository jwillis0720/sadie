# Research Summary: Germline Database Integration

*Converted from spec-kit research.md*

## Key Findings

### Stack
- **Python 3.10+** with pyhmmer, Biopython, pydantic
- **pyhmmer** for HMM building (Stockholm format)
- **BLAST+** for IgBLAST database generation
- No new dependencies required

### Architecture Patterns

1. **Adapter Pattern**: GermlineGene → G3 format transformation
   - Maintains compatibility with existing Reference system
   - Isolated in `g3_adapter.py`

2. **Builder Pattern**: LocalHMMBuilder for renumbering
   - Queries germlines for V/J genes
   - Generates HMM via pyhmmer
   - Caches results

3. **Strategy Pattern**: Feature flag switching
   - `SADIE_USE_GERMLINES_MODULE` env var
   - Runtime selection of backend

### Integration Points

| Component | File | Integration Method |
|-----------|------|-------------------|
| IgBLAST | `airr/igblast/germline.py` | Path switching |
| HMM | `renumbering/aligners/hmmer.py` | LocalHMMBuilder |
| Reference | `reference/reference.py` | G3Adapter |

### Data Flow

**AIRR**: Request → Feature flag check → germlines/igblast/ or airr/data/ → IgBLAST

**Renumbering**: Request → Feature flag check → LocalHMMBuilder (cache/build) → HMMER

**Reference**: Request → use_germlines param → GermlineManager + G3Adapter or G3 API

## Watch Out For

1. **Gapped sequence availability** — V/J genes need gapped AA or gapped NT for HMM building
2. **BLAST+ dependency** — External tool required for database building
3. **Performance parity** — Must match G3 lookup times
4. **Provider priority** — custom > ogrdb > vdjbase > imgt (NON-NEGOTIABLE)

## Full Research

See `specs/002-germline-integration/research.md` for complete research documentation.

---
*Converted: 2026-01-21*
