# Research Summary: v1.2 Reference Module Unification

## Key Findings

### Stack
The Reference module already has foundational integration with the Germlines module via `use_germlines` flag—the core adapter (`GermlineToG3Adapter`) and manager queries are complete and functional. The primary blockers are **two small changes**: expanding source validation in `models.py` from `["imgt", "custom"]` to include `ogrdb` and `vdjbase`, and adding a `use_germlines` parameter to `References.from_yaml()`.

### Features

**Table Stakes (Must Have):**
| ID | Feature | Complexity |
|----|---------|------------|
| TS-1 | Expand source validation to 4 sources | Low |
| TS-2 | Route source selection through Germlines | Medium |
| TS-3 | Add `use_g3=False` parameter (deprecate G3) | Low |
| TS-4 | Create reference-sample.yml demo | Low |
| TS-5 | Source-specific gene validation | Medium |
| TS-6 | Clear error messages for missing data | Low |

**Differentiators (Nice to Have):**
| ID | Feature | Complexity |
|----|---------|------------|
| D-1 | Auto-discovery of available sources per species | Low |
| D-2 | Dry-run validation report for reference.yml | Medium |
| D-3 | Include source column in output DataFrame | Low |
| D-4 | Reference template generator CLI | Medium |

**Out of Scope:** Multi-provider blending per analysis, automatic G3 fallback, real-time provider sync, per-gene source overrides, TR support.

### Architecture

| Integration Point | Status | Action Needed |
|-------------------|--------|---------------|
| `Reference.__init__(use_germlines)` | ✅ Complete | None |
| `Reference._get_gene()` germlines path | ✅ Complete | None |
| `GermlineToG3Adapter.to_g3_format()` | ⚠️ Partial | Add `_id`, `chimera` fields |
| `References.from_yaml()` | ❌ Missing | Add `use_germlines` param |
| `models.py` source validation | ❌ Blocking | Expand `VALID_SOURCES` |

### Watch Out For

**P1: Source Validation Hardcoded** (🔴 CRITICAL)
- `models.py` rejects `ogrdb`/`vdjbase` at parse time
- Prevention: Expand `VALID_SOURCES = ["imgt", "ogrdb", "vdjbase", "custom"]`

**P2: Missing `_id` Field in Adapter** (🔴 CRITICAL)
- Downstream code uses `_id` for deduplication/indexing
- Prevention: Generate synthetic `_id` via hash of `source:species:gene`

**P5: Provider Priority vs Explicit Source** (🟡 HIGH)
- `get_gene_by_name()` ignores the `source` field, uses priority instead
- Prevention: Pass `providers=[gene.source]` to enforce explicit source

**P9: G3 API Blocking Constructor** (🟡 HIGH)
- Network call in `__init__` fails if G3 down, even when user wants germlines
- Prevention: Defer endpoint validation; allow germlines without G3 connectivity

**P6: Batch Query Inefficiency** (🟠 MEDIUM)
- Individual gene lookups O(n) per gene; slow for ~300 V genes
- Prevention: Add `get_genes_by_names()` batch method; target <2s for human

## Recommendations

**Implementation Order (Prioritized):**

1. **models.py:VALID_SOURCES expansion** — 5 min, unblocks everything
2. **reference.py:from_yaml() parameter** — 10 min, enables germlines path
3. **g3_adapter.py:add `_id` generation** — 10 min, fixes downstream breaks
4. **reference.py:pass explicit source to lookups** — 15 min, ensures reproducibility
5. **Test dual-backend parity** — 30 min, validates equivalence
6. **Create reference-sample.yml** — 15 min, demonstrates multi-source

**Key Principles:**
- No new external dependencies (all libraries already in pyproject.toml)
- Additive changes only—existing reference.yml files must continue working
- G3 API path remains functional; `use_germlines=False` for backwards compatibility
- Default to `use_germlines=True` for new behavior

## Files

- **STACK.md** — Stack analysis: current integration, required changes, implementation approach
- **FEATURES.md** — Feature categorization: table stakes vs differentiators vs anti-features
- **ARCHITECTURE.md** — Integration design: data flow diagrams, build order, risk assessment
- **PITFALLS.md** — Risk mitigation: 11 pitfalls with prevention strategies and phase mapping

---
*Synthesized: 2026-01-23*
