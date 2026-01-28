# Phase 28 Summary: Fix Germline Priority Order

## Completion Date
2026-01-25

## Goal
Reorder germline provider priority from `["custom", "ogrdb", "vdjbase", "imgt"]` to `["vdjbase", "ogrdb", "imgt", "custom"]` for optimal data quality.

## Rationale
- **VDJbase**: Best for human/macaque - curated, validated alleles from population studies
- **OGRDB**: Good for mouse - community-curated novel alleles
- **IMGT**: Species diversity - comprehensive reference database
- **Custom**: Fill gaps - internal lab sequences for edge cases

## Changes Made

### Wave 1 - Plan 1: Core Code + Tests
| Task | File | Change |
|------|------|--------|
| 28-1-1 | `src/sadie/germlines/manager.py` | Updated DEFAULT_PROVIDERS constant and all docstrings (3 locations) |
| 28-1-2 | `src/sadie/germlines/pipeline.py` | Updated hardcoded source iteration list |
| 28-1-3 | `tests/unit/germlines/test_compliance.py` | Updated test assertions and renamed test method |
| 28-1-4 | `tests/unit/germlines/test_manager.py` | Updated test assertion and comment |

### Wave 1 - Plan 2: Module Docstrings
| Task | File | Change |
|------|------|--------|
| 28-2-1 | `src/sadie/germlines/__init__.py` | Updated 3 docstring locations (lines 12, 20, 102) |

### Wave 2 - Plan 3: Documentation Updates
| Task | File | Change |
|------|------|--------|
| 28-3-1 | `src/sadie/germlines/README.md` | Updated priority section with rationale |
| 28-3-2 | `src/sadie/germlines/sources/custom/README.md` | Updated priority explanation ("fill gaps" not "HIGHEST") |
| 28-3-3 | `docs/germlines/architecture.md` | Updated diagram and reverse iteration code sample |
| 28-3-4 | `docs/germlines/provider-guide.md` | Updated priority order and example scenario |
| 28-3-5 | `src/sadie/germlines/sources/vdjbase/README.md` | Updated to show "HIGHEST priority" |
| 28-3-6 | `src/sadie/germlines/sources/ogrdb/OGRDB_DATA.md` | Updated priority section header and order |
| 28-3-7 | `.specify/memory/constitution.md` | Updated NON-NEGOTIABLE Principle II (ignored by git) |

## Verification Results

### Must-Haves Verified
| ID | Requirement | Result |
|----|-------------|--------|
| MH-1 | DEFAULT_PROVIDERS == ["vdjbase", "ogrdb", "imgt", "custom"] | ✅ PASS |
| MH-2 | Manager initializes with new default order | ✅ PASS |
| MH-3 | Priority tests pass with new expected values | ✅ PASS |
| MH-4 | No regressions in germlines tests (185 tests) | ✅ PASS |
| MH-5 | Docstrings in manager.py reflect new priority | ✅ PASS |
| MH-6 | Docstrings in __init__.py updated | ✅ PASS |
| MH-7 | pipeline.py has consistent source order | ✅ PASS |
| MH-8 | architecture.md code sample updated | ✅ PASS |
| MH-9 | VDJbase README reflects HIGHEST priority | ✅ PASS |
| MH-10 | OGRDB_DATA.md shows new priority order | ✅ PASS |
| MH-11 | Constitution updated | ✅ PASS (locally, gitignored) |

### Test Results
```
tests/unit/germlines/ - 185 passed, 134 warnings in 11.68s
```

## Commits
| Hash | Message |
|------|---------|
| 7c20450f | fix(28): change default germline priority to vdjbase > ogrdb > imgt > custom |

## Deviations
None. Plan executed as specified.

## Impact
- Human and macaque queries now prefer VDJbase alleles (curated, validated)
- Mouse queries continue to get OGRDB alleles (VDJbase has limited mouse data)
- Custom sequences now serve as gap-fillers rather than overrides
- To override a standard allele, users add custom sequence with same gene name

## Milestone v1.3 Progress
- Phase 25: ✅ Macaque Germlines Integration
- Phase 26: ✅ Add AIRR Package Dependency
- Phase 27: ✅ Remove Deprecated G3 Tests
- Phase 28: ✅ Fix Germline Priority Order

**Milestone v1.3 Status: 100% COMPLETE (4/4 phases)**
