# Technical Concerns and Debt

**Last Updated:** 2026-01-25
**Scope:** SADIE codebase technical debt, known issues, and integration risks

---

## Critical: Reference Module Gaps

### G3 Adapter Missing IMGT Region Fields

**Location:** `src/sadie/germlines/g3_adapter.py`
**Impact:** HIGH — `--use-germlines` flag in CLI is non-functional
**Status:** Documented gap, workaround exists

The `GermlineToG3Adapter` does not provide all IMGT region fields required by `make_airr_database()`:
- Missing: `imgt.fwr1_start`, `imgt.fwr1_end`, `imgt.cdr1_start`, `imgt.cdr1_end`, etc.
- These fields are required for creating IgBLAST internal annotation files (ndm.imgt)

**Affected Command:**
```bash
sadie reference build ref.yml -o /tmp/db --use-germlines  # FAILS
```

**Workaround:** Use without `--use-germlines` flag (uses G3 API path)

**Fix Approach:**
1. Enhance `_add_regions_to_imgt()` in `g3_adapter.py` lines 196-239
2. Ensure all region start/end positions are populated from GermlineGene.region_positions
3. Add integration test verifying all 52 IMGT fields present

**References:**
- `.planning/phases/phase-21/SUMMARY.md` line 51-53
- `.planning/ROADMAP.md` line 273

---

## Deprecation: G3 API Removal

**Deadline:** 2026-06-01
**Impact:** MEDIUM — Legacy code paths will be removed

The G3 API is deprecated and scheduled for removal:
- `src/sadie/airr/igblast/germline.py:96` — Deprecation warning
- `src/sadie/germlines/utils/feature_flags.py:46` — Deprecation warning

**Migration Status:**
- ✓ Default switched to germlines module (`SADIE_USE_GERMLINES_MODULE=true`)
- ✓ 98.29% structural parity achieved with G3
- ⚠ Some legacy code paths still reference G3

**Action Required:**
- Remove G3 fallback paths after deadline
- Update `src/sadie/airr/igblast/germline.py:203` legacy path
- Clean up G3-related tests and fixtures

---

## Known TODOs in Code

### High Priority

| Location | Description | Impact |
|----------|-------------|--------|
| `src/sadie/airr/methods.py:235` | `find_best_codon` edge case handling | Data quality |
| `src/sadie/airr/igblast/igblast.py:915` | Move to Pydantic model | Type safety |
| `src/sadie/renumbering/renumbering.py:105,121` | Move checks out of aligner class | Code organization |
| `src/sadie/germlines/providers/imgt.py:16,69,122,167,343` | Implement IMGT data download/parsing | Feature completeness |

### Medium Priority

| Location | Description | Impact |
|----------|-------------|--------|
| `src/sadie/airr/models/series.py:7` | Update to AIRR-community stable news | Standards compliance |
| `src/sadie/numbering/schemes.py:429` | N-terminal FWR2/3/4 edge cases | Numbering accuracy |
| `src/sadie/renumbering/clients/g3.py:16,28,29` | rhesus/dog species fixes | Species support |
| `src/sadie/renumbering/aligners/hmmer.py:35,324,340,355` | HMM model consolidation | Performance |
| `src/sadie/germlines/providers/ogrdb.py:327` | Extract version from OGRDB API | Version tracking |

### Low Priority

| Location | Description | Impact |
|----------|-------------|--------|
| `src/sadie/renumbering/clients/g3.py:187,275` | Refactor hand arch parsing | Code quality |
| `src/sadie/typing/species.py:7` | Verify viable species tests | Test coverage |

---

## Test Coverage Gaps

### Skipped Tests

| File | Reason | Risk |
|------|--------|------|
| `tests/unit/germlines/test_multi_species.py:197,238,278,359` | Species config incomplete | Medium |
| `tests/unit/germlines/test_compliance.py:262` | Aux file not found | Medium |

### Missing Test Coverage

- **Reference CLI `--use-germlines` path** — No integration test for germlines-based database build
- **Prebuilt database validation** — Limited tests for `Airr(database=<path>)` edge cases
- **Multi-source reference.yml** — No tests combining IMGT + OGRDB + VDJbase in single build

---

## Performance Concerns

### Unmeasured Performance

**Location:** `src/sadie/airr/airr.py` database parameter
**Issue:** Phase 22 success criteria states "Performance not measured"

The prebuilt database path should be faster than runtime lookup, but this hasn't been verified.

**Recommended Action:**
1. Add benchmark comparing `Airr(database=...)` vs default path
2. Document expected performance characteristics

### Potential Bottlenecks

| Component | Concern |
|-----------|---------|
| `Reference._get_genes()` | Pulls all genes then filters (see line in reference.py: `# @Todo, add a find_genes method to G3`) |
| `GermlineManager` multi-provider lookup | Sequential provider checking could be parallelized |
| IgBLAST execution | Process spawning overhead for small batches |

---

## Technical Debt Items

### Architecture Debt

1. **Dual Data Paths**
   - G3 API (legacy, deprecated)
   - Germlines module (current default)
   - Creates maintenance burden until G3 removal

2. **Adapter Pattern Complexity**
   - `GermlineToG3Adapter` transforms to legacy G3 format
   - Could be simplified after G3 removal to native format

3. **Species Mapping Duplication**
   - `src/sadie/reference/settings.py` has species mappings
   - `src/sadie/germlines/scripts/download_imgt.py:125` duplicates mappings
   - `src/sadie/typing/species.py` has another copy
   - Should consolidate to single source of truth

### Code Quality Debt

1. **IMGT Provider Incomplete**
   - `src/sadie/germlines/providers/imgt.py` has multiple TODO markers
   - Automatic download from IMGT not implemented
   - Header metadata parsing incomplete

2. **HMM Model Management**
   - `src/sadie/renumbering/aligners/hmmer.py:35` notes merge needed with G3-created HMMs
   - No clear documentation of which HMMs are legacy vs live-built

3. **Deprecated Exception Class**
   - `src/sadie/airr/exceptions.py:24` — BadIgBLASTInput deprecated
   - Dead code should be removed

---

## Security Considerations

### No Critical Issues Found

The codebase does not appear to handle sensitive data (passwords, API keys, tokens) directly.

### Recommendations

1. **Network Requests**
   - G3 API calls use plain HTTP requests
   - Consider certificate pinning for production deployments

2. **File Path Handling**
   - Database paths are user-provided
   - Path traversal validation exists in `validate_prebuilt_database()`

3. **External Process Execution**
   - IgBLAST executed via subprocess
   - Command construction appears safe (no shell injection vectors)

---

## Integration Risks

### IgBLAST Quirks

**Documented in:** `audit/igblast-quirk.md`

The `complete_vdj` field behavior differs based on:
- Internal data file structure (combined vs V-only)
- Database size/configuration

SADIE implements post-processing recalculation (`_recalculate_complete_vdj`) to normalize behavior.

**Risk:** Future IgBLAST versions may change behavior, requiring recalculation logic updates.

### IMGT Database Version Variance

**Documented in:** `audit/parity-notes.md`

The germlines module uses current IMGT GENE-DB (40 D alleles) vs G3's legacy snapshot (34 D alleles). This causes:
- 1.71% structural difference in annotations
- D gene allele selection differences
- Orphon gene detection differences

**Risk:** Accepted as improvement, not defect. Users comparing against legacy systems may see differences.

### Dependency on External Data Sources

| Source | Risk | Mitigation |
|--------|------|------------|
| IMGT GENE-DB | Manual download required | Bundled data in `sources/imgt/` |
| OGRDB API | Network dependency | Cached in `sources/ogrdb/` |
| VDJbase | Network dependency | Cached in `sources/vdjbase/` |
| G3 API (deprecated) | Network dependency | Will be removed after 2026-06-01 |

---

## Priority Recommendations

### Immediate (Before Next Release)

1. **Fix G3 adapter IMGT region fields** — Unblocks `--use-germlines` CLI flag
2. **Add integration test for prebuilt database path** — Validates Phase 22 implementation

### Short-term (Next Quarter)

1. **Consolidate species mappings** — Single source of truth
2. **Benchmark prebuilt database performance** — Document expected improvements
3. **Complete IMGT provider implementation** — Automatic download capability

### Long-term (After G3 Deprecation)

1. **Remove G3 API code paths** — Simplify architecture after 2026-06-01
2. **Refactor adapter to native format** — Eliminate G3 format transformation
3. **Clean up legacy tests and fixtures**

---

*Document generated: 2026-01-25*
*Sources: Code analysis, phase verification files, ROADMAP.md, audit reports*
