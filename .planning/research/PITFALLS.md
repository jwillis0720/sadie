# Pitfalls Research: Reference Module Unification

> Research for v1.2 milestone: Integrating Reference module with Germlines module

## Executive Summary

The v1.2 Reference Module Unification presents **11 critical pitfalls** spanning source validation, data format consistency, error handling, and backwards compatibility. The highest-risk issues involve source validation (currently hardcoded to `imgt|custom`), gene name format differences between G3 and Germlines, and missing G3 fields in the adapter. This document provides actionable prevention strategies with phase mappings.

---

## Critical Pitfalls

### P1: Source Validation Hardcoded to `imgt|custom`

**Risk Level:** 🔴 CRITICAL

**Current Code (src/sadie/reference/models.py:36-38):**
```python
@field_validator("source")
def check_source(cls, v: str) -> str:
    if v not in ["imgt", "custom"]:
        raise ValueError(f"{v} is not a valid source, choices are 'imgt' or 'custom'")
```

**Problem:** The `GeneEntry` and `GeneEntries` models reject `ogrdb` and `vdjbase` sources at validation time. Any reference.yml using these sources will fail at YAML parse, not runtime.

**Warning Signs:**
- ValidationError on `References.from_yaml()` call
- Error message: `"{source}" is not a valid source, choices are 'imgt' or 'custom'`
- Fails silently in development if only testing with imgt sources

**Prevention Strategy:**
1. Expand validator to accept: `["imgt", "ogrdb", "vdjbase", "custom"]`
2. Create constant `VALID_SOURCES` in models.py
3. Add unit tests for all four sources
4. Validate source availability at runtime (gene fetch) vs parse time (structure only)

**Phase:** Infrastructure (Phase 1)

---

### P2: Missing `_id` Field in G3 Adapter Output

**Risk Level:** 🔴 CRITICAL

**Current Code (src/sadie/germlines/g3_adapter.py):** The adapter builds a dict with fields like `source`, `common`, `gene`, etc., but does **not** include `_id`.

**Problem:** G3 API responses include `_id` field (MongoDB ObjectId). Reference module code uses this for:
- Deduplication: `_df = _df[~_df["_id"].duplicated()]` (reference.py:364)
- DataFrame indexing in `from_dataframe()` (reference.py:259)

**Warning Signs:**
- `KeyError: '_id'` during `get_dataframe()` operations
- Deduplication logic fails silently
- DataFrame validation fails in `from_dataframe()`

**Prevention Strategy:**
1. Generate synthetic `_id` in adapter: `hashlib.md5(f"{source}:{species}:{gene}".encode()).hexdigest()`
2. Ensure uniqueness via source+species+gene composite
3. Add `_id` to adapter unit tests
4. Document that `_id` format differs between backends

**Phase:** Adapter Enhancement (Phase 2)

---

### P3: Gene Name Format Inconsistencies

**Risk Level:** 🟡 HIGH

**Problem:** Gene names may have subtle format differences between G3 and Germlines:

| Issue | G3 Format | Germlines Format |
|-------|-----------|------------------|
| Allele separator | `*01` | `*01` ✓ |
| Family separator | `-` | `-` ✓ |
| Pseudogene suffix | `IGHV1-69P` | May be `IGHV1-69*01 (P)` |
| Chimera prefix | `human\|IGHV1-69*01` | Not implemented |

**Warning Signs:**
- Gene lookups fail despite gene existing
- Parity tests show name mismatches
- YAML validation passes but runtime lookup fails

**Prevention Strategy:**
1. Normalize gene names on input (strip whitespace, uppercase)
2. Create gene name comparison utility handling variants
3. Test with edge cases: pseudogenes, novel alleles, OGRDB names
4. Log warning when names differ slightly but match on sequence

**Phase:** Testing & Validation (Phase 3)

---

### P4: G3Error Raised from Germlines Path

**Risk Level:** 🟡 HIGH

**Current Code (src/sadie/reference/reference.py:194):**
```python
if not germline_gene:
    raise G3Error(f"Gene {gene.gene} not found in germlines database for species {gene.species}")
```

**Problem:** Using `G3Error` for germlines failures is semantically incorrect and confusing for debugging. Error messages reference "G3 database" when G3 wasn't involved.

**Warning Signs:**
- Stack traces mention G3 when user is using germlines backend
- Error handling logic that catches `G3Error` may not expect germlines failures
- Confusion in logs and debugging

**Prevention Strategy:**
1. Create `GermlineError` exception class for germlines-specific failures
2. Create common base `ReferenceError` with subclasses `G3Error`, `GermlineError`
3. Preserve backwards compatibility: catch both in existing handlers
4. Include backend type in error messages

**Phase:** Error Handling (Phase 2)

---

## Integration Pitfalls

### P5: Provider Priority vs Explicit Source Selection

**Risk Level:** 🟡 HIGH

**Germlines Default Priority:** `custom > ogrdb > vdjbase > imgt`

**Reference.yml Expects:** Explicit source per gene block

**Problem:** When `use_germlines=True`, the adapter bypasses source specification and uses provider priority. This violates the reference.yml contract where source is explicitly chosen.

**Current Code (reference.py:191-192):**
```python
germline_gene = get_gene_by_name(gene.gene, gene.species)  # Ignores gene.source!
```

**Warning Signs:**
- Gene from different source than specified in YAML
- User specifies `ogrdb` but gets `imgt` gene (if OGRDB data missing)
- Silent source mismatch breaks reproducibility

**Prevention Strategy:**
1. Pass `providers=[gene.source]` to `get_gene_by_name()` call
2. Enforce single-provider-per-lookup when source is explicit
3. Fail loudly if requested source doesn't have the gene
4. Add `providers` parameter to `get_gene_by_name()` function

**Phase:** Core Integration (Phase 2)

---

### P6: Batch Query Inefficiency

**Risk Level:** 🟠 MEDIUM

**Current Code (reference.py:247-255):**
```python
for gene_name in genes.genes:
    germline_gene = get_gene_by_name(gene_name, genes.species)
    # ... process each gene individually
```

**Problem:** Individual gene lookups are O(n) per gene. For V genes (~300), this means ~300 separate file parses.

**Contrast with G3:** Single API call with `limit=-1` returns all genes.

**Warning Signs:**
- Slow reference.yml processing (>10s for human)
- High I/O during bulk operations
- Memory spikes from repeated FASTA parsing

**Prevention Strategy:**
1. Add `get_genes_by_names(names: List[str], species: str)` to manager
2. Implement batch lookup: load FASTA once, filter by name set
3. Consider caching parsed FASTA in memory
4. Benchmark: target <2s for full human reference build

**Phase:** Performance Optimization (Phase 4)

---

### P7: Region Data Incompleteness

**Risk Level:** 🟠 MEDIUM

**Problem:** G3 API always returns region data (CDR1, FWR1, etc.) for V genes. Germlines module may have incomplete region data if:
- Source FASTA lacks IMGT-gapped sequences
- Provider doesn't supply region annotations
- Parse errors during region extraction

**G3 Adapter current behavior (g3_adapter.py:119-139):**
```python
if gene.regions and gene.region_positions:
    self._add_regions_to_imgt(imgt_dict, gene)
# No fallback if regions missing!
```

**Warning Signs:**
- `imgt.fwr1`, `imgt.cdr1` columns missing or empty
- IgBLAST internal data file generation fails
- Aux file generation errors (needs J gene regions)

**Prevention Strategy:**
1. Validate region data completeness before adapter use
2. Add fallback: derive regions from gapped sequence if positions missing
3. Fail fast with clear error if V gene lacks required regions
4. Add region completeness check to provider quality gate

**Phase:** Data Validation (Phase 3)

---

### P8: YAML Structure vs Germlines API Mismatch

**Risk Level:** 🟠 MEDIUM

**YAML Structure:**
```yaml
name:           # e.g., "human"
  source:       # e.g., "imgt"
    species:    # e.g., "human"
      - genes   # list of gene names
```

**Germlines API:**
```python
get_genes(species, segment, chain)  # No "name" concept
get_gene_by_name(name, species)     # No source parameter
```

**Problem:** The "name" level in YAML (used for chimeric references) has no equivalent in germlines. The adapter must maintain name→species mapping.

**Warning Signs:**
- Chimeric references fail
- Gene prefixing (`species|gene`) not applied
- Multiple species under same name not supported

**Prevention Strategy:**
1. Preserve name handling in Reference level, not adapter
2. Pass species explicitly to adapter, apply name later
3. Test chimeric case: two species under one name
4. Document that chimeric requires separate calls per species

**Phase:** Architecture (Phase 1)

---

## Compatibility Pitfalls

### P9: G3 API Availability Assumption

**Risk Level:** 🟡 HIGH

**Current Code:** When `use_germlines=False` (default), the Reference class checks G3 API availability in `__init__`:

```python
@endpoint.setter
def endpoint(self, endpoint: str) -> None:
    while True:
        _get = requests.get(endpoint)  # Network call in constructor!
```

**Problem:** Constructor makes blocking network call. If G3 is down, initialization fails even if user wants germlines.

**Warning Signs:**
- Slow startup (5+ seconds) when G3 unavailable
- `G3Error` on import if network issues
- Cannot use germlines as fallback

**Prevention Strategy:**
1. Defer endpoint validation to first G3 call
2. Add `validate_endpoint=False` parameter
3. Allow `use_germlines=True` without G3 connectivity
4. Consider lazy initialization pattern

**Phase:** Compatibility (Phase 2)

---

### P10: Existing Tests Assume G3 Output Format

**Risk Level:** 🟠 MEDIUM

**Problem:** Integration tests in `tests/unit/reference/` may assert on G3-specific fields or exact values. Running with germlines backend could expose format differences.

**Tests at risk:**
- `test_reference.py` - Reference class tests
- `test_advanced_reference.py` - Complex reference scenarios
- `test_reference_integration.py` - Cross-module tests

**Warning Signs:**
- Test failures when `use_germlines=True`
- Assertion errors on G3-specific fields
- DataFrame column mismatch errors

**Prevention Strategy:**
1. Run existing tests with both backends (parameterized)
2. Abstract backend-specific assertions
3. Create backend-agnostic comparison helpers
4. Document expected differences between backends

**Phase:** Testing (Phase 3)

---

### P11: reference.yml Backwards Compatibility

**Risk Level:** 🟠 MEDIUM

**Problem:** Existing reference.yml files only use `imgt` and `custom` sources. Adding new sources must not break existing files.

**Current reference.yml structure validated:**
- Source must be `imgt` or `custom`
- Gene names must start with IG (position 3 must be V/D/J/C)
- Species must be valid

**Warning Signs:**
- Existing reference.yml files fail validation
- Source validation errors on previously-valid files
- Breaking change for downstream users

**Prevention Strategy:**
1. New sources are additive (never remove `imgt`, `custom`)
2. Default behavior unchanged: G3 backend, imgt source
3. Add feature flag for new sources if needed
4. Migration guide for users adopting new sources

**Phase:** Compatibility (Phase 1)

---

## Prevention Strategies Summary

### Phase 1: Infrastructure
| Pitfall | Action | Effort |
|---------|--------|--------|
| P1 | Expand source validation | S |
| P8 | Document name/species mapping | S |
| P11 | Ensure backwards compatibility | S |

### Phase 2: Core Integration
| Pitfall | Action | Effort |
|---------|--------|--------|
| P2 | Add `_id` generation to adapter | S |
| P4 | Create GermlineError exception | S |
| P5 | Pass explicit source to lookups | M |
| P9 | Defer G3 endpoint validation | M |

### Phase 3: Testing & Validation
| Pitfall | Action | Effort |
|---------|--------|--------|
| P3 | Gene name normalization | M |
| P7 | Region data validation | M |
| P10 | Dual-backend test suite | M |

### Phase 4: Performance
| Pitfall | Action | Effort |
|---------|--------|--------|
| P6 | Batch query optimization | L |

---

## Quality Gate Checklist

Before marking v1.2 complete, verify:

- [ ] P1: All four sources (imgt, ogrdb, vdjbase, custom) accepted in YAML
- [ ] P2: `_id` field present in all adapter output
- [ ] P3: Gene name comparison handles edge cases
- [ ] P4: Correct exception types for each backend
- [ ] P5: Explicit source honored (no priority fallback)
- [ ] P6: Batch operations <2s for human reference
- [ ] P7: Region data validated before IgBLAST build
- [ ] P8: Chimeric references work with germlines
- [ ] P9: Germlines usable without G3 connectivity
- [ ] P10: Existing tests pass with both backends
- [ ] P11: Existing reference.yml files still valid

---

*Research completed: 2026-01-23*
*Pitfalls: 11 identified (3 critical, 4 high, 4 medium)*
