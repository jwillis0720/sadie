# Phase 28 Research: Fix Germline Priority Order

## Summary

This phase involves a simple constant change in `GermlineManager.DEFAULT_PROVIDERS` from `["custom", "ogrdb", "vdjbase", "imgt"]` to `["vdjbase", "ogrdb", "imgt", "custom"]`. The change affects one line of code but requires updating multiple documentation files and tests that assert the expected order.

**Primary Recommendation**: Update the constant, then systematically update all documentation and tests that reference the priority order. This is a straightforward refactoring task with no external library dependencies.

---

## Standard Stack

| Component | Version | Purpose |
|-----------|---------|---------|
| Python | 3.10+ | Runtime (existing) |
| pytest | 7.x | Test verification (existing) |
| No new dependencies required | | |

---

## Architecture Patterns

### Pattern 1: Priority-First-Wins

The `GermlineManager` implements a **first-provider-wins** pattern:

```python
DEFAULT_PROVIDERS = ["vdjbase", "ogrdb", "imgt", "custom"]

# Iteration order determines priority
for provider in self.providers:
    genes = self._fetch_from_provider(provider, species, segment, chain, functional_only)
    for gene in genes:
        if self._should_include_gene(gene, all_genes, seq_to_gene):
            all_genes[gene.name] = gene  # First wins on name conflict
            seq_to_gene[gene.sequence] = gene.name  # First wins on sequence match
```

**Key insight**: The list order directly maps to priority. First = highest priority.

### Pattern 2: Deduplication Rules

The manager uses two deduplication rules in `_should_include_gene()`:
1. **Name conflict**: If gene name already exists → skip (first provider wins)
2. **Sequence match**: If exact sequence already exists → skip (prevents duplicates)

This means changing priority order changes which allele variant is used when providers have the same gene name with different sequences.

### Pattern 3: Single-Provider-Per-Run

Importantly, the priority applies to ALL segments (V, D, J, C) within a single run. Per FR-014, there's no per-segment provider mixing:

```python
# All queries use same provider list
v_genes = manager.get_genes("human", "V", "H")
d_genes = manager.get_genes("human", "D", "H")  # Same priority as V
j_genes = manager.get_genes("human", "J", "H")  # Same priority as V
```

---

## Don't Hand-Roll

| Problem | Use Instead | Location |
|---------|-------------|----------|
| Priority logic | Existing `_should_include_gene()` | `manager.py:216-238` |
| Provider iteration | Existing `for provider in self.providers` | `manager.py:187` |
| Test assertions | Existing `test_compliance.py` patterns | `tests/unit/germlines/` |

**Do not** create a new priority resolution mechanism. The existing first-wins pattern is simple and effective.

---

## Common Pitfalls

### Pitfall 1: Incomplete Documentation Updates

The priority order is documented in **8+ locations**. Missing one creates inconsistency:

| File | Line(s) | Current Text |
|------|---------|--------------|
| `src/sadie/germlines/manager.py` | 28, 41, 50, 60 | docstrings reference old order |
| `src/sadie/germlines/__init__.py` | 12, 102 | docstrings say "custom > imgt > ogrdb > vdjbase" |
| `src/sadie/germlines/README.md` | 49 | "custom > ogrdb > vdjbase > imgt" |
| `src/sadie/germlines/sources/custom/README.md` | 113 | "custom > imgt > ogrdb" (outdated) |
| `docs/germlines/architecture.md` | 48, 220, 232 | references old priority |
| `docs/germlines/provider-guide.md` | 170 | priority order documentation |

### Pitfall 2: Test Assertion Mismatches

Two test files explicitly assert the expected priority order:

```python
# tests/unit/germlines/test_compliance.py:28
expected = ["custom", "ogrdb", "vdjbase", "imgt"]  # MUST CHANGE

# tests/unit/germlines/test_manager.py:29
assert manager.provider_names == ["custom", "ogrdb", "vdjbase", "imgt"]  # MUST CHANGE
```

### Pitfall 3: Custom Provider Documentation Contradiction

The custom provider README states:
> **Custom sequences have HIGHEST priority!**

With the new order `["vdjbase", "ogrdb", "imgt", "custom"]`, custom becomes LOWEST priority. This documentation must be updated to explain the new rationale (custom fills gaps, not overrides).

### Pitfall 4: Pipeline Source Change Detection

`pipeline.py:168` iterates sources in order `["custom", "ogrdb", "vdjbase", "imgt"]` for file change detection:

```python
for source_name in ["custom", "ogrdb", "vdjbase", "imgt"]:
    if self._source_newer_than(source_name, species, latest_normalized):
        return True
```

This is for change detection only (not priority), but should be updated for consistency.

### Pitfall 5: Breaking Existing User Workflows

If users have custom sequences designed to override provider sequences, the new priority order will break their workflows. Consider:
- Adding a migration note in docs
- The change is intentional - VDJbase has best validated human/macaque alleles

---

## Code Examples

### Example 1: The Change in manager.py

```python
# BEFORE (line 50)
DEFAULT_PROVIDERS = ["custom", "ogrdb", "vdjbase", "imgt"]

# AFTER
DEFAULT_PROVIDERS = ["vdjbase", "ogrdb", "imgt", "custom"]
```

### Example 2: Updated Test Assertion

```python
# tests/unit/germlines/test_compliance.py
def test_default_priority_order(self):
    """Verify default priority is vdjbase > ogrdb > imgt > custom."""
    expected = ["vdjbase", "ogrdb", "imgt", "custom"]
    assert (
        GermlineManager.DEFAULT_PROVIDERS == expected
    ), f"Expected priority {expected}, got {GermlineManager.DEFAULT_PROVIDERS}"
```

### Example 3: Updated Docstring in manager.py

```python
class GermlineManager:
    """
    Manages multiple germline databases with priority-based lookup.

    Default priority: vdjbase > ogrdb > imgt > custom
    - VDJbase provides curated, validated alleles (best for human/macaque)
    - OGRDB adds community-curated novel alleles (excellent for mouse)
    - IMGT provides comprehensive reference coverage
    - Custom fills gaps from internal lab data
    """
```

### Example 4: Rationale Comment in Code

```python
# Default provider priority order (highest to lowest):
# 1. VDJbase: Best for human/macaque - curated, validated alleles from population studies
# 2. OGRDB: Good for mouse - community-curated novel alleles
# 3. IMGT: Species diversity - comprehensive reference database
# 4. Custom: Fill gaps - internal lab sequences for edge cases
DEFAULT_PROVIDERS = ["vdjbase", "ogrdb", "imgt", "custom"]
```

---

## Files Requiring Updates

### Primary Code Changes

| File | Change Type | Lines |
|------|-------------|-------|
| `src/sadie/germlines/manager.py` | Constant + docstrings | 28, 41, 50, 60 |

### Test Updates

| File | Change Type | Lines |
|------|-------------|-------|
| `tests/unit/germlines/test_compliance.py` | Assertions | 28, 36 |
| `tests/unit/germlines/test_manager.py` | Assertions | 28-29 |

### Documentation Updates

| File | Change Type |
|------|-------------|
| `src/sadie/germlines/__init__.py` | Docstrings |
| `src/sadie/germlines/README.md` | Priority section |
| `src/sadie/germlines/sources/custom/README.md` | Priority explanation |
| `src/sadie/germlines/pipeline.py` | Source iteration order (for consistency) |
| `docs/germlines/architecture.md` | Priority diagrams |
| `docs/germlines/provider-guide.md` | Priority documentation |

---

## Verification Checklist

- [ ] `GermlineManager.DEFAULT_PROVIDERS == ["vdjbase", "ogrdb", "imgt", "custom"]`
- [ ] `test_default_priority_order` passes with new expected value
- [ ] `test_manager_initializes_with_default_priority` passes
- [ ] Human queries return VDJbase alleles when available
- [ ] Macaque queries return VDJbase alleles when available
- [ ] Mouse queries return OGRDB alleles (VDJbase has no mouse data)
- [ ] Custom sequences still work (just at lowest priority)
- [ ] All documentation references updated consistently
- [ ] `pytest tests/unit/germlines/` passes

---

## Sources

| Source | Confidence | Notes |
|--------|------------|-------|
| `src/sadie/germlines/manager.py` | HIGH | Primary implementation |
| `tests/unit/germlines/test_compliance.py` | HIGH | Current expected values |
| `src/sadie/germlines/sources/vdjbase/README.md` | HIGH | VDJbase species support |
| `.planning/ROADMAP.md` Phase 28 | HIGH | Requirements specification |
| `src/sadie/germlines/README.md` | MEDIUM | May be outdated |
| `docs/germlines/*.md` | MEDIUM | May be outdated |

---

## Risk Assessment

**Risk Level**: LOW

- Simple constant change with clear scope
- No new dependencies
- No architectural changes
- Existing test patterns can be reused
- All locations identified through grep search

**Mitigation**: Run full germlines test suite after each file change to catch assertion mismatches early.
