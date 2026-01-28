# Phase 28: Fix Germline Priority Order — Executable Plan

## Goal Statement

Reorder germline provider priority from `["custom", "ogrdb", "vdjbase", "imgt"]` to `["vdjbase", "ogrdb", "imgt", "custom"]` for optimal data quality.

**Rationale:**
- **VDJbase**: Best for human and macaque (curated, validated alleles from population studies)
- **OGRDB**: Good for mouse (community-curated novel alleles)
- **IMGT**: Species diversity (comprehensive reference database)
- **Custom**: Fill gaps (internal lab sequences for edge cases)

## Context Summary

- **Single constant change**: `GermlineManager.DEFAULT_PROVIDERS` in `manager.py`
- **First-provider-wins pattern**: List order = priority (deduplication uses first match)
- **8+ documentation locations** reference the priority order
- **2 test files** assert expected order and must be updated
- **Custom README** says "HIGHEST priority" — needs update to explain new rationale
- **pipeline.py** has hardcoded source iteration order that must match priority

---

## Plan 1: Update Core Code and Tests

```yaml
wave: 1
depends_on: []
files_modified:
  - src/sadie/germlines/manager.py
  - src/sadie/germlines/pipeline.py
  - tests/unit/germlines/test_compliance.py
  - tests/unit/germlines/test_manager.py
autonomous: true
```

**Objective:** Change the DEFAULT_PROVIDERS constant, update pipeline.py hardcoded list, and update all tests that assert the expected order.

<task id="28-1-1">
<title>Update DEFAULT_PROVIDERS Constant and All Docstrings</title>
<file>src/sadie/germlines/manager.py</file>

**Action:** Change DEFAULT_PROVIDERS from old order to new order and update ALL docstrings referencing priority order.

**Code Change (line ~50):**
```python
# BEFORE
DEFAULT_PROVIDERS = ["custom", "ogrdb", "vdjbase", "imgt"]

# AFTER
# Default provider priority order (highest to lowest):
# 1. VDJbase: Best for human/macaque - curated, validated alleles from population studies
# 2. OGRDB: Good for mouse - community-curated novel alleles
# 3. IMGT: Species diversity - comprehensive reference database
# 4. Custom: Fill gaps - internal lab sequences for edge cases
DEFAULT_PROVIDERS = ["vdjbase", "ogrdb", "imgt", "custom"]
```

**Docstring Update (line 28):**
```python
# BEFORE
    Default priority: custom > imgt > ogrdb > vdjbase

# AFTER
    Default priority: vdjbase > ogrdb > imgt > custom
```

**Docstring Update (line 41):**
```python
# BEFORE
    >>> # Default priority: custom, imgt, ogrdb

# AFTER
    >>> # Default priority: vdjbase, ogrdb, imgt, custom
```

**Docstring Update (line 60):**
```python
# BEFORE
            Default: ["custom", "imgt", "ogrdb", "vdjbase"]

# AFTER
            Default: ["vdjbase", "ogrdb", "imgt", "custom"]
```

**Verification:**
```bash
grep -n "DEFAULT_PROVIDERS\|vdjbase > ogrdb\|vdjbase, ogrdb" src/sadie/germlines/manager.py
# Should show new priority order in all locations
```
</task>

<task id="28-1-2">
<title>Update pipeline.py Hardcoded Source List</title>
<file>src/sadie/germlines/pipeline.py</file>

**Action:** Update the hardcoded source iteration order to match new priority.

**Code Change (line 168):**
```python
# BEFORE
        for source_name in ["custom", "ogrdb", "vdjbase", "imgt"]:
            if self._source_newer_than(source_name, species, latest_normalized):

# AFTER
        for source_name in ["vdjbase", "ogrdb", "imgt", "custom"]:
            if self._source_newer_than(source_name, species, latest_normalized):
```

**Verification:**
```bash
grep -n "for source_name in" src/sadie/germlines/pipeline.py
# Should show: for source_name in ["vdjbase", "ogrdb", "imgt", "custom"]
```
</task>

<task id="28-1-3">
<title>Update test_compliance.py Assertions</title>
<file>tests/unit/germlines/test_compliance.py</file>

**Action:** Update expected priority order in test assertions.

**Code Change (line 28):**
```python
# BEFORE
def test_default_priority_order(self):
    """Verify default priority is custom > ogrdb > vdjbase > imgt."""
    expected = ["custom", "ogrdb", "vdjbase", "imgt"]

# AFTER
def test_default_priority_order(self):
    """Verify default priority is vdjbase > ogrdb > imgt > custom."""
    expected = ["vdjbase", "ogrdb", "imgt", "custom"]
```

**Code Change (line 36):**
```python
# BEFORE
assert manager.provider_names == ["custom", "ogrdb", "vdjbase", "imgt"]

# AFTER
assert manager.provider_names == ["vdjbase", "ogrdb", "imgt", "custom"]
```

**Code Change (line 42):**
```python
# BEFORE
assert manager.provider_names[0] == "custom"

# AFTER
assert manager.provider_names[0] == "vdjbase"
```

**Update test docstring (line 39-40):**
```python
# BEFORE
def test_custom_provider_overrides_others(self):
    """Verify custom provider sequences take precedence."""

# AFTER
def test_vdjbase_provider_has_highest_priority(self):
    """Verify VDJbase provider sequences take precedence."""
```

**Verification:**
```bash
pytest tests/unit/germlines/test_compliance.py::TestPriorityOrder -v
```
</task>

<task id="28-1-4">
<title>Update test_manager.py Assertions</title>
<file>tests/unit/germlines/test_manager.py</file>

**Action:** Update expected priority order in test assertion.

**Code Change (lines 28-29):**
```python
# BEFORE
def test_manager_default_providers():
    """Test default provider order."""
    from sadie.germlines.manager import GermlineManager

    manager = GermlineManager()
    # Priority: custom > ogrdb > vdjbase > imgt (novel alleles prioritized)
    assert manager.provider_names == ["custom", "ogrdb", "vdjbase", "imgt"]

# AFTER
def test_manager_default_providers():
    """Test default provider order."""
    from sadie.germlines.manager import GermlineManager

    manager = GermlineManager()
    # Priority: vdjbase > ogrdb > imgt > custom (data quality prioritized)
    assert manager.provider_names == ["vdjbase", "ogrdb", "imgt", "custom"]
```

**Verification:**
```bash
pytest tests/unit/germlines/test_manager.py::test_manager_default_providers -v
```
</task>

### Plan 1 Verification

```bash
# Run all germlines unit tests
pytest tests/unit/germlines/ -v --tb=short

# Verify constant value programmatically
python -c "
from sadie.germlines.manager import GermlineManager
expected = ['vdjbase', 'ogrdb', 'imgt', 'custom']
assert GermlineManager.DEFAULT_PROVIDERS == expected, f'Got {GermlineManager.DEFAULT_PROVIDERS}'
print('SUCCESS: DEFAULT_PROVIDERS is correct')
"
```

---

## Plan 2: Update Module Docstrings

```yaml
wave: 1
depends_on: []
files_modified:
  - src/sadie/germlines/__init__.py
autonomous: true
```

**Objective:** Update module-level docstrings to reflect new priority order.

<task id="28-2-1">
<title>Update __init__.py Docstrings</title>
<file>src/sadie/germlines/__init__.py</file>

**Action:** Update all docstring references to priority order. There are TWO locations with inconsistent old order that must be fixed.

**Code Change (line 12) — CRITICAL: currently shows wrong priority order:**
```python
# BEFORE
    >>> # Simple API - uses default priority (custom > imgt > ogrdb > vdjbase)

# AFTER
    >>> # Simple API - uses default priority (vdjbase > ogrdb > imgt > custom)
```

**Code Change (line 26):**
```python
# BEFORE
- Multiple databases are used by default (custom, IMGT, OGRDB, VDJbase)

# AFTER
- Multiple databases are used by default (VDJbase, OGRDB, IMGT, custom)
```

**Code Change (line 102) — CRITICAL: currently shows wrong priority order:**
```python
# BEFORE
        Custom provider priority order. Default: ["custom", "imgt", "ogrdb", "vdjbase"]

# AFTER
        Custom provider priority order. Default: ["vdjbase", "ogrdb", "imgt", "custom"]
```

**Verification:**
```bash
grep -n "priority\|Default:" src/sadie/germlines/__init__.py | grep -v "^#"
# Should show vdjbase > ogrdb > imgt > custom in all priority references
```
</task>

---

## Plan 3: Update Documentation Files

```yaml
wave: 2
depends_on: [28-1, 28-2]
files_modified:
  - src/sadie/germlines/README.md
  - src/sadie/germlines/sources/custom/README.md
  - src/sadie/germlines/sources/vdjbase/README.md
  - src/sadie/germlines/sources/ogrdb/OGRDB_DATA.md
  - .specify/memory/constitution.md
  - docs/germlines/architecture.md
  - docs/germlines/provider-guide.md
autonomous: true
```

**Objective:** Update all documentation files to reflect new priority order and rationale.

<task id="28-3-1">
<title>Update germlines/README.md</title>
<file>src/sadie/germlines/README.md</file>

**Action:** Update priority section (around line 49).

**Code Change:**
```markdown
# BEFORE
### Priority System

Default priority order: **custom > ogrdb > vdjbase > imgt**

# AFTER
### Priority System

Default priority order: **vdjbase > ogrdb > imgt > custom**

This order prioritizes data quality:
- **VDJbase**: Curated, validated alleles (best for human/macaque)
- **OGRDB**: Community-curated novel alleles (excellent for mouse)
- **IMGT**: Comprehensive reference coverage
- **Custom**: Internal lab sequences to fill gaps
```

**Verification:**
```bash
grep -A5 "Default priority order" src/sadie/germlines/README.md
```
</task>

<task id="28-3-2">
<title>Update custom/README.md Priority Section</title>
<file>src/sadie/germlines/sources/custom/README.md</file>

**Action:** Update the priority explanation (around line 113) from "HIGHEST priority" to explain the new rationale.

**Code Change:**
```markdown
# BEFORE
## Priority System

**Custom sequences have HIGHEST priority!**

Default priority order: `custom > imgt > ogrdb`

# AFTER
## Priority System

**Custom sequences fill gaps from standard databases.**

Default priority order: `vdjbase > ogrdb > imgt > custom`

Custom sequences are applied LAST, after validated databases. This ensures:
- **Validated data first**: VDJbase and OGRDB provide curated, community-validated alleles
- **Gap filling**: Custom sequences add lab-specific alleles not in public databases
- **Override capability**: To override a standard allele, use the same gene name in custom

To explicitly override a VDJbase/OGRDB gene, add a custom sequence with the identical gene name.
```

**Verification:**
```bash
grep -A10 "## Priority System" src/sadie/germlines/sources/custom/README.md
```
</task>

<task id="28-3-3">
<title>Update docs/germlines/architecture.md</title>
<file>docs/germlines/architecture.md</file>

**Action:** Update priority order references AND code sample. There is a code sample at line 233 showing reverse iteration order that must be updated.

**Code Changes:**
1. Update any diagrams or text showing priority order
2. Change `custom > ogrdb > vdjbase > imgt` to `vdjbase > ogrdb > imgt > custom`
3. Update rationale explanations
4. **CRITICAL: Update code sample (line 233):**
```python
# BEFORE (line 233)
    for provider in ['imgt', 'vdjbase', 'ogrdb', 'custom']:

# AFTER
    for provider in ['custom', 'imgt', 'ogrdb', 'vdjbase']:
```
   *(Note: This is REVERSE order since line 232 says "Process in reverse priority order (lowest first)")*

**Verification:**
```bash
grep -n "priority\|for provider in" docs/germlines/architecture.md
```
</task>

<task id="28-3-4">
<title>Update docs/germlines/provider-guide.md</title>
<file>docs/germlines/provider-guide.md</file>

**Action:** Update priority documentation (around line 170).

**Code Changes:**
1. Update default priority order list
2. Add rationale for new order
3. Update examples if any show old order

**Verification:**
```bash
grep -n "priority" docs/germlines/provider-guide.md
```
</task>

<task id="28-3-5">
<title>Update VDJbase README.md Priority Section</title>
<file>src/sadie/germlines/sources/vdjbase/README.md</file>

**Action:** Update the priority explanation (line 103) from "lower priority" to reflect VDJbase now having HIGHEST priority.

**Code Change (line 103):**
```markdown
# BEFORE
3. **Priority**: By default, VDJbase has lower priority than OGRDB and IMGT. Configure priority in GermlineManager if you want VDJbase sequences to take precedence.

# AFTER
3. **Priority**: By default, VDJbase has HIGHEST priority (`vdjbase > ogrdb > imgt > custom`). VDJbase provides curated, validated alleles from population studies, making it the preferred source for human and macaque data.
```

**Verification:**
```bash
grep -n "Priority" src/sadie/germlines/sources/vdjbase/README.md
# Should show: VDJbase has HIGHEST priority
```
</task>

<task id="28-3-6">
<title>Update OGRDB_DATA.md Priority Section</title>
<file>src/sadie/germlines/sources/ogrdb/OGRDB_DATA.md</file>

**Action:** Update the priority order reference (line 106).

**Code Change (line 106):**
```markdown
# BEFORE
Default priority order: `custom > ogrdb > vdjbase > imgt`

# AFTER
Default priority order: `vdjbase > ogrdb > imgt > custom`
```

**Code Change (update surrounding context, lines 104-112):**
```markdown
# BEFORE
## Priority with IMGT

Default priority order: `custom > ogrdb > vdjbase > imgt`

- If a gene exists in both OGRDB and IMGT → OGRDB version is used (higher priority)

# AFTER
## Priority with Other Sources

Default priority order: `vdjbase > ogrdb > imgt > custom`

- VDJbase has highest priority for human/macaque data
- If a gene exists in both OGRDB and IMGT → OGRDB version is used (OGRDB > IMGT)
```

**Verification:**
```bash
grep -A5 "## Priority" src/sadie/germlines/sources/ogrdb/OGRDB_DATA.md
# Should show: vdjbase > ogrdb > imgt > custom
```
</task>

<task id="28-3-7">
<title>Update Constitution Priority (NON-NEGOTIABLE)</title>
<file>.specify/memory/constitution.md</file>

**Action:** Update the NON-NEGOTIABLE priority order in Principle II (line 11).

**Code Change (line 11):**
```markdown
# BEFORE
The system MUST support explicit priority ordering for conflict resolution when multiple providers supply the same gene. Default priority: `custom > ogrdb > vdjbase > imgt`. When duplicate gene names exist with different sequences, the highest-priority provider's sequence MUST be kept; lower-priority duplicates dropped with WARNING log.

# AFTER
The system MUST support explicit priority ordering for conflict resolution when multiple providers supply the same gene. Default priority: `vdjbase > ogrdb > imgt > custom`. When duplicate gene names exist with different sequences, the highest-priority provider's sequence MUST be kept; lower-priority duplicates dropped with WARNING log.
```

**Verification:**
```bash
grep "Default priority:" .specify/memory/constitution.md
# Should show: Default priority: `vdjbase > ogrdb > imgt > custom`
```
</task>

---

## Dependencies Between Plans

```
Plan 1 (Core Code + Tests) ←─┐
                              ├── Wave 1 (parallel)
Plan 2 (Module Docstrings) ←─┘
         ↓
Plan 3 (Documentation) ← Wave 2 (depends on Wave 1)
```

Wave 1 plans can execute in parallel. Plan 3 (documentation) runs after code changes are verified.

---

## Success Criteria (Phase Level)

1. ✅ `GermlineManager.DEFAULT_PROVIDERS == ["vdjbase", "ogrdb", "imgt", "custom"]`
2. ✅ `pipeline.py` source iteration matches new priority order
3. ✅ All tests in `tests/unit/germlines/` pass
4. ✅ Human/macaque queries prefer VDJbase alleles when available
5. ✅ Mouse queries get OGRDB alleles (VDJbase has limited mouse data)
6. ✅ Fallback chain works correctly when preferred source lacks data
7. ✅ All documentation references updated consistently

---

## must_haves

Derived from phase goal using goal-backward methodology:

| ID | Requirement | Verification |
|----|-------------|--------------|
| MH-1 | DEFAULT_PROVIDERS constant is `["vdjbase", "ogrdb", "imgt", "custom"]` | `python -c "from sadie.germlines.manager import GermlineManager; assert GermlineManager.DEFAULT_PROVIDERS == ['vdjbase', 'ogrdb', 'imgt', 'custom']"` |
| MH-2 | Manager initializes with new default order | `pytest tests/unit/germlines/test_manager.py::test_manager_default_providers -v` |
| MH-3 | Priority tests pass with new expected values | `pytest tests/unit/germlines/test_compliance.py::TestPriorityOrder -v` |
| MH-4 | No regressions in germlines tests | `pytest tests/unit/germlines/ -v` |
| MH-5 | Docstrings in manager.py reflect new priority | `grep -c "vdjbase > ogrdb" src/sadie/germlines/manager.py` (should be ≥1) |
| MH-6 | Docstrings in __init__.py reflect new priority | `grep "vdjbase > ogrdb\|vdjbase, ogrdb" src/sadie/germlines/__init__.py \| wc -l` (should be ≥2) |
| MH-7 | pipeline.py has consistent source order | `grep "for source_name in" src/sadie/germlines/pipeline.py \| grep -q "vdjbase.*ogrdb.*imgt.*custom"` |
| MH-8 | architecture.md code sample updated | `grep "for provider in" docs/germlines/architecture.md \| grep -q "custom.*imgt.*ogrdb.*vdjbase"` |
| MH-9 | VDJbase README reflects HIGHEST priority | `grep "HIGHEST priority" src/sadie/germlines/sources/vdjbase/README.md` |
| MH-10 | OGRDB_DATA.md shows new priority order | `grep "vdjbase > ogrdb > imgt > custom" src/sadie/germlines/sources/ogrdb/OGRDB_DATA.md` |
| MH-11 | Constitution updated (NON-NEGOTIABLE) | `grep "vdjbase > ogrdb > imgt > custom" .specify/memory/constitution.md` |

---

## Files Modified Summary

### Code Changes (4 files)
| File | Change Type |
|------|-------------|
| `src/sadie/germlines/manager.py` | Constant + docstrings (3 locations) |
| `src/sadie/germlines/pipeline.py` | Hardcoded source list (line 168) |
| `tests/unit/germlines/test_compliance.py` | Test assertions |
| `tests/unit/germlines/test_manager.py` | Test assertions |

### Docstring Updates (1 file)
| File | Change Type |
|------|-------------|
| `src/sadie/germlines/__init__.py` | Module docstrings (lines 12, 26, 102) |

### Documentation Updates (7 files)
| File | Change Type |
|------|-------------|
| `src/sadie/germlines/README.md` | Priority section |
| `src/sadie/germlines/sources/custom/README.md` | Priority explanation |
| `src/sadie/germlines/sources/vdjbase/README.md` | Priority section (line 103) — "HIGHEST priority" |
| `src/sadie/germlines/sources/ogrdb/OGRDB_DATA.md` | Priority section (line 106) |
| `.specify/memory/constitution.md` | Principle II NON-NEGOTIABLE priority (line 11) |
| `docs/germlines/architecture.md` | Architecture docs + code sample (line 233) |
| `docs/germlines/provider-guide.md` | Provider guide |

---

## Verification Checklist

```bash
# Core verification
[ ] GermlineManager.DEFAULT_PROVIDERS == ["vdjbase", "ogrdb", "imgt", "custom"]
[ ] pipeline.py source order == ["vdjbase", "ogrdb", "imgt", "custom"]
[ ] pytest tests/unit/germlines/test_compliance.py::TestPriorityOrder passes
[ ] pytest tests/unit/germlines/test_manager.py::test_manager_default_providers passes

# Regression check
[ ] pytest tests/unit/germlines/ passes (all tests)

# Documentation consistency
[ ] manager.py docstrings say "vdjbase > ogrdb > imgt > custom" (3 locations)
[ ] __init__.py docstrings updated (lines 12, 26, 102)
[ ] README.md priority section updated
[ ] custom/README.md priority explanation updated
[ ] vdjbase/README.md says "HIGHEST priority" (line 103)
[ ] OGRDB_DATA.md says "vdjbase > ogrdb > imgt > custom" (line 106)
[ ] constitution.md Principle II updated (line 11)
[ ] docs/germlines/architecture.md updated (including line 233 code sample)
[ ] docs/germlines/provider-guide.md updated
```

---

## Risk Assessment

**Risk Level**: LOW

- Simple constant change with clear scope
- No new dependencies
- No architectural changes
- Existing test patterns can verify the change
- All locations identified through grep search

**Mitigation**: Run full germlines test suite after each file change to catch assertion mismatches early.
