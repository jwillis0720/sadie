# Phase 26: Add AIRR Package Dependency - Research

## Summary

This phase adds the official AIRR-C Python package (`airr`) from PyPI as a production dependency. Currently, it's only a dev dependency, causing `test_write_and_check_airr` to be skipped when running tests without dev dependencies installed.

**Primary Recommendation:** Move `airr = "^1.5.0"` from `[tool.poetry.group.dev.dependencies]` to `[tool.poetry.dependencies]` in pyproject.toml.

## Standard Stack

| Library | Version | Purpose | Source |
|---------|---------|---------|--------|
| airr | ^1.5.0 | AIRR-C standards validation | PyPI |

**Note:** Version `^1.5.0` (currently 1.5.1 on PyPI) is already specified in dev dependencies. Keep the same version constraint when moving to production dependencies.

### Dependency Chain

The `airr` package depends on:
- `pandas>=0.24.0` ✓ (already in SADIE dependencies)
- `pyyaml>=3.12` ✓ (already in SADIE as `PyYAML = "^6.0"`)
- `yamlordereddictloader>=0.4.0` (will be pulled in transitively)
- `setuptools>=2.0` ✓ (standard Python dependency)

**No new direct dependencies** need to be added—the `airr` package is lightweight and its dependencies overlap with existing SADIE dependencies.

## Architecture Patterns

### Naming Disambiguation

**CRITICAL:** The codebase has TWO "airr" entities:
1. **`sadie.airr`** - SADIE's internal antibody analysis module (major functionality)
2. **`airr`** (PyPI) - Official AIRR-C standards library (used only for validation)

They are completely separate. The PyPI `airr` package is only imported in tests:
```python
# In test_airr.py - line 603
airr = pytest.importorskip("airr", reason="airr package not installed")
```

### Test Pattern

The test uses the official airr package for AIRR-C compliance validation:
```python
# Writing SADIE's AIRR output
airr_table.to_airr(output_file)

# Validating with official airr package
d = airr.load_rearrangement(output_file, debug=True, validate=True)

# Catching validation errors
with pytest.raises(airr.schema.ValidationError):
    airr.load_rearrangement(bad_file, debug=True, validate=True)
```

## Don't Hand-Roll

| Problem | Use Instead | Why |
|---------|-------------|-----|
| AIRR TSV validation | `airr.load_rearrangement(..., validate=True)` | Official AIRR-C reference implementation |
| AIRR schema errors | `airr.schema.ValidationError` | Standard exception class |

**Never** implement custom AIRR validation logic—always use the official `airr` package.

## Common Pitfalls

### 1. Forgetting to Remove `pytest.importorskip`
**Problem:** Adding the dependency but leaving `pytest.importorskip("airr")` in the test.
**Consequence:** Test still works but the skip message becomes misleading.
**Solution:** Remove the entire `pytest.importorskip` line after adding the dependency.

### 2. Namespace Confusion
**Problem:** Confusing `sadie.airr` (internal module) with `airr` (PyPI package).
**Consequence:** Could accidentally import the wrong module or create circular imports.
**Solution:** The test imports `airr` (PyPI) AFTER all `sadie.airr` imports. Maintain this pattern.

### 3. Version Mismatch
**Problem:** Using an older version that lacks `load_rearrangement` or `ValidationError`.
**Consequence:** Test failures or AttributeError.
**Solution:** Use `^1.5.0` which includes all required functionality.

### 4. Poetry Lock File Not Updated
**Problem:** Changing pyproject.toml without running `poetry lock`.
**Consequence:** CI/CD or fresh installs may not pick up the new dependency.
**Solution:** Run `poetry lock --no-update` after modifying pyproject.toml (or `poetry lock` if full update is desired).

## Code Examples

### 1. Moving Dependency in pyproject.toml

**Before:**
```toml
[tool.poetry.dependencies]
# ... other dependencies ...

[tool.poetry.group.dev.dependencies]
airr = "^1.5.0"
# ... other dev dependencies ...
```

**After:**
```toml
[tool.poetry.dependencies]
# ... other dependencies ...
airr = "^1.5.0"

[tool.poetry.group.dev.dependencies]
# airr moved to main dependencies
# ... other dev dependencies ...
```

### 2. Removing pytest.importorskip in Test

**Before:**
```python
def test_write_and_check_airr(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    """Check that the offical airr can validate our airr tables"""
    airr = pytest.importorskip("airr", reason="airr package not installed")
    # ... rest of test
```

**After:**
```python
import airr  # Add at top of file with other imports

def test_write_and_check_airr(tmp_path_factory: pytest.TempPathFactory, fixture_setup: SadieFixture) -> None:
    """Check that the offical airr can validate our airr tables"""
    # pytest.importorskip removed - airr is now a required dependency
    # ... rest of test
```

### 3. Import Placement

Place the import at the top of the test file with other external imports:
```python
# ... other imports ...
import airr  # AIRR-C standards validation library
from sadie.airr import Airr, AirrSeries, AirrTable, GermlineData, LinkedAirrTable
# ... rest of imports ...
```

## Verification Checklist

- [ ] `airr` appears in `[tool.poetry.dependencies]` section
- [ ] `airr` removed from `[tool.poetry.group.dev.dependencies]` section
- [ ] `poetry lock` runs without errors
- [ ] `import airr` at top of test file (not inside function)
- [ ] `pytest.importorskip("airr")` line removed from test
- [ ] `pytest tests/unit/airr/test_airr.py::test_write_and_check_airr -v` passes (not skipped)
- [ ] Full test suite passes

## Sources

| Source | Confidence | Notes |
|--------|------------|-------|
| [PyPI - airr](https://pypi.org/project/airr/) | HIGH | Official package page, version 1.5.1 |
| [AIRR Python Docs](https://docs.airr-community.org/en/stable/packages/airr-python/overview.html) | HIGH | Official documentation |
| [GitHub - airr-standards](https://github.com/airr-community/airr-standards/tree/master/lang/python) | HIGH | Official source code |
| [Bioconda - airr](https://bioconda.github.io/recipes/airr/README.html) | MEDIUM | Alternative installation method |
| Local pyproject.toml analysis | HIGH | Direct codebase inspection |
| Local test file analysis | HIGH | Direct codebase inspection |

## Risk Assessment

**Overall Risk: LOW**

- This is a simple dependency relocation (dev → main)
- No code changes beyond import statement
- The package is already tested in dev environment
- Dependencies are compatible with existing stack
- No breaking changes expected

**Potential Issues:**
- Package size increase: ~83KB (negligible)
- Additional transitive dependency: `yamlordereddictloader` (small, well-maintained)
