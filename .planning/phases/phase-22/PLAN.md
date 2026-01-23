# Phase 22: Runtime Usage

## Overview

**Goal:** Enable Airr to use prebuilt databases directly, bypassing runtime gene lookup

**Dependencies:** Phase 21 (Build CLI) creates database structure with:
- `Ig/blastdb/{name}/` - BLAST database files
- `Ig/internal_data/{name}/` - Internal annotation files (`.ndm.imgt`)
- `aux_db/` - Auxiliary data (`{name}_gl.aux`)
- `.references_dataframe.csv.gz` - References manifest

**Impact:** Reduce startup from >1s (network/file lookup) to <100ms (prebuilt validation)

---

## Problem Analysis

### Current Behavior

`Airr.__init__()` currently:
1. Accepts `reference_name: str` for species lookup
2. Optionally accepts `references: References` for custom databases
3. Uses `GermlineData` to resolve database paths from germlines module or G3

`GermlineData.__init__()` currently:
1. Accepts optional `database_dir` but only uses it to set `base_dir`
2. Does NOT set up gene directory paths when custom `database_dir` is provided
3. Falls through to germlines module path resolution

### Required Behavior

When `Airr(database="./db")` is called:
1. Skip all germlines/G3 lookups entirely
2. Validate prebuilt structure exists
3. Use database paths directly from prebuilt location
4. Never make network calls or file discovery searches

---

## Tasks

### Task 1: Add Database Structure Validation

**File:** `src/sadie/airr/igblast/germline.py`

**Changes:**

Add validation function before `GermlineData` class:

```python
def validate_prebuilt_database(database_path: Path, name: str) -> dict:
    """Validate prebuilt database structure and return paths.
    
    Expected structure:
        database_path/
        ├── Ig/
        │   ├── blastdb/{name}/{name}_V, {name}_D, {name}_J files
        │   └── internal_data/{name}/{name}.ndm.imgt
        ├── aux_db/{name}_gl.aux
        └── .references_dataframe.csv.gz (optional)
    
    Parameters
    ----------
    database_path : Path
        Root path to prebuilt database
    name : str
        Reference name (species)
        
    Returns
    -------
    dict
        Validated paths: blast_dir, v_gene_dir, d_gene_dir, j_gene_dir, 
        c_gene_dir, aux_path, igdata
        
    Raises
    ------
    FileNotFoundError
        If required structure is missing
    """
    db = Path(database_path)
    errors = []
    
    # Required directories
    ig_dir = db / "Ig"
    internal_data = ig_dir / "internal_data" / name
    blastdb = ig_dir / "blastdb" / name
    aux_db = db / "aux_db"
    
    if not ig_dir.exists():
        errors.append(f"Missing Ig/ directory at {ig_dir}")
    if not internal_data.exists():
        errors.append(f"Missing internal_data/{name}/ at {internal_data}")
    if not blastdb.exists():
        errors.append(f"Missing blastdb/{name}/ at {blastdb}")
    if not aux_db.exists():
        errors.append(f"Missing aux_db/ directory at {aux_db}")
    
    aux_path = aux_db / f"{name}_gl.aux"
    if not aux_path.exists():
        errors.append(f"Missing auxiliary file {name}_gl.aux at {aux_path}")
    
    if errors:
        raise FileNotFoundError(
            f"Invalid prebuilt database at {database_path}:\n" + 
            "\n".join(f"  - {e}" for e in errors)
        )
    
    # Build paths
    blast_prefix = blastdb / f"{name}_"
    
    return {
        "base_dir": db,
        "blast_dir": blast_prefix,
        "v_gene_dir": Path(str(blast_prefix) + "V"),
        "d_gene_dir": Path(str(blast_prefix) + "D"),
        "j_gene_dir": Path(str(blast_prefix) + "J"),
        "c_gene_dir": Path(str(blast_prefix) + "C"),
        "aux_path": aux_path,
        "igdata": ig_dir,
    }
```

**Success Criteria:**
- Function raises `FileNotFoundError` with clear message listing missing components
- Returns valid path dict when structure is complete
- Validates V/D/J blast database prefixes exist

---

### Task 2: Update GermlineData for Prebuilt Support

**File:** `src/sadie/airr/igblast/germline.py`

**Changes:**

Modify `GermlineData.__init__()` to handle prebuilt database paths:

```python
def __init__(
    self,
    name: str,
    receptor: str = "Ig",
    database_dir: Optional[str | Path] = None,
    scheme: str = "imgt",
    prebuilt: bool = False,  # NEW PARAMETER
):
    """
    Parameters
    ----------
    name : str
        The species of interest, e.g. human
    receptor : str, optional
        the receptor type, by default "Ig"
    database_dir : Optional[str | Path]
        Custom database directory path
    scheme : str, optional
        Numbering scheme, by default "imgt"
    prebuilt : bool, optional
        If True, database_dir is a prebuilt database from `sadie reference build`.
        Validates structure and uses paths directly without germlines/G3 lookup.
        By default False.
    """
    self.name = name

    # NEW: Handle prebuilt database path
    if prebuilt and database_dir:
        paths = validate_prebuilt_database(Path(database_dir), name)
        self._base_dir = paths["base_dir"]
        self._blast_dir = paths["blast_dir"]
        self._v_gene_dir = paths["v_gene_dir"]
        self._d_gene_dir = paths["d_gene_dir"]
        self._j_gene_dir = paths["j_gene_dir"]
        self._c_gene_dir = paths["c_gene_dir"]
        self._aux_path = paths["aux_path"]
        self._igdata = paths["igdata"]
        return  # Skip all other initialization
    
    # Existing logic follows...
```

**Key Change:** When `prebuilt=True`:
1. Call `validate_prebuilt_database()` to validate and get paths
2. Set all path attributes directly (bypassing setters that may fail validation)
3. Return early - skip germlines module and legacy G3 lookups entirely

**Success Criteria:**
- `GermlineData(name, database_dir=path, prebuilt=True)` works without germlines lookup
- Standard `GermlineData(name)` behavior unchanged
- Invalid prebuilt paths raise clear error

---

### Task 3: Add `database` Parameter to Airr

**File:** `src/sadie/airr/airr.py`

**Changes:**

1. Add `database` parameter to `__init__()`:

```python
def __init__(
    self,
    reference_name: str,
    igblast_exe: Path | str = "",
    adaptable: bool = False,
    # ... existing parameters ...
    coerce: bool = False,
    database: Optional[Path | str] = None,  # NEW PARAMETER
):
    """Airr constructor
    
    Parameters
    ----------
    reference_name : str | Reference
        the reference name to run annotate against. ...
    # ... existing docstring ...
    database : Optional[Path | str]
        Path to prebuilt database from `sadie reference build`.
        When provided, uses database directly without germlines/G3 lookup.
        Expected structure: Ig/blastdb/, Ig/internal_data/, aux_db/.
        By default None (uses germlines module or G3).
    """
```

2. Add initialization logic after existing `self.references` handling:

```python
# Store database path for recursive Airr calls
self._database_path = database

# ... existing _create_temp and references handling ...

# NEW: Handle prebuilt database path
if database:
    database_path = Path(database)
    if not database_path.exists():
        raise FileNotFoundError(f"Database path not found: {database_path}")
    
    # Use prebuilt database - skip all germlines/G3 lookup
    self.germline_data = GermlineData(
        reference_name, 
        receptor, 
        database_path, 
        scheme,
        prebuilt=True
    )
elif isinstance(references, References):
    # Existing custom References handling...
else:
    # Existing germlines/G3 lookup...
```

3. Pass `database` to recursive `Airr` calls in `run_fasta()`:

```python
adaptable_api = Airr(
    self.name,
    igblast_exe=self.executable,
    # ... existing parameters ...
    coerce=self.coerce,
    database=self._database_path,  # NEW: Pass through database path
)
```

**Success Criteria:**
- `Airr(database="./db")` works without network calls
- `Airr("human")` behavior unchanged
- Database path passed to recursive penalty adaptation calls

---

## Execution Order

```
Wave 1: Task 1 (validation function - no dependencies)
Wave 2: Task 2 (GermlineData update - depends on Task 1)
Wave 3: Task 3 (Airr update - depends on Task 2)
```

---

## Verification Plan

### Unit Tests

**File:** `tests/unit/airr/test_database_parameter.py`

```python
import pytest
from pathlib import Path
from sadie.airr import Airr
from sadie.airr.igblast.germline import GermlineData, validate_prebuilt_database

class TestValidatePrebuiltDatabase:
    """Test database structure validation."""
    
    def test_missing_ig_dir_raises(self, tmp_path):
        """Missing Ig/ directory raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="Missing Ig/"):
            validate_prebuilt_database(tmp_path, "human")
    
    def test_missing_blastdb_raises(self, tmp_path):
        """Missing blastdb raises FileNotFoundError."""
        (tmp_path / "Ig" / "internal_data" / "human").mkdir(parents=True)
        (tmp_path / "aux_db").mkdir()
        (tmp_path / "aux_db" / "human_gl.aux").touch()
        with pytest.raises(FileNotFoundError, match="Missing blastdb"):
            validate_prebuilt_database(tmp_path, "human")
    
    def test_missing_aux_raises(self, tmp_path):
        """Missing aux_db raises FileNotFoundError."""
        (tmp_path / "Ig" / "internal_data" / "human").mkdir(parents=True)
        (tmp_path / "Ig" / "blastdb" / "human").mkdir(parents=True)
        with pytest.raises(FileNotFoundError, match="Missing aux_db"):
            validate_prebuilt_database(tmp_path, "human")


class TestGermlineDataPrebuilt:
    """Test GermlineData with prebuilt=True."""
    
    def test_prebuilt_skips_germlines_lookup(self, valid_prebuilt_db):
        """Prebuilt database skips germlines module lookup."""
        gd = GermlineData("human", database_dir=valid_prebuilt_db, prebuilt=True)
        assert gd.base_dir == valid_prebuilt_db
        assert "human_V" in str(gd.v_gene_dir)


class TestAirrDatabaseParameter:
    """Test Airr(database=...) parameter."""
    
    def test_database_path_not_found_raises(self):
        """Non-existent database path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            Airr("human", database="/nonexistent/path")
    
    def test_database_parameter_uses_prebuilt(self, valid_prebuilt_db):
        """Airr with database= uses prebuilt paths."""
        # This test requires a fully valid prebuilt database
        # with actual BLAST database files
        pass  # Integration test
```

### Integration Tests

**File:** `tests/integration/test_prebuilt_database.py`

```python
import pytest
from sadie.airr import Airr

class TestPrebuiltDatabaseIntegration:
    """Integration tests for prebuilt database usage."""
    
    @pytest.fixture
    def prebuilt_human_db(self, tmp_path):
        """Build a real prebuilt database for testing."""
        from sadie.reference import References
        refs = References.from_yaml()
        refs.make_airr_database(tmp_path)
        return tmp_path
    
    def test_annotation_identical_prebuilt_vs_runtime(self, prebuilt_human_db):
        """Results identical with prebuilt vs runtime database."""
        seq = "CAGGTGCAGCTGGTGGAGTCTGGGGGAGGC..."  # Test sequence
        
        # Runtime lookup
        airr_runtime = Airr("human")
        result_runtime = airr_runtime.run_single("test", seq)
        
        # Prebuilt database
        airr_prebuilt = Airr("human", database=prebuilt_human_db)
        result_prebuilt = airr_prebuilt.run_single("test", seq)
        
        # Compare key columns
        assert result_runtime["v_call"].iloc[0] == result_prebuilt["v_call"].iloc[0]
        assert result_runtime["j_call"].iloc[0] == result_prebuilt["j_call"].iloc[0]
    
    def test_no_network_calls_with_prebuilt(self, prebuilt_human_db, mocker):
        """No network/file discovery with prebuilt database."""
        # Mock germlines lookup to ensure it's not called
        mock_germlines = mocker.patch(
            "sadie.airr.igblast.germline._use_germlines_module"
        )
        
        Airr("human", database=prebuilt_human_db)
        
        mock_germlines.assert_not_called()
```

### Performance Test

```python
import time

def test_startup_performance_prebuilt(prebuilt_human_db):
    """Prebuilt database startup under 100ms."""
    start = time.perf_counter()
    Airr("human", database=prebuilt_human_db)
    elapsed = time.perf_counter() - start
    
    assert elapsed < 0.1, f"Startup took {elapsed:.3f}s, expected <0.1s"
```

---

## Success Criteria Summary

| Requirement | Implementation | Verification |
|-------------|----------------|--------------|
| RUN-01: Add `database` param | Task 3: Airr.__init__() | Unit test + docstring |
| RUN-02: Skip germlines lookup | Task 2+3: prebuilt=True | Mock test |
| RUN-03: Validate structure | Task 1: validate_prebuilt_database() | Error message test |

### Acceptance Criteria

1. ✅ `Airr(database="./db")` uses prebuilt database instead of default
2. ✅ No network calls or germlines lookups when database path provided
3. ✅ Clear error if database path missing required structure
4. ✅ Annotation results identical whether using prebuilt or runtime-built
5. ✅ Performance: <100ms startup with prebuilt vs >1s with runtime lookup

---

## Files Created/Modified

| Action | File |
|--------|------|
| MODIFY | `src/sadie/airr/igblast/germline.py` |
| MODIFY | `src/sadie/airr/airr.py` |
| CREATE | `tests/unit/airr/test_database_parameter.py` |
| CREATE | `tests/integration/test_prebuilt_database.py` |

---

## Risks

| Risk | Mitigation |
|------|------------|
| Phase 21 structure differs | Verify Phase 21 output structure first |
| Recursive Airr calls miss database | Pass `_database_path` to all recursive calls |
| Setter validation fails on prebuilt | Use direct `_attribute` assignment |

---

*Plan created: 2026-01-23*
