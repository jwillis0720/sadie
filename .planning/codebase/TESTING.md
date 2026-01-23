# G3 Testing Patterns for Segment Position Discovery

## Overview

This document describes how SADIE tests G3 backend segment position discovery, particularly for J gene segments and CDR3/junction boundary calculations.

---

## 1. Test Organization

### Directory Structure
```
tests/
├── unit/
│   ├── airr/
│   │   ├── test_cdr3.py           # CDR3 field population tests
│   │   ├── test_igblast.py        # IgBLAST setup and execution
│   │   └── test_methods.py        # AIRR methods
│   ├── germlines/
│   │   ├── test_airr_integration.py  # AIRR + germlines integration
│   │   └── test_compliance.py     # AIRR compliance
│   └── renumbering/
│       └── test_g3.py             # G3 Stockholm alignment
├── integration/
│   └── airr/
│       └── test_airr_integration.py
└── conftest.py                    # Shared fixtures

src/sadie/germlines/tests/
├── test_g3_regression.py          # G3 gene presence regression
├── test_normalized_output.py      # Normalized output structure
└── test_germline_data_legacy.py   # Legacy API compatibility
```

---

## 2. CDR3 Position Validation Tests

### Test File: `tests/unit/airr/test_cdr3.py`

### Test Class: `TestCDR3Population`

#### Test: CDR3 Field Not None
```python
def test_cdr3_field_not_none(self, cdr3_known_bugs_fasta: Path) -> None:
    """Test that CDR3 field is not None for sequences.
    Note: Seq5 is expected to have no CDR3 due to truncated J region.
    """
    airr_api = Airr("human")
    airr_table = airr_api.run_fasta(str(cdr3_known_bugs_fasta))
    
    assert not airr_table.empty, "AIRR annotation returned empty table"
    
    cdr3_nulls = airr_table["cdr3"].isna()
    sequences_with_null_cdr3 = airr_table[cdr3_nulls]["sequence_id"].tolist()
    expected_null_cdr3 = ["Seq5"]  # Known exception
    unexpected_null_cdr3 = [seq for seq in sequences_with_null_cdr3 
                           if seq not in expected_null_cdr3]
    
    assert not unexpected_null_cdr3
```

#### Test: CDR3 Start/End Position Validity
```python
def test_cdr3_start_end_positions(self, cdr3_known_bugs_fasta: Path) -> None:
    """Test that CDR3 start and end positions are valid when CDR3 is present."""
    airr_table = airr_api.run_fasta(str(cdr3_known_bugs_fasta))
    sequences_with_cdr3 = airr_table[airr_table["cdr3"].notna()]
    
    for idx, row in sequences_with_cdr3.iterrows():
        # Position not null
        assert pd.notna(row["cdr3_start"])
        assert pd.notna(row["cdr3_end"])
        
        # Positions are numeric
        assert isinstance(row["cdr3_start"], (int, float))
        assert isinstance(row["cdr3_end"], (int, float))
        
        # End > Start
        assert row["cdr3_end"] > row["cdr3_start"]
        
        # Length consistency
        expected_length = int(row["cdr3_end"] - row["cdr3_start"] + 1)
        actual_length = len(row["cdr3"])
        assert actual_length == expected_length
```

#### Test: CDR3 Extraction Accuracy
```python
def test_cdr3_extraction_from_sequence(self, cdr3_known_bugs_fasta: Path) -> None:
    """Test that CDR3 sequence is correctly extracted from the full sequence."""
    for idx, row in sequences_with_cdr3.iterrows():
        if pd.notna(row["cdr3_start"]) and pd.notna(row["cdr3_end"]):
            start = int(row["cdr3_start"]) - 1  # Convert to 0-based
            end = int(row["cdr3_end"])          # End is inclusive
            
            extracted_cdr3 = row["sequence"][start:end]
            cdr3_no_gaps = row["cdr3"].replace("-", "")
            
            assert extracted_cdr3 == cdr3_no_gaps
```

---

## 3. IgBLAST Setup and Aux File Tests

### Test File: `tests/unit/airr/test_igblast.py`

#### Test: Aux Path Validation
```python
def test_antibody_igblast_setup() -> None:
    """Test manual IgBLAST setup including aux files."""
    for name in ["human", "mouse", "rat", "dog"]:
        aux_ref = os.path.join(germline_ref, "aux_db/imgt/")
        aux_ref = os.path.join(aux_ref, f"{name}_gl.aux")
        
        ig_blast.aux_path = aux_ref  # Will raise if invalid
        ig_blast.pre_check()
```

#### Test: Bad Aux Path Raises Exception
```python
def test_aux_path_exceptions() -> None:
    with pytest.raises(airr_exceptions.BadIgBLASTArgument):
        ig_blast.aux_path = "nonexistent_aux_file"
```

---

## 4. G3 Stockholm Alignment Tests

### Test File: `tests/unit/renumbering/test_g3.py`

#### Test: Stockholm Alignment Consistency
```python
def test_stockholm_pairs(fixture_setup):
    """Test G3 Stockholm alignments match numbering module."""
    g3 = G3()
    
    species_list = ["human", "mouse", "rat", "rabbit"]
    chains = ["H", "L", "K"]
    
    for species in species_list:
        for chain in chains:
            stockholm_pairs = g3.get_stockholm_pairs(species=species, chain=chain)
            
            for name, align in stockholm_pairs:
                # Check end of V to J alignment is conserved
                assert numbering_align[-5:] == g3_align[-5:]
                assert numbering_align[-35:].count("-") == g3_align[-35:].count("-")
```

---

## 5. G3 Regression Tests

### Test File: `src/sadie/germlines/tests/test_g3_regression.py`

#### Test: Required Genes Present
```python
def test_germlines_returns_ighj4(germline_manager):
    """Verify IGHJ4 is returned by germlines manager."""
    genes = germline_manager.get_genes("human", "J", "H")
    gene_names = [g.name for g in genes]
    assert any("IGHJ4" in name for name in gene_names)

def test_sequence_identity_ighj4(germline_manager):
    """Verify IGHJ4*01 sequence matches expected."""
    genes = germline_manager.get_genes("human", "J", "H")
    ighj4 = next((g for g in genes if "IGHJ4*01" in g.name), None)
    
    assert ighj4 is not None
    expected = "actactttgactactggggccaaggaaccctggtcaccgtctcctcag"
    actual = ighj4.sequence.replace(".", "").lower()
    assert actual == expected
```

---

## 6. AIRR Integration Tests

### Test File: `tests/unit/germlines/test_airr_integration.py`

#### Test: Gene Call Presence
```python
def test_airr_annotation_with_germlines(self, monkeypatch, tmp_path):
    """Test AIRR annotation using germlines backend."""
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
    airr = Airr(reference_name="human")
    result = airr.run_single("test_seq", self.TEST_SEQ)
    
    # Verify gene calls were made
    assert "v_call" in result.columns
    assert "d_call" in result.columns
    assert "j_call" in result.columns
    
    v_call = result["v_call"].iloc[0]
    assert v_call is not None and v_call != ""
```

#### Test: Offline Operation
```python
def test_offline_operation(self, monkeypatch, tmp_path):
    """Test AIRR annotation works offline with germlines module."""
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
    
    gd = GermlineData("human")
    assert gd.base_dir.exists()
    
    result = airr.run_single("test_seq", self.TEST_SEQ)
    assert not result.empty
```

---

## 7. Aux File Validation Tests

### Test File: `src/sadie/germlines/tests/test_normalized_output.py`

#### Test: Aux File Structure
```python
def test_aux_file_has_correct_columns():
    """Verify aux file has exactly 5 tab-separated columns."""
    aux_file = Path("igblast/aux_db/human_gl.aux")
    lines = aux_file.read_text().strip().split("\n")
    
    for line in lines:
        fields = line.split("\t")
        assert len(fields) == 5, f"Expected 5 columns, got {len(fields)}"
```

---

## 8. Test Fixtures

### Test Fixture: SadieFixture
Location: `tests/conftest.py`

```python
class SadieFixture(AirrSequences, AirrTables, ReferenceFixtures, NumberingFixtures):
    def __init__(self, tmp_path_factory):
        # Initialize test data paths
        self.base_datadir = Path("tests/data/fixtures/")
        self.numbering_data = self.base_datadir / "anarci/curated_alignments"
```

### Test Fixture: GermlineManager
Location: `src/sadie/germlines/tests/test_g3_regression.py`

```python
@pytest.fixture
def germline_manager():
    from sadie.germlines.manager import GermlineManager
    return GermlineManager(providers=["imgt"])
```

---

## 9. Test Data Files

### CDR3 Known Bugs Test Data
Location: `tests/data/fixtures/cdr3_known_bugs.fasta`
Purpose: Test sequences with known CDR3 edge cases

### Regression Test Sequences
Location: `src/sadie/germlines/tests/data/regression/FR025b_sequences.fasta`
Purpose: Verify specific gene alleles are returned

---

## 10. Testing Patterns

### Pattern: Expected Failures
Document known test failures in assertions:
```python
expected_null_cdr3 = ["Seq5"]  # Has truncated J region
unexpected_null_cdr3 = [seq for seq in sequences_with_null_cdr3 
                       if seq not in expected_null_cdr3]
```

### Pattern: Feature Flag Testing
Use monkeypatch for feature flags:
```python
def test_with_germlines_enabled(self, monkeypatch):
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
    # Test with new backend
    
def test_with_germlines_disabled(self, monkeypatch):
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")
    # Test with legacy backend
```

### Pattern: Path Mocking
Mock filesystem paths for unit tests:
```python
def test_missing_database_error(self, monkeypatch, tmp_path):
    empty_igblast = tmp_path / "igblast" / "database" / "human"
    empty_igblast.mkdir(parents=True)
    
    monkeypatch.setattr("sadie.germlines.get_germlines_base_dir", 
                        lambda: tmp_path)
```

### Pattern: Alignment Tolerance
Compare alignments with positional tolerance:
```python
# Only check last positions for conservation
assert numbering_align[-5:] == g3_align[-5:]
assert numbering_align[-35:].count("-") == g3_align[-35:].count("-")
```

---

## 11. Running Tests

### Run All CDR3 Tests
```bash
pytest tests/unit/airr/test_cdr3.py -v
```

### Run G3 Regression Tests
```bash
pytest src/sadie/germlines/tests/test_g3_regression.py -v
```

### Run AIRR Integration Tests
```bash
pytest tests/unit/germlines/test_airr_integration.py -v
```

### Run All Segment-Related Tests
```bash
pytest -v -k "cdr3 or j_gene or aux or segment"
```

---

## 12. Test Coverage Requirements

### Required Coverage Areas
1. **Aux file format**: 5-column validation
2. **CDR3 positions**: Start/end validity, extraction accuracy
3. **J gene reference**: All known alleles have correct data
4. **Gene call presence**: V, D, J calls populated
5. **Position consistency**: CDR3 length matches positions
6. **Offline operation**: All data available locally
7. **Backwards compatibility**: Legacy paths still work
