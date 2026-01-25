# SADIE Testing Strategy

## Overview

This document describes testing strategy for SADIE, with emphasis on the Reference module and Germlines integration (v1.2).

---

## 1. Test Organization

### Directory Structure
```
tests/
├── conftest.py                           # Global fixtures
├── __init__.py
├── data/
│   └── fixtures/                         # Test data files
│       ├── fasta_inputs/
│       ├── airr_tables/
│       └── reference/
│           ├── short_reference.yml
│           ├── duplicated_in_source.yml
│           └── duplicated_inter_source.yml
├── unit/
│   ├── __init__.py
│   ├── reference/
│   │   ├── test_reference.py             # Reference class tests
│   │   └── test_advanced_reference.py
│   ├── germlines/
│   │   ├── __init__.py
│   │   ├── conftest.py                   # Germlines-specific fixtures
│   │   ├── test_reference_integration.py # Reference + Germlines integration
│   │   ├── test_compliance.py            # AIRR compliance tests
│   │   ├── test_manager.py               # GermlineManager tests
│   │   ├── test_integration.py
│   │   ├── test_airr_integration.py
│   │   └── test_offline_operation.py
│   └── airr/
│       ├── test_airr.py
│       ├── test_methods.py
│       ├── test_cdr3.py
│       └── test_igblast.py
└── integration/
    ├── airr/
    │   └── test_airr_intergration.py
    └── reference/
        └── test_reference_integration.py
```

---

## 2. Pytest Configuration

### Location: `pytest.ini`
```ini
[pytest]
addopts = -x
testpaths = tests
norecursedirs = antibody_objects_old .git __pycache__ .eggs *.egg-info dist build docs
log_format = %(asctime)s.%(msecs)03d %(levelname)s: %(message)s
log_date_format = %Y-%m-%d %H:%M:%S
log_cli_format = %(asctime)s.%(msecs)03d %(levelname)s: %(message)s
log_cli_date_format = %Y-%m-%d %H:%M:%S
log_cli=True
```

### Configuration Details
| Option | Value | Purpose |
|--------|-------|---------|
| `addopts` | `-x` | Stop on first failure |
| `testpaths` | `tests` | Test discovery location |
| `log_cli` | `True` | Enable logging during tests |

### Coverage Configuration (pyproject.toml)
```toml
[tool.coverage.run]
source = ["sadie"]
omit = ["*/tests/*", "*/test_*.py"]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "def __repr__",
    "if TYPE_CHECKING:",
    "@abstractmethod",
]
```

---

## 3. Test Fixtures

### Global Fixture: SadieFixture
Location: `tests/conftest.py`

```python
class SadieFixture(AirrSequences, AirrTables, ReferenceFixtures, NumberingFixtures):
    """Master fixture combining all test data categories."""

    def __init__(self, tmp_path_factory):
        tmp_path = tmp_path_factory.mktemp("sadie_fixture")
        base_datadir = Path("tests/data/fixtures/")

        AirrSequences.__init__(self, tmp_path, base_datadir)
        AirrTables.__init__(self, tmp_path, base_datadir)
        ReferenceFixtures.__init__(self, tmp_path, base_datadir)
        NumberingFixtures.__init__(self, tmp_path, base_datadir)

@pytest.fixture(scope="session", autouse=True)
def fixture_setup(tmp_path_factory: pytest.TempPathFactory):
    return SadieFixture(tmp_path_factory)
```

### Reference Fixtures
Location: `tests/conftest.py` (ReferenceFixtures class)

```python
class ReferenceFixtures:
    """Reference-related test fixtures."""

    def get_reference_dataset_csv(self) -> Path:
        """Reference dataset as CSV."""
        return self.reference_data / "reference_object_dataframe.csv.gz"

    def get_duplicated_yaml(self) -> Path:
        """YAML file with duplicated entries (for error testing)."""
        return self.reference_data / "duplicated_in_source.yml"

    def get_shortened_yaml(self) -> Path:
        """Short YAML file for quick tests."""
        return self.reference_data / "short_reference.yml"
```

### Germlines Fixtures
Location: `tests/unit/germlines/conftest.py`

```python
@pytest.fixture
def network_disabled(monkeypatch):
    """Disable network for offline operation tests."""
    def block_socket(*args, **kwargs):
        raise OSError("Network is disabled for this test")
    monkeypatch.setattr(socket, "socket", block_socket)
    yield

@pytest.fixture
def enable_germlines(monkeypatch):
    """Enable germlines module via feature flag."""
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
    yield

@pytest.fixture
def disable_germlines(monkeypatch):
    """Disable germlines module (use G3 API)."""
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")
    yield
```

---

## 4. Reference Module Tests

### Test File: `tests/unit/reference/test_reference.py`

#### Test: YAML Loading
```python
def test_yaml(tmp_path_factory, fixture_setup) -> None:
    """Test YAML reference loading."""
    yaml_object = YamlRef()
    assert yaml_object.get_names() == {"clk", "dog", "human", "mouse", "rabbit", "se09", "rat", "macaque"}
    assert len(yaml_object.get_genes("human", "imgt", "human")) == 479
```

#### Test: Pydantic Model Validation
```python
def test_check_models(fixture_setup) -> None:
    """Test Pydantic model validation."""
    # Valid entry
    GeneEntry(species="human", gene="IGHV3-10*01", source="custom")

    # Invalid gene position
    with pytest.raises(ValidationError):
        GeneEntry(species="human", gene="IGHZ3-10*01", source="custom")

    # Invalid source
    with pytest.raises(ValidationError):
        GeneEntry(species="human", gene="IGHV1-69*01", source="invalid_source")
```

#### Test: All Four Sources Valid (v1.2)
```python
def test_source_validation_all_providers() -> None:
    """Test all four germline sources are valid."""
    for source in ["imgt", "ogrdb", "vdjbase", "custom"]:
        entry = GeneEntry(species="human", gene="IGHV1-69*01", source=source)
        assert entry.source == source
```

#### Test: Reference with Germlines Backend
```python
def test_reference_use_germlines() -> None:
    """Test Reference class with use_germlines=True."""
    ref = Reference(use_germlines=True)
    assert ref.use_germlines is True

    ref.add_gene({"species": "human", "gene": "IGHV1-69*01", "source": "imgt"})
    assert len(ref.data) == 1

    # Verify deterministic _id field
    df = ref.get_dataframe()
    assert "_id" in df.columns
    assert len(df["_id"].iloc[0]) == 24  # SHA-256 truncated
```

---

## 5. Reference Integration Tests

### Test File: `tests/unit/germlines/test_reference_integration.py`

#### Test Class: TestReferenceIntegration
```python
class TestReferenceIntegration:
    """Test Reference system with germlines module backend."""

    def test_reference_with_germlines_backend(self, monkeypatch):
        """T027: Reference initializes with use_germlines=True."""
        ref = Reference(use_germlines=True)
        assert ref.use_germlines is True

        gene_dict = {"species": "human", "gene": "IGHV1-69*01", "source": "imgt"}
        ref.add_gene(gene_dict)

        df = ref.get_dataframe()
        expected_fields = ["source", "common", "gene", "sequence"]
        for field in expected_fields:
            assert field in df.columns

    def test_output_format_consistency(self, monkeypatch):
        """T029: Output format consistency between backends."""
        ref = Reference(use_germlines=True)
        ref.add_gene({"species": "human", "gene": "IGHV1-69*01", "source": "imgt"})

        df = ref.get_dataframe()

        # G3-compatible fields
        required_fields = ["source", "common", "gene", "label", "gene_segment", "receptor", "sequence"]
        for field in required_fields:
            assert field in df.columns

        # IMGT nested fields
        imgt_fields = ["imgt.sequence_gapped", "imgt.sequence_gapped_aa", "imgt.imgt_functional"]
        for field in imgt_fields:
            assert field in df.columns
```

#### Test Class: TestGermlineToG3Adapter
```python
class TestGermlineToG3Adapter:
    """Test G3 adapter functionality."""

    def test_adapter_transform(self, monkeypatch):
        """Test adapter transforms GermlineGene to G3 format."""
        from sadie.germlines import get_gene_by_name
        from sadie.germlines.g3_adapter import GermlineToG3Adapter

        gene = get_gene_by_name("IGHV1-69*01", "human")
        adapter = GermlineToG3Adapter()
        g3_dict = adapter.to_g3_format(gene)

        assert "gene" in g3_dict
        assert "source" in g3_dict
        assert "imgt" in g3_dict
        assert "sequence_gapped" in g3_dict["imgt"]
```

---

## 6. Compliance Tests

### Test File: `tests/unit/germlines/test_compliance.py`

#### Test: Priority Order (FR-004)
```python
class TestPriorityOrder:
    """Test default provider priority order."""

    def test_default_priority_order(self):
        """Verify: custom > ogrdb > vdjbase > imgt."""
        expected = ["custom", "ogrdb", "vdjbase", "imgt"]
        assert GermlineManager.DEFAULT_PROVIDERS == expected
```

#### Test: No G3 Fallback (NFR-002)
```python
def test_no_g3_fallback_when_germlines_enabled(self, monkeypatch):
    """Verify no silent fallback to G3 when germlines selected."""
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
    clear_feature_flag_cache()

    with pytest.raises(ValueError) as exc_info:
        GermlineData("nonexistent_species_xyz")

    assert "not found" in str(exc_info.value).lower()
```

#### Test: Species/Chain/Segment Parity (FR-010)
```python
class TestParity:
    """Test species/chain/segment parity."""

    def test_human_coverage(self):
        """Human has V, D, J, C for H, K, L chains."""
        manager = GermlineManager()
        for chain in ["H", "K", "L"]:
            for segment in ["V", "D", "J", "C"]:
                if segment == "D" and chain in ["K", "L"]:
                    continue  # No D genes for light chains
                genes = manager.get_genes("human", segment, chain)
                assert len(genes) > 0
```

---

## 7. Test Patterns

### Pattern: Feature Flag Testing
```python
def test_with_germlines_enabled(self, monkeypatch):
    """Test with germlines module."""
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
    clear_feature_flag_cache()
    # Test new backend

def test_with_germlines_disabled(self, monkeypatch):
    """Test with legacy G3 API."""
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")
    clear_feature_flag_cache()
    # Test legacy backend
```

### Pattern: Skip Tests Based on Availability
```python
def _macaque_available() -> bool:
    """Check if macaque germlines are available."""
    macaque_path = get_germlines_base_dir() / "igblast" / "Ig" / "internal_data" / "macaque"
    return macaque_path.exists()

skip_no_macaque = pytest.mark.skipif(
    not _macaque_available(),
    reason="macaque germlines not available"
)

@skip_no_macaque
def test_macaque_annotation():
    ...
```

### Pattern: Validation Error Testing
```python
def test_invalid_source_rejected():
    """Test invalid source raises ValidationError."""
    with pytest.raises(ValidationError) as exc_info:
        GeneEntry(species="human", gene="IGHV1-69*01", source="invalid")

    assert "not a valid source" in str(exc_info.value)
    assert "choices" in str(exc_info.value)  # Helpful error message
```

### Pattern: Expected Failures
```python
def test_cdr3_field_not_none(self, fasta_path) -> None:
    """Test CDR3 population with known exceptions."""
    result = airr_api.run_fasta(str(fasta_path))
    cdr3_nulls = result["cdr3"].isna()

    expected_null_cdr3 = ["Seq5"]  # Known exception: truncated J
    unexpected = [seq for seq in result[cdr3_nulls]["sequence_id"]
                  if seq not in expected_null_cdr3]

    assert not unexpected, f"Unexpected null CDR3: {unexpected}"
```

### Pattern: Parametrized Tests
```python
@pytest.mark.parametrize(
    "codon, aa, expected",
    [
        ("", "M", "ATG"),
        ("ATG", "M", "ATG"),
        ("FAKE", "M", "ATG"),
    ]
)
def test_find_best_codon(codon: str, aa: str, expected: str) -> None:
    assert find_best_codon(codon, aa) == expected
```

### Pattern: Test Classes for Organization
```python
class TestReferenceIntegration:
    """Test Reference + Germlines integration."""

    def test_reference_with_germlines_backend(self):
        """T027: Test Reference with germlines."""
        ...

    def test_output_format_consistency(self):
        """T029: Test output format."""
        ...
```

---

## 8. Test Data Management

### Fixture Data Location
```
tests/data/fixtures/
├── fasta_inputs/           # FASTA test files
│   ├── PG9_H.fasta
│   ├── catnap_nt_heavy.fasta
│   └── OAS_subsample_1000.fasta
├── airr_tables/            # Pre-computed AIRR results
│   ├── catnap_heavy_airrtable.feather
│   └── dog_igh.tsv.gz
└── reference/              # Reference YAML files
    ├── duplicated_in_source.yml     # Error case: duplicates
    ├── duplicated_inter_source.yml  # Error case: cross-source dups
    └── short_reference.yml          # Quick test file
```

### Test Data Access Pattern
```python
def test_yaml_duplicates(fixture_setup: SadieFixture):
    """Test duplicate detection in YAML."""
    with pytest.raises(ValueError):
        YamlRef(fixture_setup.get_duplicated_yaml())
```

---

## 9. Running Tests

### Run All Tests
```bash
pytest tests/ -v
```

### Run Unit Tests Only
```bash
pytest tests/unit/ -v
```

### Run Reference Module Tests
```bash
pytest tests/unit/reference/ tests/unit/germlines/test_reference_integration.py -v
```

### Run Compliance Tests
```bash
pytest tests/unit/germlines/test_compliance.py -v
```

### Run with Coverage
```bash
pytest tests/ --cov=sadie --cov-report=html
```

### Run Specific Test Class
```bash
pytest tests/unit/germlines/test_reference_integration.py::TestReferenceIntegration -v
```

### Run Tests Matching Pattern
```bash
pytest -v -k "reference or germline"
```

---

## 10. Coverage Requirements

### Minimum Coverage
Target: **80%** coverage for core modules

### Required Coverage Areas

| Area | Description |
|------|-------------|
| Reference models | GeneEntry, GeneEntries validation |
| Source validation | All 4 sources (imgt, ogrdb, vdjbase, custom) |
| YAML loading | YamlRef parsing and duplicate detection |
| Germlines integration | Reference(use_germlines=True) |
| G3 adapter | GermlineGene to G3 format conversion |
| Error handling | ValidationError messages |
| Feature flags | Backend switching |

### Coverage Exclusions
```toml
[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "def __repr__",
    "if TYPE_CHECKING:",
    "@abstractmethod",
]
```

---

## 11. Test Naming Conventions

### Test Function Names
```python
# Pattern: test_<what>_<condition>
def test_yaml_loads_correctly():
    ...

def test_invalid_source_raises_error():
    ...

def test_reference_with_germlines_backend():
    ...
```

### Test IDs (for Traceability)
```python
def test_reference_with_germlines_backend(self):
    """T027: Test Reference with germlines."""

def test_output_format_consistency(self):
    """T029: Verify output format between backends."""
```

---

## 12. Integration Test Patterns

### Network-Dependent Tests
```python
def test_g3_api_connection():
    """Test G3 API connectivity (requires network)."""
    ref = Reference(use_germlines=False)
    # Will connect to G3 API on endpoint validation
```

### Database-Dependent Tests
```python
def test_human_genes_available():
    """Test human germline database exists."""
    manager = GermlineManager()
    genes = manager.get_genes("human", "V", "H")
    assert len(genes) > 100  # Expect many V genes
```

---

## 13. Debugging Test Failures

### Enable Verbose Logging
```bash
pytest tests/unit/germlines/test_reference_integration.py -v --log-cli-level=DEBUG
```

### Run Single Test
```bash
pytest tests/unit/germlines/test_reference_integration.py::TestReferenceIntegration::test_reference_with_germlines_backend -v
```

### Show Captured Output
```bash
pytest tests/ -v -s  # -s shows print statements
```

### Drop into Debugger on Failure
```bash
pytest tests/ --pdb
```
