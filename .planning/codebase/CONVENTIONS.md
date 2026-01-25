# SADIE Coding Conventions

## Overview

This document describes coding conventions for the SADIE codebase, with emphasis on the Reference and Germlines modules (v1.2 unification).

---

## 1. Python Style

### Type Hints
Always use type hints for function signatures and class attributes:

```python
# Function signatures
def get_gene_by_name(name: str, species: str) -> Optional[GermlineGene]:
    ...

def build_for_species(self, species: str, source_dir: Path, output_file: Path) -> bool:
    ...

# Class attributes with Field() for Pydantic
class GermlineGene(BaseModel):
    name: str = Field(..., description="Gene name (e.g., IGHV1-69*01)")
    species: str = Field(..., description="Species (e.g., human)")
    sequence: str = Field(..., description="Ungapped nucleotide sequence")
    sequence_gapped: Optional[str] = Field(None, description="IMGT-gapped nucleotide")
```

### Import Conventions
```python
# Standard library first
from __future__ import annotations
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

# Third-party
import pandas as pd
import pytest
from Bio.Seq import Seq
from pydantic import BaseModel, Field, field_validator

# Local imports
from sadie.reference.models import GeneEntry, GeneEntries
from sadie.germlines.models import GermlineGene
```

### Line Length
- Maximum line length: **120 characters** (configured in pyproject.toml)
- Target versions: Python 3.10, 3.11, 3.12, 3.13

---

## 2. Naming Conventions

### Functions and Variables
- Use `snake_case` for functions and variables
- Use descriptive names that indicate purpose

```python
# Good
def get_j_gene_data(allele_name: str, chain: str) -> tuple:
    ...

def validate_aux_file(self, aux_file: Path) -> bool:
    ...

# Variables
species_available: List[str] = []
gene_names = [g.name for g in genes]
```

### Classes
- Use `PascalCase` for class names
- Suffix exception classes with `Error` or keep domain-specific names

```python
class GermlineGene(BaseModel):
    ...

class GermlineManager:
    ...

class G3Error(Exception):
    ...

class BadDataSet(Error):
    ...
```

### Constants
- Use `UPPER_SNAKE_CASE` for module-level constants

```python
# Valid sources
VALID_SOURCES = ["imgt", "ogrdb", "vdjbase", "custom"]

# Reference data
HUMAN_J_GENE_DATA = {
    "IGHJ1*01": (0, 17, 1),
    "IGHJ2*01": (1, 18, 1),
}

# Chain type mapping
CHAIN_TYPE_MAP = {
    "H": "JH",
    "K": "JK",
    "L": "JL",
}
```

---

## 3. Docstrings

### NumPy Style (Preferred)
Use NumPy-style docstrings with Parameters, Returns, and Examples sections:

```python
def get_j_gene_data(allele_name: str, chain: str) -> tuple:
    """
    Get J gene reference data for an allele.

    Parameters
    ----------
    allele_name : str
        Full allele name (e.g., "IGHJ1*01")
    chain : str
        Chain type (H, K, or L)

    Returns
    -------
    tuple
        (reading_frame, chain_type, cdr3_end, is_functional)

    Examples
    --------
    >>> data = get_j_gene_data("IGHJ4*01", "H")
    >>> data[0]  # reading_frame
    2
    """
```

### Class Docstrings
Include class-level description with Attributes section:

```python
class GermlineGene(BaseModel):
    """
    Unified germline gene model across all providers.

    This model represents a single germline gene sequence with all
    associated metadata. It normalizes data from IMGT, OGRDB, and
    custom sources into a common format.

    Attributes
    ----------
    name : str
        Gene name in IMGT format (e.g., "IGHV1-69*01")
    species : str
        Species name (e.g., "human", "mouse")
    segment : str
        Segment type: "V", "D", "J", or "C"
    """
```

---

## 4. Error Handling Patterns

### Exception Hierarchy
Location: `src/sadie/airr/exceptions.py`

```python
class Error(Exception):
    """Base class for exceptions in this module."""

class BadDataSet(Error):
    """Exception raised for invalid species/dataset."""
    def __init__(self, requested_type: str, accepted_types: List[str]):
        self.requested_type = requested_type
        self.accepted_types = accepted_types

    def __str__(self) -> str:
        return f"{self.requested_type} dataset, avail datasets{sorted(self.accepted_types)}"
```

### Reference Module Exceptions
Location: `src/sadie/reference/reference.py`

```python
class G3Error(Exception):
    """Exception for G3 API errors."""
    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
```

### Validation Pattern
```python
# Raise specific exceptions with helpful messages
if species not in available_species:
    raise ValueError(
        f"Species '{species}' not found. "
        f"Available species: {sorted(available_species)}"
    )

# Use pytest.raises for validation error testing
with pytest.raises(ValidationError) as exc_info:
    GeneEntry(species="human", gene="IGHV1-69*01", source="invalid_source")
assert "not a valid source" in str(exc_info.value)
```

---

## 5. Pydantic Model Conventions

### Model Definition (Pydantic v2)
Location: `src/sadie/reference/models.py`, `src/sadie/germlines/models.py`

```python
from pydantic import BaseModel, Field, field_validator

class GeneEntry(BaseModel):
    """V, D, or J Gene Entry with validation."""
    species: str
    gene: str
    source: str

    @field_validator("source")
    @classmethod
    def check_source(cls, v: str) -> str:
        if v not in VALID_SOURCES:
            raise ValueError(f"{v} is not a valid source, choices are {VALID_SOURCES}")
        return v

    @field_validator("gene")
    @classmethod
    def check_vgene(cls, v: str) -> str:
        if v[3] not in ["V", "D", "J"]:
            raise ValueError(f"gene must contain V, D or J at 3rd index, got {v[3]} in {v}")
        return v
```

### Field Validators
- Use `@field_validator` decorator (Pydantic v2 syntax)
- Include `@classmethod` decorator
- Normalize values in validators (e.g., uppercase)

```python
@field_validator("segment")
@classmethod
def validate_segment(cls, v: str) -> str:
    """Validate segment is V, D, J, or C."""
    v = v.upper()
    if v not in ["V", "D", "J", "C"]:
        raise ValueError(f"Segment must be V, D, J, or C, got: {v}")
    return v

@field_validator("sequence")
@classmethod
def validate_sequence(cls, v: str) -> str:
    """Validate sequence contains only valid nucleotides."""
    v = v.upper()
    valid_chars = set("ACGTN")
    invalid = set(v) - valid_chars
    if invalid:
        raise ValueError(f"Sequence contains invalid characters: {invalid}")
    return v
```

### Model String Representations
```python
def __str__(self) -> str:
    """String representation."""
    return f"{self.name} ({self.source})"

def __repr__(self) -> str:
    """Detailed representation."""
    return (
        f"GermlineGene(name='{self.name}', "
        f"species='{self.species}', "
        f"segment='{self.segment}')"
    )
```

---

## 6. Reference Module Conventions

### YAML Schema
Location: `src/sadie/reference/data/reference.yml`

```yaml
# Structure: name -> source -> species -> [genes]
human:
  imgt:
    human:
      - IGHV1-2*01
      - IGHV1-2*02
      - IGHV1-69*01
      - IGHD3-10*01
      - IGHJ4*02
  ogrdb:
    human:
      - IGHV1-69*i01  # Novel allele
```

### Valid Sources (v1.2)
```python
# src/sadie/reference/models.py
VALID_SOURCES = ["imgt", "ogrdb", "vdjbase", "custom"]
```

### Reference Class Pattern
```python
class Reference:
    """Reference class to handle reference databases."""

    _endpoint = "https://g3.jordanrwillis.com/api/v1/genes"

    def __init__(self, endpoint: str = _endpoint, use_germlines: bool = False):
        self.data: List[Dict[str, str]] = []
        self.use_germlines = use_germlines
        
        if not use_germlines:
            self.endpoint = endpoint
```

### Gene Entry Validation
```python
def add_gene(self, gene: Dict[str, str]) -> None:
    """Add a single gene to the reference data."""
    gene_valid = GeneEntry(**gene)  # Validates via Pydantic
    self.data.append(self._get_gene(gene_valid))
```

---

## 7. Germlines Module Conventions

### Provider Priority Order
```python
# Default: custom > ogrdb > vdjbase > imgt (novel alleles prioritized)
DEFAULT_PROVIDERS = ["custom", "ogrdb", "vdjbase", "imgt"]
```

### GermlineManager Pattern
```python
class GermlineManager:
    """Priority-based germline database manager."""

    DEFAULT_PROVIDERS = ["custom", "ogrdb", "vdjbase", "imgt"]

    def __init__(self, providers: List[str] = None):
        self.provider_names = providers or self.DEFAULT_PROVIDERS
```

### Feature Flag Pattern
```python
# Environment variable controls backend
# SADIE_USE_GERMLINES_MODULE=true (default) -> germlines module
# SADIE_USE_GERMLINES_MODULE=false -> G3 API

def use_germlines_module() -> bool:
    """Check if germlines module should be used."""
    return os.environ.get("SADIE_USE_GERMLINES_MODULE", "true").lower() in ("true", "1", "yes")
```

### G3 Adapter Pattern
Location: `src/sadie/germlines/g3_adapter.py`

```python
class GermlineToG3Adapter:
    """Transform GermlineGene objects to G3 API format."""

    def to_g3_format(self, gene: GermlineGene) -> Dict[str, Any]:
        """Convert GermlineGene to G3 API response format."""
        return {
            "_id": self._generate_id(gene.source, gene.species, gene.name),
            "source": gene.source,
            "common": gene.species,
            "gene": gene.name,
            "sequence": gene.sequence,
            "imgt": {
                "sequence_gapped": gene.sequence_gapped or "",
                "imgt_functional": gene.functionality,
            }
        }
```

---

## 8. Logging Conventions

### Logger Setup
```python
import logging

# Module-level logger
logger = logging.getLogger(__name__)

# Or named logger
logger = logging.getLogger("Reference")
```

### Log Levels
```python
# Debug: Detailed diagnostic information
logger.debug(f"Processing gene: {gene_name}")

# Info: General operational messages
logger.info(f"Building auxiliary file for {species}")

# Warning: Something unexpected but handled
logger.warning(f"No gapped sequence for {gene_name}, using ungapped")

# Error: Something failed
logger.error(f"Failed to parse {fasta_path}: {e}")
```

---

## 9. Aux File Format Convention

### 5-Column Tab-Separated Format
```
<gene_name>\t<reading_frame>\t<chain_type>\t<cdr3_end>\t<is_functional>
```

| Column | Name | Values |
|--------|------|--------|
| 1 | gene_name | `IGHJ1*01`, `IGKJ2*01` |
| 2 | reading_frame | `0`, `1`, `2` |
| 3 | chain_type | `JH`, `JK`, `JL` |
| 4 | cdr3_end | Integer (1-based position) |
| 5 | is_functional | `1` (functional), `0` (non-functional) |

### Example
```
IGHJ1*01	0	JH	17	1
IGHJ2*01	1	JH	18	1
IGKJ1*01	1	JK	6	1
```

---

## 10. File Organization

### Source Structure
```
src/sadie/
├── reference/
│   ├── __init__.py
│   ├── models.py          # Pydantic models (GeneEntry, GeneEntries)
│   ├── reference.py       # Reference class
│   ├── yaml.py            # YAML reference loader
│   └── data/
│       └── reference.yml  # Reference database YAML
├── germlines/
│   ├── __init__.py
│   ├── models.py          # GermlineGene model
│   ├── manager.py         # GermlineManager (priority-based)
│   ├── g3_adapter.py      # G3 format adapter
│   ├── providers/         # Data providers
│   │   ├── base.py
│   │   ├── imgt.py
│   │   ├── ogrdb.py
│   │   ├── vdjbase.py
│   │   └── custom.py
│   └── builders/          # File builders
│       ├── aux.py
│       ├── blast.py
│       └── j_gene_data.py
└── airr/
    ├── exceptions.py      # Exception hierarchy
    └── igblast/
        └── germline.py    # GermlineData class
```

### Test Structure
```
tests/
├── conftest.py                    # Shared fixtures (SadieFixture)
├── unit/
│   ├── reference/
│   │   ├── test_reference.py      # Reference class tests
│   │   └── test_advanced_reference.py
│   ├── germlines/
│   │   ├── conftest.py            # Germlines fixtures
│   │   ├── test_reference_integration.py  # Reference + Germlines
│   │   ├── test_compliance.py     # AIRR compliance tests
│   │   └── test_manager.py        # GermlineManager tests
│   └── airr/
│       └── test_airr.py
└── integration/
    └── reference/
        └── test_reference_integration.py
```
