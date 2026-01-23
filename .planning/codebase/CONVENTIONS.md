# G3 Coding Conventions for Segment Handling

## Overview

This document describes the coding conventions used in SADIE for G3 backend segment discovery, particularly for J gene segments and CDR3/junction boundary calculations.

---

## 1. Aux File Format Convention

### File Structure
IgBLAST auxiliary files follow a strict **5-column, tab-separated format**:

```
<gene_name>\t<reading_frame>\t<chain_type>\t<cdr3_end>\t<is_functional>
```

### Column Definitions

| Column | Name | Type | Description | Values |
|--------|------|------|-------------|--------|
| 1 | `gene_name` | string | J gene allele name | `IGHJ1*01`, `IGKJ2*01`, etc. |
| 2 | `reading_frame` | int | Codon frame offset (0-2) | `0`, `1`, `2` |
| 3 | `chain_type` | string | Chain marker | `JH`, `JK`, `JL` |
| 4 | `cdr3_end` | int | CDR3 end position (1-based, nucleotide) | Chain-specific values |
| 5 | `is_functional` | int | Functional status | `1` (functional), `0` (non-functional) |

### Example Aux File Content
```
IGHJ1*01	0	JH	17	1
IGHJ2*01	1	JH	18	1
IGHJ3*01	1	JH	15	1
IGHJ4*01	2	JH	13	1
IGKJ1*01	1	JK	6	1
IGLJ1*01	1	JL	6	1
```

### Location Convention
- Production aux files: `src/sadie/germlines/igblast/aux_db/{species}_gl.aux`
- Legacy aux files: `src/sadie/airr/data/germlines/aux_db/imgt/{species}_gl.aux`

---

## 2. J Gene Reference Data Convention

### Data Structure
J gene reference data is stored in `src/sadie/germlines/builders/j_gene_data.py` as a dictionary:

```python
# Format: {allele: (reading_frame, cdr3_end, is_functional)}
HUMAN_J_GENE_DATA = {
    "IGHJ1*01": (0, 17, 1),
    "IGHJ2*01": (1, 18, 1),
    "IGHJ3*01": (1, 15, 1),
    # ...
}
```

### Chain Type Mapping
```python
CHAIN_TYPE_MAP = {
    "H": "JH",  # Heavy chain J
    "K": "JK",  # Kappa chain J
    "L": "JL",  # Lambda chain J
}
```

### Default Fallback Values
When a J gene allele is not found in reference data, use chain-specific defaults:

| Chain | Reading Frame | Chain Type | CDR3 End | Is Functional |
|-------|--------------|------------|----------|---------------|
| H | 1 | JH | 15 | 1 |
| K | 1 | JK | 6 | 1 |
| L | 1 | JL | 6 | 1 |

---

## 3. Segment Position Validation Conventions

### CDR3 Position Fields (AIRR Format)
- `cdr3_start`: 1-based inclusive start position
- `cdr3_end`: 1-based inclusive end position
- Positions are relative to the query sequence

### Validation Rules
1. **Position integrity**: `cdr3_end > cdr3_start` must always hold
2. **Length match**: `len(cdr3) == cdr3_end - cdr3_start + 1`
3. **Extraction validation**: Extracted CDR3 from sequence must match reported CDR3

### Position Extraction Pattern
```python
# AIRR positions are 1-based, convert to 0-based for Python
start = int(row["cdr3_start"]) - 1  # Convert to 0-based
end = int(row["cdr3_end"])          # End is inclusive in AIRR
extracted_cdr3 = row["sequence"][start:end]
```

---

## 4. Aux File Builder Conventions

### Builder Location
`src/sadie/germlines/builders/aux.py`

### Key Methods

| Method | Purpose |
|--------|---------|
| `build_for_species()` | Build complete aux file for a species |
| `_process_segment()` | Process J segment FASTA to aux entries |
| `_create_aux_entry()` | Create single aux line for J gene |
| `validate_aux_file()` | Validate aux file format compliance |

### Building Process
1. Iterate chains: H, K, L (only J segments)
2. Parse gapped FASTA from `normalized/{species}/gapped/IG{chain}J.fasta`
3. Look up reference data via `get_j_gene_data(gene_name, chain)`
4. Write 5-column tab-separated output

### Validation Rules
```python
def validate_aux_file(self, aux_file: Path) -> bool:
    lines = aux_file.read_text().strip().split("\n")
    for line in lines:
        fields = line.split("\t")
        if len(fields) != 5:  # Must have exactly 5 columns
            return False
    return True
```

---

## 5. Germline Data Integration Conventions

### Feature Flag
The germlines module is controlled by `SADIE_USE_GERMLINES_MODULE` environment variable:
- Default: `true` (uses germlines module)
- Legacy: `false` (uses G3 API paths)

### Path Resolution
```python
# Germlines module path (new)
base_dir = germlines_module_path / "igblast"

# Legacy path
base_dir = sadie.airr.data.germlines / "Ig"
```

### GermlineData Class Pattern
```python
class GermlineData:
    def __init__(self, species: str, database_dir: Path = None):
        self.species = species
        self.base_dir = self._resolve_base_dir(database_dir)
        self.v_gene_dir = self.base_dir / f"blastdb/{species}/{species}_V"
        self.j_gene_dir = self.base_dir / f"blastdb/{species}/{species}_J"
        self.aux_path = self.base_dir / f"aux_db/imgt/{species}_gl.aux"
```

---

## 6. Coding Style Conventions

### Naming Conventions
- Functions: `snake_case` (e.g., `get_j_gene_data`, `validate_aux_file`)
- Classes: `PascalCase` (e.g., `AuxFileBuilder`, `GermlineData`)
- Constants: `UPPER_SNAKE_CASE` (e.g., `HUMAN_J_GENE_DATA`, `CHAIN_TYPE_MAP`)

### Type Hints
Always use type hints for function signatures:
```python
def get_j_gene_data(allele_name: str, chain: str) -> tuple:
    ...

def validate_aux_file(self, aux_file: Path) -> bool:
    ...
```

### Docstrings
Use NumPy-style docstrings with Parameters and Returns sections:
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
    """
```

### Logging
Use Python's logging module for diagnostic output:
```python
import logging
logger = logging.getLogger(__name__)

logger.info(f"Building auxiliary file for {species}")
logger.debug(f"No file: {fasta_path}")
logger.error(f"Failed to parse {fasta_path}: {e}")
```

---

## 7. Error Handling Conventions

### Exception Hierarchy
AIRR/IgBLAST exceptions are in `sadie/airr/exceptions.py`:
- `BadIgBLASTArgument`: Invalid argument to IgBLAST
- `BadIgBLASTExe`: IgBLAST executable not found
- `BadIgDATA`: Invalid IGDATA path
- `IgBLASTRunTimeError`: Runtime error during execution
- `MissingIgBLASTArgument`: Required argument missing

### Validation Pattern
```python
if not aux_file.exists():
    logger.error(f"Aux file doesn't exist: {aux_file}")
    return False

for line in lines:
    fields = line.split("\t")
    if len(fields) != 5:
        logger.error(f"Invalid aux line (expected 5 columns): {line}")
        return False
```

---

## 8. File Organization

### Source Structure
```
src/sadie/germlines/
├── builders/
│   ├── aux.py           # Aux file builder
│   ├── j_gene_data.py   # J gene reference data
│   └── gapper.py        # Sequence gapping
├── igblast/
│   └── aux_db/          # Production aux files
└── tests/
    └── test_*.py        # Module-level tests
```

### Test Structure
```
tests/
├── unit/
│   ├── airr/
│   │   ├── test_cdr3.py      # CDR3 field tests
│   │   └── test_igblast.py   # IgBLAST integration
│   ├── germlines/
│   │   └── test_airr_integration.py
│   └── renumbering/
│       └── test_g3.py        # G3 alignment tests
└── integration/
    └── airr/
        └── test_airr_integration.py
```
