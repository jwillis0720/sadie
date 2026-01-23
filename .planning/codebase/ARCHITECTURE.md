# G3 Backend Architecture for AIRR Segment Discovery

## Overview

The SADIE G3 backend provides immunoglobulin repertoire annotation by wrapping IgBLAST and processing its output into AIRR-compliant format. This document focuses on how V, D, J, and C region segment positions are discovered and CDR3/junction boundaries are determined.

## Core Architecture Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        AIRR Annotation Pipeline                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Input Sequence                                                              │
│       │                                                                      │
│       ▼                                                                      │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │  Airr Class (src/sadie/airr/airr.py)                                    │ │
│  │  - Orchestrates annotation pipeline                                      │ │
│  │  - Configures IgBLAST parameters                                         │ │
│  │  - Handles adaptive penalty adjustment                                   │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│       │                                                                      │
│       ▼                                                                      │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │  GermlineData Class (src/sadie/airr/igblast/germline.py)                │ │
│  │  - Resolves paths for V/D/J/C BLAST databases                           │ │
│  │  - Locates auxiliary file for J gene CDR3 boundaries                    │ │
│  │  - Sets IGDATA environment path for internal_data                       │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│       │                                                                      │
│       ▼                                                                      │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │  IgBLASTN Class (src/sadie/airr/igblast/igblast.py)                     │ │
│  │  - Constructs IgBLAST command-line arguments                             │ │
│  │  - Executes igblastn subprocess                                          │ │
│  │  - Returns AIRR format (outfmt=19) DataFrame                             │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│       │                                                                      │
│       ▼                                                                      │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │  AirrTable Class (src/sadie/airr/airrtable/airrtable.py)                │ │
│  │  - Validates AIRR compliance                                             │ │
│  │  - Calculates liability flags (J gene issues)                            │ │
│  │  - Computes mutation frequencies                                         │ │
│  │  - Assembles VDJ recombination strings                                   │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│       │                                                                      │
│       ▼                                                                      │
│  Output: AirrTable with segment positions, CDR/FWR boundaries, sequences    │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Segment Position Discovery

### 1. V Gene Position Discovery

**Mechanism**: BLAST alignment against V gene database

**Key Fields Produced**:
- `v_sequence_start`: Start position in query sequence (1-based)
- `v_sequence_end`: End position in query sequence (1-based)
- `v_germline_start`, `v_germline_end`: Positions in germline reference
- `v_cigar`: CIGAR string describing alignment

**Database Path Resolution** (GermlineData):
```python
self.v_gene_dir = Path(self.blast_dir.__str__() + "V")
# e.g., /path/igblast/Ig/internal_data/human/human_V
```

### 2. D Gene Position Discovery (Heavy Chain Only)

**Mechanism**: BLAST alignment with minimum match constraint

**Key Fields Produced**:
- `d_sequence_start`, `d_sequence_end`: Query positions
- `d_germline_start`, `d_germline_end`: Germline positions
- `d_cigar`: CIGAR string
- `d_call`: Best D gene match

**Configuration**:
```python
self.min_d_match = 5  # Minimum consecutive nucleotide matches for D
```

### 3. J Gene Position Discovery

**Mechanism**: BLAST alignment + auxiliary file lookup for CDR3 boundary

**Key Fields Produced**:
- `j_sequence_start`, `j_sequence_end`: Query positions
- `j_germline_start`, `j_germline_end`: Germline positions
- `j_cigar`: CIGAR string

**Critical Dependency**: Auxiliary file required for FWR4/CDR3 boundary calculation

### 4. C Region Position Discovery

**Mechanism**: BLAST alignment against C gene database (optional)

**Key Fields Produced**:
- `c_sequence_start`, `c_sequence_end`: Query positions
- `c_call`: Constant region gene match

**Database Configuration**:
```python
self.germline_db_c = self.germline_data.c_gene_dir
# Optional - warns if not found
```

## Auxiliary File Architecture

### Purpose
IgBLAST aux files provide metadata for J genes to determine CDR3 endpoint boundaries. Without correct aux data, CDR3 and FWR4 annotations fail.

### Format (5-column TSV)
```
<gene_name>\t<reading_frame>\t<chain_type>\t<cdr3_end>\t<is_functional>
```

### Example Entries
```
IGHJ1*01	0	JH	17	1
IGHJ2*01	1	JH	18	1
IGKJ1*01	1	JK	6	1
```

### Column Definitions

| Column | Description | Values |
|--------|-------------|--------|
| gene_name | Full allele name | IGHJ1*01, IGKJ2*03, etc. |
| reading_frame | Frame offset (0-2) | 0, 1, or 2 |
| chain_type | Chain identifier | JH, JK, JL |
| cdr3_end | Nucleotides from J start to CDR3 end | Integer (e.g., 17, 6) |
| is_functional | Functional gene flag | 1 (functional) or 0 (pseudogene) |

### Path Resolution
```python
# Legacy G3 path:
self.aux_path = self.base_dir / f"aux_db/{scheme}/{name}_gl.aux"

# New germlines module path:
self.aux_path = germlines_igblast / "aux_db" / f"{name}_gl.aux"
```

### CDR3 Boundary Calculation

IgBLAST uses aux file to calculate:
1. **CDR3 Start**: From V gene conserved cysteine (defined in internal_data/*.ndm files)
2. **CDR3 End**: J gene start position + `cdr3_end` value from aux file

**Why Aux File Matters**: The `cdr3_end` column tells IgBLAST how many nucleotides into the J gene the CDR3 extends. Without this, the CDR3/FWR4 boundary cannot be determined.

## Internal Data Architecture (NDM Files)

### Purpose
NDM (Numbering Data Map) files define framework and CDR boundaries for V genes based on IMGT or Kabat numbering schemes.

### Location
```
igdata/internal_data/{species}/{species}.ndm.{scheme}
```

### Contains
- V gene allele names
- Framework/CDR boundary positions
- Conserved residue positions (e.g., Cys at position 23 for CDR3 start anchor)

## Liability Detection

The AirrTable performs J gene liability detection to identify annotation failures:

```python
def _check_j_gene_liability(self, X: pd.Series) -> bool:
    """
    Liability conditions:
    1. Have CDR1/CDR2 but missing CDR3 → Likely aux file issue
    2. Have FWR3 + CDR3 but missing FWR4 → J gene not properly resolved
    3. Only FWR3 + CDR3 + FWR4 → Short read, acceptable
    """
```

**Liability Flag Usage**: When `liable=True`, SADIE can trigger adaptive penalty adjustment to retry annotation with different V/D/J penalties.

## Feature Flag: Germlines Module

```python
def _use_germlines_module() -> bool:
    env_value = os.environ.get("SADIE_USE_GERMLINES_MODULE", "true").lower()
    # "true" = Use new germlines module (default)
    # "false" = Use legacy G3 paths
```

This controls whether segment databases come from:
- **New path**: `src/sadie/germlines/igblast/`
- **Legacy path**: `src/sadie/airr/data/germlines/`

## IgBLAST Command Construction

Key arguments affecting segment discovery:

```python
cmd = [
    str(self.executable),
    "-germline_db_V", path_to_v_db,
    "-germline_db_D", path_to_d_db,
    "-germline_db_J", path_to_j_db,
    "-c_region_db", path_to_c_db,      # Optional
    "-auxiliary_data", path_to_aux,     # Critical for CDR3
    "-organism", species,
    "-domain_system", "imgt",           # or "kabat"
    "-outfmt", "19",                    # AIRR format
    "-V_penalty", v_penalty,
    "-D_penalty", d_penalty,
    "-J_penalty", j_penalty,
    "-min_D_match", min_d_match,
    "-extend_align5end",
    "-extend_align3end",
    "-query", input_fasta
]
```

## Data Flow Summary

```
1. Input FASTA → Airr.run_fasta()
2. GermlineData resolves all paths (V/D/J/C databases + aux file)
3. IgBLASTN constructs command with all paths
4. IgBLAST executes:
   - BLAST V/D/J/C alignments → segment positions
   - Reads aux file → CDR3 end position for J genes
   - Reads internal_data NDM → FWR/CDR boundaries for V genes
5. Returns AIRR format TSV (outfmt=19)
6. AirrTable validates and enriches:
   - Verifies CDR3/FWR4 presence (liability check)
   - Computes VDJ recombination string
   - Calculates mutation frequencies
7. Returns AirrTable with full segment annotations
```

## Critical Files for Segment Discovery

| Component | File | Purpose |
|-----------|------|---------|
| Aux file | `{species}_gl.aux` | J gene CDR3 end positions |
| NDM file | `{species}.ndm.imgt` | V gene FR/CDR boundaries |
| V database | `{species}_V.*` | V gene BLAST index |
| D database | `{species}_D.*` | D gene BLAST index |
| J database | `{species}_J.*` | J gene BLAST index |
| C database | `{species}_C.*` | C region BLAST index |
