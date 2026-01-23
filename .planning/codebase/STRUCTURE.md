# G3 Backend Code Structure

## Directory Layout

```
src/sadie/
├── airr/                           # Main AIRR annotation module
│   ├── airr.py                     # Airr class - orchestrator
│   ├── airrtable/                  # Output table handling
│   │   ├── airrtable.py            # AirrTable, AirrSeries, LinkedAirrTable
│   │   ├── constants.py            # AIRR column type definitions
│   │   └── genbank.py              # GenBank format export
│   ├── igblast/                    # IgBLAST integration
│   │   ├── igblast.py              # IgBLASTN wrapper class
│   │   └── germline.py             # GermlineData path resolution
│   ├── models/                     # Pydantic models
│   │   └── series.py               # AirrSeriesModel validation
│   ├── methods.py                  # Utility methods
│   ├── exceptions.py               # Custom exceptions
│   ├── data/germlines/             # Legacy G3 germline data
│   │   ├── aux_db/imgt/            # Auxiliary files (J gene data)
│   │   └── Ig/                     # BLAST databases
│   └── bin/                        # IgBLAST executables
│
├── germlines/                      # New germline management module
│   ├── __init__.py                 # Public API
│   ├── manager.py                  # GermlineManager class
│   ├── pipeline.py                 # Build pipeline orchestrator
│   ├── models.py                   # GermlineGene, ProviderMetadata
│   ├── builders/                   # Database builders
│   │   ├── aux.py                  # AuxFileBuilder class
│   │   ├── j_gene_data.py          # J gene reference data
│   │   ├── blast.py                # BLAST database builder
│   │   └── gapper.py               # Sequence gapping
│   ├── providers/                  # Data sources
│   │   ├── imgt.py                 # IMGT provider
│   │   ├── ogrdb.py                # OGRDB provider
│   │   ├── vdjbase.py              # VDJbase provider
│   │   └── custom.py               # Custom reference provider
│   ├── igblast/                    # Output databases
│   │   ├── Ig/internal_data/       # NDM files + BLAST DBs
│   │   └── aux_db/                 # Auxiliary files
│   └── scripts/                    # CLI utilities
│       ├── build_aux_files.py
│       └── build_internal_data.py
│
└── receptor/                       # Receptor chain models
    └── rearrangment.py             # ReceptorChain, AlignmentPositions, etc.
```

## Key Classes and Their Responsibilities

### AIRR Module (`src/sadie/airr/`)

#### `Airr` (airr.py:67)
Main entry point for AIRR annotation.

**Responsibilities**:
- Configure IgBLAST parameters (penalties, alignment options)
- Initialize GermlineData with correct paths
- Execute annotation via `run_fasta()`, `run_single()`, `run_records()`
- Handle adaptive penalty adjustment for failed annotations
- Apply allele coercion when exact matches unavailable

**Key Methods**:
```python
def run_fasta(file: Path, scfv: bool = False) -> AirrTable
def run_single(seq_id: str, seq: str) -> AirrTable
def run_records(seqrecords: List[SeqRecord]) -> AirrTable
def _apply_allele_coercion(result: AirrTable) -> AirrTable
```

**Segment Discovery Interface**:
```python
self.igblast.germline_db_v = self.germline_data.v_gene_dir
self.igblast.germline_db_d = self.germline_data.d_gene_dir
self.igblast.germline_db_j = self.germline_data.j_gene_dir
self.igblast.germline_db_c = self.germline_data.c_gene_dir
self.igblast.aux_path = self.germline_data.aux_path  # Critical for CDR3
```

---

#### `IgBLASTN` (igblast/igblast.py:188)
Low-level IgBLAST subprocess wrapper.

**Responsibilities**:
- Construct IgBLAST command-line arguments
- Execute igblastn and capture AIRR output
- Validate all required arguments before execution

**Key Properties**:
```python
@property
def aux_path(self) -> Path:
    """Auxiliary data for J gene CDR3 boundaries."""

@property
def germline_db_v(self) -> Path:
    """V gene BLAST database prefix."""

@property  
def germline_db_j(self) -> Path:
    """J gene BLAST database prefix."""
```

**Command Construction** (line 886):
```python
@property
def cmd(self) -> List[str]:
    _cmd = [str(self.executable)]
    for blast_arg in self.arguments:
        kv = blast_arg.get_formatted_blast_arg()
        if kv:
            _cmd += kv
    return _cmd
```

---

#### `GermlineData` (igblast/germline.py:28)
Resolves paths to germline databases and auxiliary files.

**Responsibilities**:
- Locate V/D/J/C BLAST database prefixes
- Find auxiliary file for species/scheme
- Set IGDATA environment for IgBLAST
- Support both legacy G3 and new germlines module paths

**Path Resolution Logic**:
```python
# Feature flag check
if _use_germlines_module():
    # New germlines module paths
    self.aux_path = germlines_igblast / "aux_db" / f"{name}_gl.aux"
else:
    # Legacy G3 paths
    self.aux_path = self.base_dir / f"aux_db/{scheme}/{name}_gl.aux"
```

**Critical Properties**:
```python
@property
def aux_path(self) -> Path:
    """Path to J gene auxiliary data for CDR3 reconstruction."""
    
@property
def j_gene_dir(self) -> Path:
    """J gene BLAST database prefix path."""
```

---

#### `AirrTable` (airrtable/airrtable.py:192)
DataFrame subclass with AIRR-specific functionality.

**Responsibilities**:
- Validate AIRR column compliance
- Compute liability flags (J gene annotation issues)
- Calculate mutation frequencies
- Build VDJ recombination strings
- Export to AIRR TSV, FASTA, GenBank formats

**Liability Detection** (line 650):
```python
def _check_j_gene_liability(self, X: pd.Series) -> bool:
    """
    Check CDR3/FWR4 annotation completeness.
    Returns True if J gene annotation likely failed.
    """
```

**Segment Position Fields**:
```python
# V segment
"v_sequence_start", "v_sequence_end"

# D segment  
"d_sequence_start", "d_sequence_end"

# J segment
"j_sequence_start", "j_sequence_end"

# C region
"c_sequence_start", "c_sequence_end"

# CDR3 boundaries
"cdr3_start", "cdr3_end"
```

---

### Germlines Module (`src/sadie/germlines/`)

#### `AuxFileBuilder` (builders/aux.py:24)
Builds IgBLAST auxiliary files from gapped sequences.

**Aux File Format** (5-column TSV):
```
IGHJ1*01	0	JH	17	1
<gene>\t<reading_frame>\t<chain_type>\t<cdr3_end>\t<is_functional>
```

**Key Methods**:
```python
def build_for_species(species: str, source_dir: Path, output_file: Path) -> None
def _create_aux_entry(record, chain: str, segment: str) -> Optional[str]
def validate_aux_file(aux_file: Path) -> bool
```

---

#### `get_j_gene_data` (builders/j_gene_data.py:58)
Lookup J gene reference data for aux file generation.

**Returns**: `(reading_frame, chain_type, cdr3_end, is_functional)`

**Reference Data Dictionary**:
```python
HUMAN_J_GENE_DATA = {
    "IGHJ1*01": (0, 17, 1),   # RF=0, CDR3 ends 17nt into J
    "IGHJ4*02": (2, 13, 1),   # RF=2, CDR3 ends 13nt into J
    "IGKJ1*01": (1, 6, 1),    # RF=1, CDR3 ends 6nt into J
}
```

---

## Critical Files for Segment Discovery

### Auxiliary Files (J Gene CDR3 Boundaries)

| Location | Purpose |
|----------|---------|
| `src/sadie/airr/data/germlines/aux_db/imgt/{species}_gl.aux` | Legacy G3 aux files |
| `src/sadie/germlines/igblast/aux_db/{species}_gl.aux` | New germlines aux files |

### Internal Data (NDM Files)

| Location | Purpose |
|----------|---------|
| `igblast/Ig/internal_data/{species}/{species}.ndm.imgt` | V gene FR/CDR boundaries |

### BLAST Databases

| Database | Location Pattern | Purpose |
|----------|------------------|---------|
| V genes | `{species}_V.n*` | V gene alignment |
| D genes | `{species}_D.n*` | D gene alignment |
| J genes | `{species}_J.n*` | J gene alignment |
| C region | `{species}_C.n*` | Constant region alignment |

---

## Where to Add New Code

### Adding a New Species

1. **Aux file**: Add entries to `src/sadie/germlines/builders/j_gene_data.py`
2. **Build databases**: Run `sadie.germlines.update_databases(species)`
3. **Verify**: Check `GermlineData.get_available_datasets()`

### Fixing J Gene CDR3 Issues

1. **Identify**: Check aux file format at `src/sadie/germlines/igblast/aux_db/{species}_gl.aux`
2. **Fix**: Update `HUMAN_J_GENE_DATA` in `builders/j_gene_data.py`
3. **Rebuild**: `AuxFileBuilder().build_for_species(species, source_dir, output_file)`

### Adding C Region Support

1. **Database**: Ensure C gene BLAST database exists at `{species}_C.*`
2. **GermlineData**: Verify `c_gene_dir` property resolves correctly
3. **IgBLAST**: Check `germline_db_c` argument is being passed

### Debugging Segment Position Issues

1. **Enable debug mode**: `Airr(species, debug=True)`
2. **Check command**: IgBLAST command is logged with all paths
3. **Verify aux file**: Ensure J gene allele exists with correct `cdr3_end` value
4. **Check liability**: `airr_table["liable"]` column indicates annotation failures

---

## Testing Entry Points

| Test File | Purpose |
|-----------|---------|
| `tests/data/fixtures/reference/igblast_aux/` | Aux file test fixtures |
| `src/sadie/germlines/tests/test_g3_regression.py` | G3 compatibility tests |
| `src/sadie/germlines/tests/test_integration.py` | Full pipeline tests |

---

## Configuration

### Environment Variables

```bash
SADIE_USE_GERMLINES_MODULE=true   # Use new germlines module (default)
SADIE_USE_GERMLINES_MODULE=false  # Use legacy G3 paths
```

### Airr Class Parameters Affecting Segment Discovery

```python
Airr(
    reference_name="human",
    v_gene_penalty=-1,          # V alignment penalty
    d_gene_penalty=-1,          # D alignment penalty  
    j_gene_penalty=-2,          # J alignment penalty
    min_d_match=5,              # Min D nucleotide matches
    extend_align5end=True,      # Extend V alignment 5'
    extend_align3end=True,      # Extend J alignment 3'
    scheme="imgt",              # Numbering scheme for boundaries
    coerce=False,               # Accept alternate alleles
)
```
