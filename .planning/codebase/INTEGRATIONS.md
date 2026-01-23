# G3 Backend Integrations

> How the G3 backend (non-germlines module) integrates with IgBLAST, databases, and auxiliary files for AIRR annotation.

## IgBLAST Integration

### Overview

G3 backend wraps NCBI IgBLAST (`igblastn`) for antibody sequence annotation. The wrapper handles:
- Command construction with proper arguments
- Environment variable setup (IGDATA)
- Output parsing to pandas DataFrames

### Integration Points

**File**: `src/sadie/airr/igblast/igblast.py`

```python
class IgBLASTN:
    """
    Python wrapper for igblastn command-line tool.
    
    Key integration:
    - Sets IGDATA environment variable for internal_data discovery
    - Constructs command with -germline_db_V, -germline_db_D, -germline_db_J, -c_region_db
    - Passes -auxiliary_data for J gene CDR3 boundaries
    """
```

### IgBLAST Parameters

| Parameter | IgBLAST Flag | SADIE Default | Purpose |
|-----------|--------------|---------------|---------|
| `germline_db_v` | `-germline_db_V` | species-specific | V gene BLAST database |
| `germline_db_d` | `-germline_db_D` | species-specific | D gene BLAST database |
| `germline_db_j` | `-germline_db_J` | species-specific | J gene BLAST database |
| `germline_db_c` | `-c_region_db` | species-specific | C region BLAST database |
| `aux_path` | `-auxiliary_data` | `{species}_gl.aux` | J gene CDR3 boundaries |
| `organism` | `-organism` | e.g., "human" | Species for internal_data lookup |
| `outfmt` | `-outfmt` | 19 | AIRR rearrangement format |
| `nomenclature` | `-domain_system` | "imgt" | Numbering scheme |
| `v_penalty` | `-V_penalty` | -1 | V gene mismatch penalty |
| `d_penalty` | `-D_penalty` | -1 | D gene mismatch penalty |
| `j_penalty` | `-J_penalty` | -2 | J gene mismatch penalty |

### IGDATA Environment Variable

IgBLAST requires `IGDATA` to point to the directory containing `internal_data/{species}/`:

```python
# src/sadie/airr/igblast/igblast.py
def run_file(self, file: Union[Path, str]) -> pd.DataFrame:
    local_env = os.environ.copy()
    local_env["IGDATA"] = str(self.igdata)  # Points to Ig/ directory
    process = subprocess.run(cmd, env=local_env, capture_output=True)
```

**G3 Legacy IGDATA**: `src/sadie/airr/data/germlines/Ig/`

---

## J Gene Segment Position Discovery

### How J Gene Boundaries Are Determined

1. **Auxiliary File Lookup**: IgBLAST reads the `.aux` file to find CDR3 end position
2. **Alignment**: J gene aligned to query sequence
3. **CDR3 Calculation**: CDR3 end = J alignment start + cdr3_end from aux file

### Auxiliary File Format

**Location (G3)**: `src/sadie/airr/data/germlines/aux_db/imgt/{species}_gl.aux`

```
# Format: gene_name	reading_frame	chain_type	cdr3_end	is_functional
IGHJ1*01	0	JH	17	1
IGHJ2*01	1	JH	18	1
```

### Column Definitions

| Column | Description |
|--------|-------------|
| `gene_name` | J gene allele name (e.g., `IGHJ1*01`) |
| `reading_frame` | Position within codon (0, 1, or 2) |
| `chain_type` | Chain identifier (JH, JK, JL) |
| `cdr3_end` | Nucleotide position where CDR3 ends (relative to J gene start) |
| `is_functional` | 1 = functional, 0 = pseudogene |

### J Gene Reference Data

**File**: `src/sadie/germlines/builders/j_gene_data.py` (germlines module reference)

Contains hardcoded CDR3 end positions sourced from IMGT:

```python
HUMAN_J_GENE_DATA = {
    "IGHJ1*01": (0, 17, 1),  # (reading_frame, cdr3_end, is_functional)
    "IGHJ2*01": (1, 18, 1),
    "IGKJ1*01": (1, 6, 1),
    # ...
}
```

### CDR3 End Position Logic

The `cdr3_end` value represents the nucleotide position relative to the start of the J gene where the CDR3 region ends:

- **Heavy chain (JH)**: Typically 13-28 nt (varies by allele)
- **Kappa chain (JK)**: Typically 6-7 nt
- **Lambda chain (JL)**: Typically 6 nt

---

## Constant Region Position Discovery

### How C Region Is Discovered

1. **BLAST Database**: Separate C region database (`{species}_C.*`)
2. **IgBLAST Query**: Uses `-c_region_db` parameter
3. **Alignment**: C region aligned after J gene

### C Region Database Files

**Location (G3)**: `src/sadie/airr/data/germlines/Ig/blastdb/{species}/{species}_C.*`

```
human_C.fasta    # Source sequences
human_C.ndb      # BLAST database index
human_C.nhi      # Header index
human_C.nhr      # Header file
human_C.nin      # Index file
human_C.nsq      # Sequence file
...
```

### C Region Output Fields

From `src/sadie/airr/airrtable/constants.py`:

```python
CONSTANTS_AIRR = {
    "c_call": "object",           # C gene with allele
    "c_cigar": "object",          # Alignment CIGAR string
    "c_germline_start": "Int16",  # Start in germline
    "c_germline_end": "Int16",    # End in germline
    "c_sequence_start": "Int16",  # Start in query sequence
    "c_sequence_end": "Int16",    # End in query sequence
    "c_identity": "float32",      # Identity percentage
    "c_score": "float32",         # Alignment score
    # ...
}
```

---

## V Gene Region Boundaries

### NDM File Integration

IgBLAST uses `.ndm.imgt` files for V gene FWR/CDR boundaries.

**Location (G3)**: `src/sadie/airr/data/germlines/Ig/internal_data/{species}/{species}.ndm.imgt`

**Format** (12 columns):
```
IGHV1-18*01	1	75	76	99	100	150	151	174	175	288	VH	0
```

| Column | Description |
|--------|-------------|
| 1 | Gene name |
| 2 | FWR1 start |
| 3-4 | CDR1 start-end |
| 5-6 | FWR2 start-end |
| 7-8 | CDR2 start-end |
| 9-10 | FWR3 start-end |
| 11 | V gene end (total length) |
| 12 | Chain type (VH, VK, VL) |
| 13 | Flags |

---

## G3 API Integration

### REST API (Legacy)

**File**: `src/sadie/reference/reference.py`

```python
class Reference:
    _endpoint = "https://g3.jordanrwillis.com/api/v1/genes"
    
    def _g3_get(self, query: str) -> Tuple[int, List[Dict[str,str]]]:
        """Query G3 REST API for gene data"""
        response = requests.get(query)
        return response.status_code, response.json()
```

### Query Format

```
GET https://g3.jordanrwillis.com/api/v1/genes?source=imgt&common=human&gene=IGHV1-69*01
```

### Response Fields

G3 API returns gene data including:
- Gene name and allele
- Nucleotide sequence
- Gapped sequence (IMGT numbering)
- CDR/FWR positions
- Species information

---

## Data Flow Summary

```
┌─────────────────────────────────────────────────────────────────┐
│                    User Query (FASTA file)                      │
└─────────────────────────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│                     Airr Class (airr.py)                        │
│  - Selects GermlineData based on species                        │
│  - Configures IgBLASTN with paths                               │
└─────────────────────────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│               GermlineData Class (germline.py)                  │
│  - Resolves paths based on SADIE_USE_GERMLINES_MODULE flag      │
│  - Returns: v_gene_dir, j_gene_dir, c_gene_dir, aux_path        │
└─────────────────────────────────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────┐
│                  IgBLASTN Class (igblast.py)                    │
│  - Sets IGDATA environment variable                             │
│  - Constructs igblastn command with all parameters              │
│  - Executes subprocess and parses TSV output                    │
└─────────────────────────────────────────────────────────────────┘
                               │
         ┌─────────────────────┼─────────────────────┐
         ▼                     ▼                     ▼
┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
│  BLAST DBs      │  │   NDM Files     │  │   AUX Files     │
│  {species}_V    │  │ V gene CDR/FWR  │  │ J gene CDR3 end │
│  {species}_D    │  │   boundaries    │  │   boundaries    │
│  {species}_J    │  │                 │  │                 │
│  {species}_C    │  │                 │  │                 │
└─────────────────┘  └─────────────────┘  └─────────────────┘
```

---

## Key Integration Files

| File | Integration Purpose |
|------|---------------------|
| `src/sadie/airr/igblast/igblast.py` | IgBLAST command execution |
| `src/sadie/airr/igblast/germline.py` | Database path resolution |
| `src/sadie/airr/airr.py` | High-level annotation API |
| `src/sadie/reference/reference.py` | G3 API client |
| `src/sadie/reference/yaml.py` | Reference YAML parser |
| `src/sadie/airr/airrtable/constants.py` | AIRR output column definitions |

---

## Adaptive Penalty Correction

When initial annotation fails (marked "liable"), the G3 backend can retry with different penalty settings:

```python
# src/sadie/airr/airr.py
if self.adapt_penalty:
    for v_penalty in range(-2, -4, -1):
        for j_penalty in range(-1, -3, -1):
            # Retry annotation with relaxed penalties
            adaptable_api = Airr(
                self.name,
                v_gene_penalty=v_penalty,
                j_gene_penalty=j_penalty,
                adaptable=False,  # Prevent infinite recursion
            )
```

This helps recover annotations when default penalties are too strict.

---

## Allele Coercion

When exact allele match not found in auxiliary files:

```python
# src/sadie/airr/airr.py
def _apply_allele_coercion(self, result: AirrTable) -> AirrTable:
    """
    Accept highest-scored available allele when exact match not found.
    
    1. Read NDM file for available alleles
    2. For each call, check if top allele is available
    3. If not, substitute with best available from same gene family
    4. Add comment documenting the coercion
    """
```
