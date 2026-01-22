# Data Model: Germline Database Integration

**Feature**: 002-germline-integration
**Date**: 2026-01-20

## Overview

This integration leverages existing data models from the germlines module (001-germline-completion). No new entity types are created; instead, we define integration contracts and adapter patterns.

## Existing Entities (from germlines module)

### GermlineGene

**Source**: `src/sadie/germlines/models.py`

**Purpose**: Unified representation of germline gene across all providers

**Key Fields**:
```python
class GermlineGene(BaseModel):
    # Core identifiers
    name: str                                    # e.g., "IGHV1-69*01"
    species: str                                 # e.g., "human"
    segment: str                                 # "V", "D", or "J"
    chain: str                                   # "H", "K", or "L"

    # Sequences
    sequence: str                                # Ungapped nucleotide (REQUIRED)
    sequence_gapped: Optional[str]               # IMGT-gapped nucleotide (REQUIRED for V/J if sequence_aa_gapped missing)
    sequence_aa: Optional[str]                   # Ungapped amino acid
    sequence_aa_gapped: Optional[str]            # IMGT-gapped amino acid (REQUIRED for V/J unless sequence_gapped available)

    # Functional annotation
    is_functional: bool = True
    functionality: str = "F"                     # "F", "ORF", or "P"

    # IMGT regions (CDR/FWR)
    regions: Optional[Dict[str, str]]            # e.g., {"cdr1": "GYTF...", ...}
    region_positions: Optional[Dict[str, Tuple[int, int]]]

    # Source tracking
    source: str                                  # "imgt", "ogrdb", "vdjbase", "custom"
    source_version: Optional[str]

    # Metadata
    allele: Optional[str]
    gene_family: Optional[str]
    accession: Optional[str]
```

**Validation Rules**:
- `segment` must be "V", "D", or "J"
- `chain` must be "H", "K", or "L"
- `sequence` must contain only valid nucleotides (ACGTN)
- `functionality` must be "F", "ORF", or "P"

**Gapped Sequence Availability (V/J genes)**:
- For V and J segment genes used in HMM building, at least ONE of `sequence_aa_gapped` OR `sequence_gapped` MUST be present
- If `sequence_aa_gapped` is missing, `sequence_gapped` is used for runtime translation via `LocalHMMBuilder._translate_gapped_nt_to_aa()`
- Validation occurs at HMM build time in `LocalHMMBuilder._get_vj_alignment_pairs()`
- Clear error message if neither is available: "Gene {name} lacks gapped sequence data for HMM building"

**Used By**:
- IgBLAST integration (nucleotide sequences)
- HMM builder (gapped amino acid sequences)
- Reference system adapter (all fields)

### GermlineManager

**Source**: `src/sadie/germlines/manager.py`

**Purpose**: Central access point for querying germline data

**Key Methods**:
```python
class GermlineManager:
    def get_genes(
        self,
        species: str,
        segment: str,
        chain: str,
        provider: Optional[str] = None
    ) -> List[GermlineGene]:
        """Query genes with optional provider filter."""

    def get_gene_by_name(
        self,
        species: str,
        name: str
    ) -> Optional[GermlineGene]:
        """Lookup single gene by name."""
```

**Used By**:
- LocalHMMBuilder (query V/J genes for HMM generation)
- GermlineToG3Adapter (query genes for Reference system)
- IgBLAST integration (indirect via database files)

## Integration-Specific Data Structures

### G3APIFormat (Dict)

**Purpose**: G3 API-compatible dictionary format for Reference system

**Structure**:
```python
G3APIFormat = Dict[str, Any]  # Shape defined below

{
    "source": str,              # Provider name
    "common": str,              # Common species name
    "latin": str,               # Scientific name (e.g., "Homo_sapiens")
    "gene": str,                # Gene name (e.g., "IGHV1-69*01")
    "label": str,               # Segment label (e.g., "V-REGION")
    "gene_segment": str,        # "V", "D", or "J"
    "receptor": str,            # "IG" or "TR"
    "sequence": str,            # Ungapped nucleotide
    "species": str,             # Common name (for compatibility)

    # Nested IMGT structure
    "imgt": {
        "sequence": str,
        "sequence_gapped": str,
        "sequence_gapped_aa": str,
        "imgt_functional": str,              # "F", "ORF", "P"
        "contrived_functional": str,

        # Regions (optional, V genes only)
        "fwr1": str, "fwr1_aa": str,
        "fwr1_start": int, "fwr1_end": int,
        "cdr1": str, "cdr1_aa": str,
        "cdr1_start": int, "cdr1_end": int,
        "fwr2": str, "fwr2_aa": str,
        "fwr2_start": int, "fwr2_end": int,
        "cdr2": str, "cdr2_aa": str,
        "cdr2_start": int, "cdr2_end": int,
        "fwr3": str, "fwr3_aa": str,
        "fwr3_start": int, "fwr3_end": int,
        "cdr3": str, "cdr3_aa": str,
        "cdr3_start": int, "cdr3_end": int,
        "fwr4": str, "fwr4_aa": str,         # J genes only
        "fwr4_start": int, "fwr4_end": int
    }
}
```

**Transformation**: `GermlineGene` → `G3APIFormat` via `GermlineToG3Adapter`

**Used By**: Reference system (`reference/reference.py`)

### IgBLASTDatabasePaths

**Purpose**: File paths for IgBLAST database components

**Structure**:
```python
@dataclass
class IgBLASTDatabasePaths:
    base_dir: Path              # Base directory for databases
    v_gene_dir: Path            # Prefix for V gene BLAST DB
    d_gene_dir: Path            # Prefix for D gene BLAST DB
    j_gene_dir: Path            # Prefix for J gene BLAST DB
    c_gene_dir: Path            # Prefix for C gene BLAST DB (optional)
    aux_path: Path              # Auxiliary file for J genes
    igdata: Path                # Internal data directory
```

**Source**: Derived from `GermlineData` class in `airr/igblast/germline.py`

**Conditional Logic**:
```python
if use_germlines_module():
    paths = germlines_module_paths(species)
else:
    paths = legacy_g3_paths(species)
```

**Used By**: IgBLAST wrapper in AIRR annotation

### HMMCacheKey

**Purpose**: Cache key for HMM models

**Structure**:
```python
@dataclass(frozen=True)
class HMMCacheKey:
    species: str
    chain: str
    source: str = "imgt"

    def to_filename(self) -> str:
        return f"{self.species}_{self.chain}.hmm"
```

**Used By**: LocalHMMBuilder for HMM caching

## Data Flow Diagrams

### AIRR Annotation Flow

```
User Request (sequences)
    ↓
AIRR.annotate()
    ↓
Check SADIE_USE_GERMLINES_MODULE
    ↓
├─[true]─→ GermlineData(species)
│             ↓
│         germlines/igblast/database/
│             ↓
│         IgBLAST execution
│
└─[false]─→ GermlineData(species)
              ↓
          airr/data/germlines/
              ↓
          IgBLAST execution
```

### Renumbering Flow

```
User Request (sequences)
    ↓
HMMER.get_hmm_models()
    ↓
Check use_germlines_module()
    ↓
├─[true]─→ LocalHMMBuilder.get_hmm(species, chain)
│             ↓
│         Check cache: germlines/hmms/{species}_{chain}.hmm
│         │
│         ├─[exists]─→ Load cached HMM
│         │
│         └─[missing]─→ GermlineManager.get_genes(species, V/J, chain)
│                         ↓
│                     Build Stockholm alignment
│                         ↓
│                     pyhmmer.build_msa()
│                         ↓
│                     Cache HMM binary
│
└─[false]─→ G3.get_hmm() [legacy]
```

### Reference System Flow

```
Reference.add_gene(species, gene_name)
    ↓
Check use_germlines parameter
    ↓
├─[true]─→ GermlineManager.get_gene_by_name(species, name)
│             ↓
│         GermlineToG3Adapter.to_g3_format(gene)
│             ↓
│         Return G3APIFormat dict
│
└─[false]─→ HTTP GET to G3 API [legacy]
```

## Entity Relationships

```
GermlineManager
    ├─→ [manages] → GermlineGene (many)
    └─→ [queries] → Providers (IMGT, OGRDB, VDJbase, custom)

GermlineGene
    ├─→ [transformed by] → GermlineToG3Adapter
    │                           ↓
    │                      G3APIFormat (dict)
    │                           ↓
    │                      Reference system
    │
    ├─→ [used by] → LocalHMMBuilder
    │                   ↓
    │              HMM binaries
    │                   ↓
    │              HMMER (renumbering)
    │
    └─→ [used by] → BlastDBBuilder (from 001-germline-completion)
                        ↓
                   BLAST databases
                        ↓
                   IgBLAST (AIRR annotation)
```

## Validation Rules

### Integration-Level Validation

**IgBLAST Integration**:
- Database files must exist at expected paths
- Prefixes validated via `ensure_prefix_to()` (checks .nhr, .nin, .nsq)
- Auxiliary file must exist for J gene CDR3 annotation
- Clear error message if databases not built

**HMM Integration**:
- Minimum 3 V or J genes required for HMM building
- Gapped amino acid sequences required
- Stockholm format validation
- Graceful fallback to G3 on build failure

**Reference Integration**:
- Species must be available in germlines module
- Gene name must exist in selected provider
- G3 format validation (all required fields present)
- Provider priority order respected

## State Transitions

### Feature Flag State Machine

```
[Default State: G3 Mode]
    ↓
SADIE_USE_GERMLINES_MODULE=true
    ↓
[Germlines Mode]
    ├─→ IgBLAST uses germlines/igblast/ paths
    ├─→ HMM uses LocalHMMBuilder
    └─→ Reference uses GermlineManager + adapter

    ↓ (on error)

[Fallback State]
    └─→ Log warning, fall back to G3 (if available)
```

### HMM Cache State Machine

```
[HMM Requested]
    ↓
Check germlines/hmms/{species}_{chain}.hmm
    ↓
├─[exists]─→ [Cached State] → Load and return
│
└─[missing]─→ [Build State]
               ↓
           Query germlines
               ↓
           Build Stockholm
               ↓
           Generate HMM
               ↓
           [Cached State]
```

## Summary

This integration leverages existing data models from the germlines module without introducing new entities. The key integration patterns are:

1. **Adapter**: GermlineGene → G3APIFormat transformation
2. **Builder**: LocalHMMBuilder for on-demand HMM generation
3. **Strategy**: Feature flag for runtime backend selection

All data validation is inherited from the germlines module. Integration adds path validation, format conversion, and error handling layers.

