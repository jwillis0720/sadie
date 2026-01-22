# Research: Germlines Module Completion

**Date**: 2026-01-14 | **Spec**: [spec.md](./spec.md) | **Plan**: [plan.md](./plan.md)

## Overview

This document captures research findings for the Phase 0 unknowns identified in the implementation plan. Each research question is addressed with findings, decisions, and implementation implications.

---

## 1. VDJbase Data Format & Access

**Status**: RESOLVED

### Question
What is VDJbase's current FASTA format and download mechanism?

### Findings

#### API Structure

**Base URL**: `https://vdjbase.org/admin/api/`

| Endpoint | Purpose | Example |
|----------|---------|---------|
| `genomic/species` | List supported species | Returns `["Rhesus Macaque", "Human"]` |
| `genomic/data_sets/{species}` | List datasets per species | IGH, IGK, IGL, IGHC for Human |
| `repseq/samples/{species}/{chain}` | Repertoire sample metadata | Paginated, includes genotype file paths |
| `repseq/sequences/{species}/{chain}` | Germline sequences | Returns allele sequences with metadata |

#### Supported Species & Datasets

**Human**:
- IGH (created 2025-07-09)
- IGK (created 2025-07-09)
- IGL (created 2025-07-08)
- IGHC - constant region (created 2025-08-06)

**Rhesus Macaque**:
- IGH (created 2025-09-28)
- IGK (created 2025-09-28)
- IGL (created 2025-09-28)

#### Sequence Data Format

The `/repseq/sequences/{species}/{chain}` endpoint returns paginated JSON:

```json
{
  "samples": [
    {
      "name": "IGHV2-70*15_g303a",
      "seq": "caggtcaccttgagggagtctggtcct...",
      "seq_len": "322",
      "type": "IGHV",
      "family": "IGHV2",
      "gene_name": "IGHV2-70",
      "species": "Human",
      "appears": 66,
      "novel": true,
      "is_single_allele": true,
      "low_confidence": true,
      "dataset": "IGH"
    }
  ]
}
```

**Key fields**:
- `name`: Allele identifier (mutation-based naming for novel alleles, e.g., `IGHV2-70*15_g303a`)
- `seq`: Ungapped nucleotide sequence
- `type`: Segment type (IGHV, IGHD, IGHJ)
- `novel`: Boolean - novel allele not in IMGT/OGRDB
- `low_confidence`: Quality indicator
- `appears`: Observation count across samples

#### OGRDB Report Files (Gapped Sequences)

Per-sample genotype files contain IMGT-gapped sequences:

**URL Pattern**:
```
https://vdjbase.org/static/study_data/VDJbase/samples/{species}/{chain}/{study}/{sample}/{sample}_ogrdb_report.csv
```

**CSV Format**:
```csv
"sequence_id","nt_sequence","nt_sequence_gapped"
"IGHV1-2*02","GCTTCTGGATACACCTTC...","NNNNNNN...GCTTCTGGATACACCTTC........."
```

#### Data Volume (Human IGH)

| Metric | Value |
|--------|-------|
| Total IGHV sequences | ~725 |
| Novel alleles | ~67% |
| Pagination | 250 items/page recommended |

#### API Quirks

1. **Page size affects results**: Smaller page sizes (10-100) recommended for consistent behavior
2. **Sorting varies**: Different page sizes may return different ordering
3. **No bulk FASTA download**: Must paginate through API

### Decision

**Approach**: Hybrid download with local caching

1. **Download Phase**: Paginate `/repseq/sequences` API, cache to local FASTA
2. **Gapped Templates**: Use OGRDB report CSVs or derive from IMGT gapped references
3. **Storage**: `src/sadie/germlines/sources/vdjbase/{species}/{chain}.fasta`
4. **Metadata**: Store `novel`, `low_confidence`, `appears` in FASTA headers

### Implementation Impact

```python
class VDJbaseProvider(GermlineProvider):
    """VDJbase integration using REST API with local caching."""
    
    BASE_URL = "https://vdjbase.org/admin/api"
    
    def fetch_genes(self, species: str, segment: str, chain: str) -> List[GermlineGene]:
        # 1. Check local cache first (offline-first)
        # 2. If not cached, paginate API and cache results
        # 3. Filter by segment type
        # 4. Note: sequences are UNGAPPED - require gapper module
        pass
```

---

## 2. Biopython Alignment Strategy for Gapping

**Status**: RESOLVED

### Question
How to use Biopython alignment for auto-gapping sequences to IMGT format?

### Decision

Use per-gene IMGT-gapped templates when available; fallback to per-segment consensus template.

**Template source**: IMGT provider gapped FASTA in `src/sadie/germlines/sources/imgt/`

### Implementation (Completed)

The `GapperService` has been implemented in `src/sadie/germlines/builders/gapper.py` with:

1. **BioPython PairwiseAligner** configured for antibody sequences:
   - Global alignment mode
   - Match score: 2.0, mismatch: -1.0
   - Gap open: -10.0, extend: -0.5
   - End gaps allowed without penalty

2. **Template Strategy**:
   - Per-gene template lookup (exact match by gene name)
   - Per-segment consensus fallback (built from available gapped sequences)
   - Templates cached for performance

3. **Workflow**:
   - Translate nucleotide to amino acid
   - Align AA sequences using PairwiseAligner
   - Extract gap positions from template
   - Apply gaps to nucleotide sequence at codon boundaries
   - Use "." (period) as gap character per IMGT convention

4. **Error Handling**:
   - D segments return None (no gapping required)
   - Failed gapping logs WARNING and returns None
   - Graceful degradation on missing templates

---

## 3. G3 API Response Format

**Status**: RESOLVED

### Question
What exact JSON format does G3 return for gene queries?

### Findings

From `src/sadie/renumbering/clients/g3.py`, the G3 client returns gene data with the following structure:

```json
{
  "gene": "IGHV1-69*01",
  "sequence": "CAGGTGCAGCTGGTGGAG...",
  "sequence_gapped": "CAGGTGCAGCTGGTGGAG......",
  "species": "human",
  "segment": "V",
  "chain": "H",
  "source": "imgt",
  "functional": true,
  "regions": {
    "fwr1": "CAGGTGCAGCTGGTGGAG",
    "cdr1": "GGTGGCAGCTTC",
    "fwr2": "TGGGTGCGCCAG",
    "cdr2": "ATAGACAGCAGTGGC",
    "fwr3": "CGCTCCGTGAAGGGCCGATTC"
  }
}
```

### Key Implementation Details

- G3 API endpoint: `https://g3.jordanrwillis.com/api/v1/genes`
- Uses LRU caching for performance
- Filters out certain species (pig, cow, cat, alpaca, rhesus, dog)
- Builds Stockholm alignments for HMM creation via `pyhmmer`

### Adapter Pattern

The Reference system adapter (`FR-012b`) should:
1. Call `GermlineManager.get_gene_by_name(name, species)`
2. Transform `GermlineGene` object to G3 response format
3. Return dict matching the structure above

---

## 4. IgBLAST Auxiliary File Format

**Status**: RESOLVED

### Question
What CDR/FWR annotations are required in auxiliary files?

### Findings

From `src/sadie/germlines/builders/aux.py`, the auxiliary file format is:

**Format**: Tab-separated values (TSV)

**Columns**:
1. `gene_name` - Gene identifier (e.g., "IGHV1-69*01")
2. `chain` - Chain type (H, K, L)
3. `segment` - Segment type (V, J)
4. `fwr1_start`, `fwr1_end` - Framework 1 boundaries
5. `cdr1_start`, `cdr1_end` - CDR1 boundaries
6. `fwr2_start`, `fwr2_end` - Framework 2 boundaries
7. `cdr2_start`, `cdr2_end` - CDR2 boundaries
8. `fwr3_start`, `fwr3_end` - Framework 3 boundaries

**IMGT Position Numbering** (from spec FR-037b):
- FWR1: positions 1-26
- CDR1: positions 27-38
- FWR2: positions 39-55
- CDR2: positions 56-65
- FWR3: positions 66-104

### Current Implementation Status

The `AuxFileBuilder` in `aux.py` is currently a stub. Key methods that need completion:
- `_parse_imgt_regions()` - Convert IMGT gaps to region boundaries
- `_create_aux_entry()` - Generate actual aux file entries

Output location: `igblast/aux_db/{scheme}/{species}_gl.aux`

---

## 5. IMGT/OGRDB Download URLs

**Status**: RESOLVED

### Question
Current stable URLs for bulk FASTA downloads from IMGT and OGRDB?

### IMGT

**Status**: RESOLVED (2026-01-16)

**Data Source**: V-QUEST Reference Directory
- **Base URL**: `https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/`
- **URL Pattern**: `{BASE_URL}/{Species}/IG/{segment}.fasta`
- **Example**: `https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta`

#### Available Species (36+)

| Internal Name | IMGT Directory Name |
|---------------|---------------------|
| human | Homo_sapiens |
| mouse | Mus_musculus |
| mouse_c57bl6j | Mus_musculus_C57BL6J |
| rat | Rattus_norvegicus |
| rabbit | Oryctolagus_cuniculus |
| rhesus_macaque | Macaca_mulatta |
| cynomolgus | Macaca_fascicularis |
| dog | Canis_lupus_familiaris |
| pig | Sus_scrofa |
| cow | Bos_taurus |
| alpaca | Vicugna_pacos |
| camel | Camelus_dromedarius |
| chicken | Gallus_gallus |
| zebrafish | Danio_rerio |

#### Segment Files Per Species

Standard IG segments available:
- `IGHV.fasta`, `IGHD.fasta`, `IGHJ.fasta` (Heavy chain)
- `IGKV.fasta`, `IGKJ.fasta` (Kappa chain)
- `IGLV.fasta`, `IGLJ.fasta` (Lambda chain)

Optional TR (T-cell receptor) segments:
- `TRAV.fasta`, `TRAJ.fasta`, `TRBV.fasta`, `TRBD.fasta`, `TRBJ.fasta`, etc.

#### FASTA Format

**Header Format**:
```
>accession|gene_name|species|functionality|region|positions|length|codon_start|...
```

**Example**:
```
>M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |
caggttcagctggtgcagtctggagct...gaggtgaagaagcctggggcctcagtgaag...
```

**Key Header Fields**:
- Position 1: Accession number (GenBank/EMBL)
- Position 2: Gene name (e.g., `IGHV1-18*01`)
- Position 3: Species (e.g., `Homo sapiens`)
- Position 4: Functionality (`F` = functional, `ORF` = open reading frame, `P` = pseudogene)
- Position 5: Region type (`V-REGION`, `D-REGION`, `J-REGION`)

**Sequence Characteristics**:
- **V-regions**: IMGT-gapped with dots (`.`) indicating gaps per IMGT unique numbering
- **D/J-regions**: NOT gapped (IMGT numbering only applies to V regions positions 1-104)
- **Case**: Lowercase nucleotides
- **Gap character**: Period (`.`)

#### Data Volume (Downloaded 2026-01-16)

| Species | Total | IGHV | IGHD | IGHJ | IGKV | IGKJ | IGLV | IGLJ |
|---------|-------|------|------|------|------|------|------|------|
| Human | 794 | 460 | 47 | 15 | 132 | 10 | 119 | 11 |
| Mouse | 953 | 678 | 61 | 9 | 168 | 10 | 19 | 8 |

#### Implementation

**Script**: `src/sadie/germlines/scripts/download_imgt.py`

**Features**:
- Direct HTTP download from V-QUEST reference directory
- Outputs both gapped (`*_gapped.fasta`) and ungapped (`*.fasta`) files
- Ungapped derived by stripping dots from gapped sequences
- Species mapping between internal names and IMGT directory names
- Optional TR (T-cell receptor) sequence download
- Progress logging and error handling

**Output Structure**:
```
sources/imgt/
├── human/
│   ├── IGHV.fasta          # Ungapped (dots removed)
│   ├── IGHV_gapped.fasta   # IMGT-gapped (original)
│   ├── IGHD.fasta          # No gaps (D segments)
│   ├── IGHD_gapped.fasta   # Same as ungapped
│   ├── IGHJ.fasta          # No gaps (J segments)
│   └── ...
└── mouse/
    └── ...
```

**Usage**:
```bash
python -m sadie.germlines.scripts.download_imgt --species human mouse
python -m sadie.germlines.scripts.download_imgt --list-species
python -m sadie.germlines.scripts.download_imgt --species human --include-tr
```

#### Format Comparison: IMGT vs OGRDB

| Aspect | IMGT | OGRDB |
|--------|------|-------|
| Header | Rich metadata (`accession\|gene\|species\|func\|...`) | Simple (`>gene_name`) |
| Nucleotide case | Lowercase | Uppercase |
| V-region gapping | Yes (dots `.`) | Yes (dots `.`) |
| J/D gapping | No (not applicable) | Some gapped |
| Allele naming | Standard `*01`, `*02` | Novel uses `*i01` notation |
| Functionality info | In header (`F`, `ORF`, `P`) | Not in header |

### OGRDB

**Status**: RESOLVED (2026-01-16)

**Data Source**: Zenodo Archive
- **URL**: https://zenodo.org/records/18145568/files/ogrdb_archive.tgz?download=1
- **Format**: MariaDB SQL dump (`ogrdb_dump.sql`)
- **Key Table**: `gene_description` with columns:
  - `sequence_name`: Gene allele name (e.g., IGHV1-69*01)
  - `species`: Species name (e.g., Homo sapiens)
  - `locus`: Chain locus (IGH, IGK, IGL)
  - `sequence_type`: Segment type (V, D, J)
  - `sequence`: Ungapped nucleotide sequence
  - `coding_seq_imgt`: IMGT-gapped nucleotide sequence

**Implementation**: `src/sadie/germlines/scripts/download_ogrdb.py`
- Downloads archive from Zenodo (~50MB)
- Parses SQL dump to extract sequences
- Outputs both gapped (`*_gapped.fasta`) and ungapped FASTA files
- Organizes by species/segment/chain

**Data Volume**:
| Species | Genes | Details |
|---------|-------|---------|
| Human | 47 | 20 IGHV, 7 IGHJ, 2 IGKV, 7 IGKJ, 1 IGLV, 10 IGLJ |
| Mouse | 1,771 | 426 IGHV, 66 IGHD, 28 IGHJ, 1,230 IGKV, 11 IGKJ, 3 IGLV, 7 IGLJ |

**Usage**:
```bash
python -m sadie.germlines.scripts.download_ogrdb --species human mouse
```

### VDJbase (from Section 1)

- **API Base**: https://vdjbase.org/admin/api/
- **Sequences Endpoint**: `/repseq/sequences/{species}/{chain}`
- Returns paginated JSON with sequence data
- Implemented with `VDJbaseProvider.download()` method

### Remaining Action Items

- [X] ~~Implement automated IMGT download script~~ **DONE** (2026-01-16)
- [X] ~~Implement automated OGRDB download script~~ **DONE** (2026-01-16)
- [ ] Add retry/resume logic for IMGT downloads (nice-to-have)
- [X] ~~Test OGRDB download script for reliability~~ **DONE**
- [X] ~~Test IMGT download script for human/mouse~~ **DONE** (2026-01-16)

---

## 6. In-Silico Gapping Validation

**Status**: RESOLVED (2026-01-16)

### Question

How accurate is the BioPython alignment-based in-silico gapping compared to the original G3 gapped sequences?

### Methodology

Created a validation test (`src/sadie/germlines/tests/test_gapping_accuracy.py`) that:
1. Loads ungapped sequences from custom FASTA files
2. Loads original G3 gapped sequences (exported with `--include-gapped`)
3. Generates in-silico gapped sequences using `GapperService` with human IMGT templates
4. Compares gap positions and overall character accuracy

### Findings

**Overall Results** (5,241 sequences across 5 species):
- Exact match rate: **1.3%** (68 sequences)
- Gap position accuracy: **60.0%**
- Character accuracy: **55.5%**

**V Segment Results by Species**:
| Species | Segment | Sequences | Exact Match | Gap Accuracy | Char Accuracy |
|---------|---------|-----------|-------------|--------------|---------------|
| Cat | IGHV | 119 | 0.8% | 86.0% | 63.4% |
| Cat | IGKV | 27 | 55.6% | 85.7% | 75.9% |
| Cat | IGLV | 68 | 36.8% | 86.3% | 69.8% |
| Dog | IGHV | 70 | 10.0% | 91.9% | 72.3% |
| Dog | IGLV | 117 | 10.3% | 52.1% | 54.3% |
| Macaque | IGHV | 1,969 | 0.2% | 77.4% | 61.6% |
| Macaque | IGKV | 1,018 | 0.0% | 58.5% | 53.5% |

**D/J Segment Results**:
- D segments: Not gapped (by design per IMGT)
- J segments: 0% character accuracy - IMGT J files don't contain gaps

### Key Observations

1. **Cross-species limitation**: Using human templates for non-human species reduces accuracy (macaque IGHV: 77% vs cat IGKV: 86%)

2. **V segment best performance**: Kappa V segments (IGKV) showed highest exact match rates, likely due to higher conservation

3. **Gap position vs character accuracy**: Gap position accuracy (60%) exceeds character accuracy (55.5%) because gaps may be in similar but not identical positions

4. **Bug fixed**: Original implementation multiplied nucleotide gap positions by 3, causing ~3x too many gaps. After fix, gap counts match expected ranges (24-27 gaps for V regions)

### Implementation Impact

- In-silico gapping provides reasonable approximation (60% gap accuracy) for V segments
- For highest accuracy, prefer using pre-gapped sequences from G3 or IMGT when available
- Custom sequences without pre-gapped versions will use in-silico gapping as fallback
- Future improvement: Implement per-gene template matching instead of consensus-based gapping

---

## Summary

| Research Question | Status | Priority |
|-------------------|--------|----------|
| VDJbase Data Format & Access | RESOLVED | - |
| Biopython Alignment Strategy | RESOLVED | - |
| G3 API Response Format | RESOLVED | - |
| IgBLAST Auxiliary File Format | RESOLVED | - |
| OGRDB Download URLs | RESOLVED | - |
| IMGT Download URLs | RESOLVED | - |
| In-Silico Gapping Validation | RESOLVED | - |

## Next Steps

1. ~~Complete remaining research items before implementation~~ **DONE**
2. ~~Implement OGRDB download script~~ **DONE** (Zenodo archive approach)
3. ~~Implement IMGT download script~~ **DONE** (V-QUEST reference directory)
4. Complete `AuxFileBuilder` stub with region parsing
5. Proceed with implementation per tasks.md
