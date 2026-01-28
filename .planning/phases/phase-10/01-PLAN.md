# Phase 10 Plan: Species Expansion

## Goal

Populate IgBLAST databases for all 33 IMGT-supported species to enable multi-species AIRR annotation and renumbering analysis.

## Context

**Current State**:
- Human IMGT data downloaded and IgBLAST databases built
- Mouse IMGT data downloaded but NO IgBLAST databases built
- 31 other species: no IMGT data, no databases

**Target State**:
- All 33 species have IMGT germline data downloaded
- All species with available data have IgBLAST BLAST databases built
- Auxiliary files (*.aux) exist for J gene CDR3 start positions
- organism.yaml configured for all species
- Multi-species annotation verified working

**External Dependency**: BLAST+ tools (makeblastdb) must be installed

## Species Groups (33 total)

| Group | Species | Count |
|-------|---------|-------|
| Primates | human, rhesus_macaque, cynomolgus, gorilla, chimpanzee, orangutan_sumatran, orangutan_bornean, lemur, owl_monkey | 9 |
| Rodents | mouse, mouse_c57bl6j, rat, naked_mole_rat | 4 |
| Carnivores | dog, cat, ferret, mink | 4 |
| Ungulates | rabbit, pig, cow, sheep, goat, horse, alpaca, camel | 8 |
| Birds | chicken | 1 |
| Fish | zebrafish, atlantic_salmon, rainbow_trout, atlantic_cod, channel_catfish | 5 |
| Marine Mammals | dolphin | 1 |
| Monotremes | platypus | 1 |

## Success Criteria

1. ✅ IMGT data downloaded for all 33 species in SPECIES_MAP
2. ✅ IgBLAST BLAST databases built for each species with available data
3. ✅ Auxiliary files (*.aux) created with J gene CDR3 positions
4. ✅ organism.yaml updated with all species configurations
5. ✅ Multi-species AIRR annotation tests pass (mouse, rhesus, chicken/zebrafish)
6. ✅ Multi-species renumbering tests pass (mouse, rabbit)

---

## Task Breakdown

### Wave 1: Infrastructure & Tooling (Blocking)

These tasks create the tools needed for automated species processing.

#### T073: Create auxiliary file generator script
**File**: `src/sadie/germlines/scripts/build_aux_files.py`

**Purpose**: Generate *.aux files for IgBLAST J gene CDR3 start position detection

**Implementation**:
```python
# Aux file format (tab-separated):
# GENE_ID  1  CDR3_START_POS  (additional CDR/FR positions)
#
# The script must:
# 1. Read J gene sequences from normalized data
# 2. Identify CDR3 start position from IMGT numbering (position 105)
# 3. Write aux file in IgBLAST expected format
```

**Acceptance**: Script generates valid aux files; human_gl.aux matches existing format

---

#### T074: Create organism.yaml update script
**File**: `src/sadie/germlines/scripts/update_organism_yaml.py`

**Purpose**: Automatically add species entries to organism.yaml after database building

**Implementation**:
```python
# For each species with built databases:
# organisms:
#   {species}:
#     aux_file: ../aux_db/{species}_gl.aux
#     database_path: ../database/{species}
#     segments:
#     - V
#     - D
#     - J
```

**Acceptance**: Script correctly updates organism.yaml; YAML remains valid

---

### Wave 2: Data Acquisition (Parallel by Group)

Download IMGT data for all species. Each group can run in parallel.

#### T075: [P] Download IMGT data for primate species (8 new)
**Command**:
```bash
python -m sadie.germlines.scripts.download_imgt \
  --species rhesus_macaque cynomolgus gorilla chimpanzee \
           orangutan_sumatran orangutan_bornean lemur owl_monkey
```
**Output**: `src/sadie/germlines/sources/imgt/{species}/` directories

---

#### T076: [P] Download IMGT data for rodent species (3 new)
**Command**:
```bash
python -m sadie.germlines.scripts.download_imgt \
  --species mouse_c57bl6j rat naked_mole_rat
```
**Output**: `src/sadie/germlines/sources/imgt/{species}/` directories

---

#### T077: [P] Download IMGT data for carnivore species (4)
**Command**:
```bash
python -m sadie.germlines.scripts.download_imgt \
  --species dog cat ferret mink
```
**Output**: `src/sadie/germlines/sources/imgt/{species}/` directories

---

#### T078: [P] Download IMGT data for ungulate species (8)
**Command**:
```bash
python -m sadie.germlines.scripts.download_imgt \
  --species rabbit pig cow sheep goat horse alpaca camel
```
**Output**: `src/sadie/germlines/sources/imgt/{species}/` directories

---

#### T079: [P] Download IMGT data for birds, fish, marine mammals, monotremes (8)
**Command**:
```bash
python -m sadie.germlines.scripts.download_imgt \
  --species chicken zebrafish atlantic_salmon rainbow_trout \
           atlantic_cod channel_catfish dolphin platypus
```
**Output**: `src/sadie/germlines/sources/imgt/{species}/` directories

---

### Wave 3: Normalization Pipeline (Sequential per Group)

Run the normalization pipeline to process downloaded IMGT data.

#### T080: Run normalization pipeline for all downloaded species
**Implementation**:
```python
from sadie.germlines.pipeline import normalize_species

for species in ALL_SPECIES:
    normalize_species(species, providers=["imgt"])
```

**Output**: `src/sadie/germlines/normalized/{species}/` directories

**Note**: This creates the ungapped FASTA files needed for BLAST database building

---

### Wave 4: Database Building (Parallel per Species)

Build IgBLAST databases for all species with normalized data.

#### T081: [P] Build BLAST databases for primate species (9)
**Implementation**:
```python
from sadie.germlines.builders.blast import BlastDBBuilder
from pathlib import Path

builder = BlastDBBuilder()
species_list = ["human", "rhesus_macaque", "cynomolgus", "gorilla",
                "chimpanzee", "orangutan_sumatran", "orangutan_bornean",
                "lemur", "owl_monkey"]

for species in species_list:
    source_dir = Path(f"src/sadie/germlines/normalized/{species}/ungapped")
    output_dir = Path(f"src/sadie/germlines/igblast/database/{species}")
    if source_dir.exists():
        builder.build_for_species(species, source_dir, output_dir)
```

**Output**: BLAST database files (*.nhr, *.nin, *.nsq, *.ndb, *.not, *.ntf, *.nto)

---

#### T082: [P] Build BLAST databases for rodent species (4)
**Species**: mouse, mouse_c57bl6j, rat, naked_mole_rat
**Output**: `src/sadie/germlines/igblast/database/{species}/` directories

---

#### T083: [P] Build BLAST databases for carnivore species (4)
**Species**: dog, cat, ferret, mink
**Output**: `src/sadie/germlines/igblast/database/{species}/` directories

---

#### T084: [P] Build BLAST databases for ungulate species (8)
**Species**: rabbit, pig, cow, sheep, goat, horse, alpaca, camel
**Output**: `src/sadie/germlines/igblast/database/{species}/` directories

---

#### T085: [P] Build BLAST databases for other species (8)
**Species**: chicken, zebrafish, atlantic_salmon, rainbow_trout, atlantic_cod, channel_catfish, dolphin, platypus
**Output**: `src/sadie/germlines/igblast/database/{species}/` directories

---

### Wave 5: Auxiliary Files Generation

#### T086: Generate auxiliary files for all species
**Implementation**:
```bash
python -m sadie.germlines.scripts.build_aux_files --all-species
```

**Output**: `src/sadie/germlines/igblast/aux_db/{species}_gl.aux` files

---

### Wave 6: Configuration Update

#### T087: Update organism.yaml with all species
**Implementation**:
```bash
python -m sadie.germlines.scripts.update_organism_yaml --all-species
```

**Output**: Updated `src/sadie/germlines/igblast/internal_data/organism.yaml`

---

### Wave 7: Verification

#### T088: Verify BLAST database integrity for all species
**Implementation**:
```python
def verify_blast_databases():
    required_extensions = [".nhr", ".nin", ".nsq"]
    for species in ALL_SPECIES:
        db_dir = Path(f"igblast/database/{species}")
        for segment in ["V", "D", "J"]:
            db_path = db_dir / f"{species}_{segment}"
            for ext in required_extensions:
                assert (db_path.with_suffix(ext)).exists(), f"Missing {db_path}{ext}"
```

**Test File**: `tests/unit/germlines/test_multi_species.py::test_blast_database_integrity`

---

#### T089: Verify organism.yaml contains all species
**Test**: Assert all 33 species have entries in organism.yaml with correct paths

---

### Wave 8: Integration Testing

#### T090: Test AIRR annotation with mouse species
**File**: `tests/unit/germlines/test_multi_species.py`
```python
def test_airr_annotation_mouse():
    """Verify AIRR annotation works with mouse germlines."""
    result = run_airr(sequence, species="mouse", germline_backend="germlines")
    assert result.v_call is not None
    assert "IGHV" in result.v_call or "IGKV" in result.v_call or "IGLV" in result.v_call
```

---

#### T091: Test AIRR annotation with non-human primate (rhesus_macaque)
**File**: `tests/unit/germlines/test_multi_species.py`
```python
def test_airr_annotation_rhesus():
    """Verify AIRR annotation works with rhesus macaque germlines."""
    result = run_airr(sequence, species="rhesus_macaque", germline_backend="germlines")
    assert result.v_call is not None
```

---

#### T092: Test AIRR annotation with non-mammalian species (chicken or zebrafish)
**File**: `tests/unit/germlines/test_multi_species.py`
```python
def test_airr_annotation_chicken():
    """Verify AIRR annotation works with chicken germlines."""
    # Note: May need to adjust if chicken has limited IMGT data
    result = run_airr(sequence, species="chicken", germline_backend="germlines")
    # Assert appropriate behavior based on available data
```

---

#### T093: Test renumbering HMM generation for mouse
**File**: `tests/unit/germlines/test_multi_species.py`
```python
def test_renumbering_mouse():
    """Verify HMM generation works for mouse germlines."""
    from sadie.germlines.renumbering_integration import LocalHMMBuilder
    builder = LocalHMMBuilder(species="mouse")
    hmm = builder.get_hmm("H")
    assert hmm is not None
```

---

#### T094: Test renumbering HMM generation for rabbit
**File**: `tests/unit/germlines/test_multi_species.py`
```python
def test_renumbering_rabbit():
    """Verify HMM generation works for rabbit germlines."""
    builder = LocalHMMBuilder(species="rabbit")
    hmm = builder.get_hmm("H")
    assert hmm is not None
```

---

#### T095: Create comprehensive multi-species integration test suite
**File**: `tests/unit/germlines/test_multi_species.py`

**Tests**:
- `test_species_data_availability`: Check which species have data
- `test_blast_database_integrity`: Verify BLAST files exist
- `test_aux_files_exist`: Verify auxiliary files for all species
- `test_organism_yaml_complete`: Verify config completeness
- `test_airr_annotation_mouse`: Mouse annotation test
- `test_airr_annotation_rhesus`: Rhesus annotation test
- `test_airr_annotation_chicken`: Non-mammalian test
- `test_renumbering_mouse`: Mouse renumbering test
- `test_renumbering_rabbit`: Rabbit renumbering test
- `test_multi_species_batch`: Batch processing across species

---

## Dependencies

```
Wave 1: Infrastructure (T073, T074)
         ↓
Wave 2: Data Acquisition (T075-T079) [PARALLEL]
         ↓
Wave 3: Normalization (T080)
         ↓
Wave 4: Database Building (T081-T085) [PARALLEL]
         ↓
Wave 5: Aux Files (T086)
         ↓
Wave 6: Config Update (T087)
         ↓
Wave 7: Verification (T088-T089) [PARALLEL]
         ↓
Wave 8: Integration Tests (T090-T095)
```

## Parallel Execution Opportunities

| Wave | Tasks | Parallel? | Notes |
|------|-------|-----------|-------|
| 1 | T073, T074 | Yes | Independent scripts |
| 2 | T075-T079 | Yes | Different species groups |
| 3 | T080 | No | Sequential per species |
| 4 | T081-T085 | Yes | Different species groups |
| 5 | T086 | No | Single script |
| 6 | T087 | No | Single script |
| 7 | T088-T089 | Yes | Independent checks |
| 8 | T090-T095 | Partial | Some tests can parallelize |

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Some species lack IMGT data | Medium | Low | Skip unavailable species gracefully |
| makeblastdb not installed | Low | High | Check at start; document dependency |
| Large download size | Medium | Low | Download in batches; checkpoint progress |
| Some species lack J genes | Medium | Medium | Handle missing segments in aux generation |
| Test sequences unavailable | Medium | Medium | Use V gene fragments for testing |

## Verification Checklist

Before marking Phase 10 complete:

- [ ] All IMGT downloads completed (check sources/imgt/{species}/)
- [ ] All normalizations completed (check normalized/{species}/)
- [ ] All BLAST databases built (check igblast/database/{species}/)
- [ ] All aux files generated (check igblast/aux_db/{species}_gl.aux)
- [ ] organism.yaml updated with all species
- [ ] `pytest tests/unit/germlines/test_multi_species.py` passes
- [ ] Mouse AIRR annotation verified working
- [ ] At least one non-human primate annotation verified
- [ ] At least one non-mammalian species annotation verified
- [ ] Mouse renumbering verified working
- [ ] Rabbit renumbering verified working

## Task Summary

| ID | Description | Wave | Parallel |
|----|-------------|------|----------|
| T073 | Create auxiliary file generator script | 1 | Yes |
| T074 | Create organism.yaml update script | 1 | Yes |
| T075 | Download IMGT data for primates (8 new) | 2 | Yes |
| T076 | Download IMGT data for rodents (3 new) | 2 | Yes |
| T077 | Download IMGT data for carnivores (4) | 2 | Yes |
| T078 | Download IMGT data for ungulates (8) | 2 | Yes |
| T079 | Download IMGT data for birds/fish/other (8) | 2 | Yes |
| T080 | Run normalization pipeline for all species | 3 | No |
| T081 | Build BLAST databases for primates (9) | 4 | Yes |
| T082 | Build BLAST databases for rodents (4) | 4 | Yes |
| T083 | Build BLAST databases for carnivores (4) | 4 | Yes |
| T084 | Build BLAST databases for ungulates (8) | 4 | Yes |
| T085 | Build BLAST databases for other species (8) | 4 | Yes |
| T086 | Generate auxiliary files for all species | 5 | No |
| T087 | Update organism.yaml with all species | 6 | No |
| T088 | Verify BLAST database integrity | 7 | Yes |
| T089 | Verify organism.yaml completeness | 7 | Yes |
| T090 | Test AIRR annotation with mouse | 8 | Partial |
| T091 | Test AIRR annotation with rhesus | 8 | Partial |
| T092 | Test AIRR with chicken/zebrafish | 8 | Partial |
| T093 | Test renumbering for mouse | 8 | Partial |
| T094 | Test renumbering for rabbit | 8 | Partial |
| T095 | Create multi-species test suite | 8 | No |

**Total Tasks**: 23 (T073-T095)

---

## Mapping to Original ROADMAP Tasks

| ROADMAP Task | Phase 10 Tasks | Notes |
|--------------|----------------|-------|
| T061 (Download IMGT) | T075-T079 | Split by species groups |
| T062 (Aux generator) | T073 | Infrastructure task |
| T063 (Build databases) | T080-T085 | Split by species groups |
| T064 (Generate aux files) | T086 | Execution task |
| T065 (Update organism.yaml) | T074, T087 | Script + execution |
| T066 (Verify databases) | T088-T089 | Verification tasks |
| T067 (Test mouse AIRR) | T090 | Direct mapping |
| T068 (Test rhesus AIRR) | T091 | Direct mapping |
| T069 (Test non-mammalian) | T092 | Direct mapping |
| T070 (Test mouse renumber) | T093 | Direct mapping |
| T071 (Test rabbit renumber) | T094 | Direct mapping |
| T072 (Multi-species tests) | T095 | Direct mapping |

---

*Created: 2026-01-21*
