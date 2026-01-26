# Roadmap: Germline Database Integration

**Milestone:** v1.1 Audit
**Phases:** 13-18 (continuing from v1.0)

## Phase 13: Backend Parity Audit

**Goal:** Validate germlines backend produces identical AIRR results to G3 backend

**Status:** ✓ Complete

**Requirements:**
- AUDIT-01: Run AIRR annotation with germlines backend
- AUDIT-02: Run AIRR annotation with G3 backend
- AUDIT-03: Compare results for column-level identity
- AUDIT-04: Document discrepancies with root cause analysis

**Findings:**
- 72.19% parity achieved
- Critical issue: Missing C gene constant region data in germlines module
- 15 columns missing (all C gene related)
- CDR3/junction annotation fails for 95% of sequences

**Deliverables:**
- `audit/audit.ipynb` — Comparison notebook
- `audit/audit.md` — Detailed audit report
- `audit/20260112_HCV_DB_example.csv` — Test data

---

## Phase 14: C Region Data Integration

**Goal:** Add C gene constant region data to germlines module sources and IgBLAST databases

**Depends on:** Phase 13

**Requirements:**
- CREG-01: Update germlines sources to pull C region data from IMGT/OGRDB/VDJbase
- CREG-02: Generate IgBLAST C gene databases in germlines module
- CREG-03: Verify C gene columns present in AIRR output
- CREG-04: Re-run audit to validate parity improvement

**Success Criteria:**
1. No "C gene directory not found" warnings
2. All 129 columns present in germlines output (matching G3)
3. CDR3/junction fields populated for productive sequences
4. `complete_vdj` flag matches G3 backend
5. Parity approaches 100%

**Files to modify:**
- `src/sadie/germlines/` — Source fetching to include C genes
- `src/sadie/germlines/igblast/` — Database generation
- `src/sadie/airr/igblast/germline.py` — Path resolution for C genes

**Status:** ✓ Complete

---

## Phase 15: J Gene Matching & CDR3 Annotation Fix

**Goal:** Fix J gene matching in IgBLAST to enable CDR3/junction annotation

**Depends on:** Phase 14

**Discovery:** Phase 14 revealed that CDR3 annotation failure is a pre-existing issue unrelated to C genes:
- J genes not being matched (`j_call = NaN` for 99% of sequences)
- CDR3, junction, fwr4 all return NaN
- C gene integration successful but masked this underlying issue

**Requirements:**
- JFIX-01: Investigate IgBLAST J gene database configuration
- JFIX-02: Verify aux file format and content
- JFIX-03: Check internal_data directory structure
- JFIX-04: Debug IgBLAST execution and parameters
- JFIX-05: Re-run audit to validate CDR3 annotation

**Success Criteria:**
1. J genes matched for productive sequences (`j_call` populated)
2. CDR3/junction fields populated
3. `complete_vdj` = True for valid sequences
4. Parity with G3 backend approaches 100%

**Files to investigate:**
- `src/sadie/airr/igblast/igblast.py` — IgBLAST execution
- `src/sadie/germlines/igblast/aux_db/` — Auxiliary files
- `src/sadie/germlines/igblast/Ig/internal_data/` — Internal data structure
- `src/sadie/germlines/igblast/database/` — BLAST databases

**Status:** ✓ Complete

---

## Phase 16: Fix ndm.imgt FWR3 End Position

**Goal:** Fix the ndm.imgt file generation to use correct FWR3 end position (IMGT position 312) instead of full V gene sequence length

**Depends on:** Phase 15

**Discovery:** Phase 15 research revealed that backend parity is limited to 77.6% because:
- Column 11 in ndm.imgt uses full V gene length (e.g., 296nt) instead of FWR3 end (288nt)
- This causes FWR3 to be ~8 nucleotides too long
- CDR3 is ~8 nucleotides too short
- All junction/region boundary calculations are offset

**Requirements:**
- NDM-01: Fix `build_internal_data.py` to calculate ungapped position of IMGT position 312
- NDM-02: Regenerate ndm.imgt files for human
- NDM-03: Re-run audit to validate parity improvement

**Success Criteria:**
1. ndm.imgt column 11 matches G3 values (e.g., 288 for IGHV1-69*01)
2. FWR3/CDR3 boundaries match G3
3. Backend parity improves from 77.6% toward 95%+

**Files to modify:**
- `src/sadie/germlines/scripts/build_internal_data.py` — Fix column 11 calculation
- `src/sadie/germlines/igblast/Ig/internal_data/human/human.ndm.imgt` — Regenerate

**Status:** ✓ Complete

---

## Phase 17: Fix complete_vdj IgBLAST Quirk

**Goal:** Ensure complete_vdj=True for sequences with valid VDJ alignments

**Depends on:** Phase 16

**Status:** ✓ Complete

**Discovery:** Audit revealed 22 sequences (28%) have complete_vdj=False in germlines but True in G3, despite IDENTICAL alignment coordinates. Root cause is IgBLAST internal behavior dependent on database configuration/size, not alignment quality.

**Key Findings:**
- All position fields identical between backends (v_germline_start/end, j_germline_start/end, etc.)
- Same V gene calls, same productive status
- Issue occurs even with IMGT-only database rebuild
- G3 internal_data has combined V+D+J+C file (non-standard but works)
- Germlines uses V-only symlinks (closer to NCBI standard but triggers quirk)

**Requirements:**
- VDJ-01: Investigate post-processing solution to recalculate complete_vdj from position data ✓
- VDJ-02: OR match G3's internal_data combined file structure (skipped - post-processing chosen)
- VDJ-03: Verify complete_vdj matches G3 for all sequences ✓
- VDJ-04: Document the IgBLAST quirk in audit/igblast-quirk.md ✓

**Results:**
- Implemented AIRR-standard-based recalculation in `Airr._recalculate_complete_vdj()`
- complete_vdj differences reduced: 22 → 4
- Direction reversed: germlines now MORE correct than G3 (SADIE=True/correct, G3=False/incorrect)
- All 22 original false negatives fixed
- Pure structural parity: 98.29% (complete_vdj is allele-dependent, not structural)

**Files Modified:**
- `src/sadie/germlines/builders/j_gene_data.py` — Added J_GENE_LENGTHS dict
- `src/sadie/airr/airr.py` — Added _recalculate_complete_vdj, integrated in run_fasta and _run_scfv
- `audit/igblast-quirk.md` — Documented quirk and resolution

---

## Phase 18: Document D-region IMGT Version Variance

**Goal:** Document remaining D-region boundary differences as acceptable IMGT version variance

**Depends on:** Phase 17

**Status:** ✓ Complete

**Discovery:** Research revealed the 1.71% structural difference (108 values) is caused by IMGT database version differences, not a bug.

**Key Findings:**
- Germlines uses current IMGT GENE-DB: 40 D alleles
- G3 uses legacy snapshot: 34 D alleles
- Germlines includes 5 OR15 orphon genes not in G3
- Germlines includes 3 newer *03 alleles from BK063800
- Germlines correctly deduplicates 2 identical sequences

**Requirements:**
- DREG-01: Investigate D gene allele selection differences ✓
- DREG-02: Compare D gene database content between germlines and G3 ✓
- DREG-03: Analyze np1/np2 calculation logic ✓
- DREG-04: Document as acceptable variance ✓

**Decision:** Accept 98.29% structural parity as final. The germlines module produces MORE accurate annotations due to current IMGT data.

**Deliverables:**
- `audit/parity-notes.md` — Parity explanation documentation
- `.planning/phases/phase-18/RESEARCH.md` — Detailed research
- `.planning/phases/phase-18/SUMMARY.md` — Phase summary

---
*Created: 2026-01-22*
*Completed: 2026-01-23*

---

## Milestone: v1.2 Reference Module Unification

**Phases:** 19-24 (continuing from v1.1)
**Goal:** Enable reference.yml to select alleles from all germline sources (imgt, ogrdb, vdjbase, custom), using germlines module as data provider instead of G3 API. Full workflow: build CLI + runtime usage of prebuilt databases.

---

## Phase 19: Source Validation

**Goal:** Expand source validation to accept all germline database providers

**Depends on:** Phase 18 (v1.1 complete)

**Requirements:**
- SRC-01: Expand VALID_SOURCES to include `ogrdb`, `vdjbase`
- SRC-02: Validate source exists in germlines before processing

**Success Criteria:**
1. `models.py` accepts `source: ogrdb` and `source: vdjbase` without validation errors
2. Error message when source not available in germlines for requested species
3. Existing `imgt` and `custom` sources continue to work unchanged
4. Unit tests cover all four source types

**Files to modify:**
- `src/sadie/reference/models.py` — Expand VALID_SOURCES list in `check_source` validators

---

## Phase 20: Integration Foundation

**Goal:** Route reference.yml processing through germlines module with explicit source selection

**Depends on:** Phase 19

**Status:** ✓ Complete

**Requirements:**
- INT-01: Add `use_germlines=True` parameter to `References.from_yaml()`
- INT-02: Route source selection through GermlineManager (explicit source, no priority)
- INT-03: Generate synthetic `_id` field in adapter

**Success Criteria:**
1. ✓ `References.from_yaml(use_germlines=True)` loads genes from germlines module
2. ✓ Source field from YAML explicitly passed to GermlineManager (not using priority fallback)
3. ✓ All returned gene dicts contain `_id` field (hash of `source:species:gene`)
4. ✓ Downstream code using `_id` for deduplication/indexing works correctly
5. ✓ G3 API path still works with `use_germlines=False` for backwards compatibility

**Files modified:**
- `src/sadie/reference/reference.py` — Add `use_germlines` param, explicit source routing
- `src/sadie/germlines/g3_adapter.py` — Add `_id` field generation using SHA-256 hash

---

## Phase 21: Build CLI

**Goal:** Add CLI command to build IgBLAST database from reference.yml

**Depends on:** Phase 20

**Status:** ✓ Complete

**Requirements:**
- CLI-01: Add `sadie reference build <yaml> --output <path>` command
- CLI-02: Build generates complete IgBLAST database structure
- CLI-03: Progress output during build

**Success Criteria:**
1. ✓ `sadie reference build reference.yml --output ./db` creates database directory
2. ✓ Output contains: `Ig/blastdb/`, `Ig/internal_data/`, `aux_db/`, `.references_dataframe.csv.gz`
3. ✓ Progress output shows: "Loading YAML...", "Fetching genes...", "Building databases...", "Complete"
4. ✓ Exit code 0 on success, non-zero with error message on failure
5. ✓ Resulting database structure identical to `References.make_airr_database()` output

**Files modified:**
- `src/sadie/app.py` — Added `@reference.command("build")` with Click decorators

**Note:** `--use-germlines` gap closed in Phase 24

---

## Phase 22: Runtime Usage

**Goal:** Enable Airr to use prebuilt databases directly, bypassing runtime gene lookup

**Depends on:** Phase 21

**Status:** ✓ Complete

**Requirements:**
- RUN-01: Add `Airr(database=<path>)` parameter to use prebuilt database
- RUN-02: Skip germlines/G3 lookup when using prebuilt database
- RUN-03: Validate database structure on load

**Success Criteria:**
1. ✓ `Airr(database="./db")` uses prebuilt database instead of default
2. ✓ No network calls or germlines lookups when database path provided
3. ✓ Clear error if database path missing required structure (blastdb, internal_data, aux_db)
4. ✓ Annotation results identical whether using prebuilt or runtime-built database
5. ⚠ Performance not measured

**Files modified:**
- `src/sadie/airr/airr.py` — Add `database` parameter, pass to recursive calls
- `src/sadie/airr/igblast/germline.py` — Add `validate_prebuilt_database()`, `prebuilt` param

---

## Phase 23: Documentation

**Goal:** Document multi-source reference.yml usage and build workflow

**Depends on:** Phase 22

**Status:** ✓ Complete

**Requirements:**
- DOC-01: Create reference-sample.yml (mouse=imgt, human=ogrdb, macaque=vdjbase)
- DOC-02: Document build → use workflow

**Success Criteria:**
1. ✓ `examples/reference-sample.yml` demonstrates multi-source configuration
2. ✓ Sample includes: mouse V/D/J from IMGT, human V/D/J from OGRDB, macaque subset from VDJbase
3. ✓ docs/reference-workflow.md explains: write YAML → build database → use in Airr
4. ✓ Code examples show complete workflow from YAML to annotation

**Files created:**
- `examples/reference-sample.yml` — Multi-source reference configuration
- `docs/reference-workflow.md` — Build and usage documentation
- `mkdocs.yml` — Added navigation entry

---

## Phase 24: v1.2 Gap Closure (Validation + IMGT Regions)

**Goal:** Close remaining gaps so `use_germlines` flows work end-to-end and milestone can be verified

**Depends on:** Phase 23

**Status:** ✓ Complete

**Requirements:**
- SRC-02: Validate source exists in germlines before processing

**Gaps closed:**
- ✓ SRC-02: Validate source exists in germlines before processing
- ✓ `--use-germlines` build failure (missing IMGT region fields)
- ✓ Missing phase verification artifacts (19-23)

**Success Criteria:**
1. ✓ `sadie reference build ... --use-germlines` completes with full database structure
2. ✓ `References.from_yaml(use_germlines=True)` builds internal data without missing IMGT fields
3. ✓ Provider/species validation fails fast with clear error
4. ✓ Phase 19 SUMMARY and phases 19-23 VERIFICATION.md files exist

---

*Milestone v1.2 created: 2026-01-23*
*Milestone v1.2 completed: 2026-01-25*

---

## Milestone: v1.3 Test Infrastructure & Species Expansion

**Phases:** 25-28
**Goal:** Fix skipped tests by adding macaque germlines, airr package dependency, removing deprecated G3 tests, and fix germline priority order

---

## Phase 25: Macaque Germlines Integration

**Goal:** Build macaque germline databases to enable skipped macaque tests

**Depends on:** Phase 24 (v1.2 complete)

**Status:** ✓ Complete

**Context:** 6 tests currently skipped with "macaque germlines not available":
- `test_five_and_three_prime_extension`
- `test_hard_igl_seqs`
- `test_hard_igl_seqs_linked`
- `test_airr_constant_region_macaque`
- `test_run_five_prime_buffer`
- `test_run_three_prime_buffer`

**Requirements:**
- MAC-01: Build macaque IgBLAST databases in germlines module
- MAC-02: Generate macaque internal_data and aux files
- MAC-03: Verify macaque AIRR annotation works
- MAC-04: Remove skip markers from macaque tests

**Success Criteria:**
1. `GermlineData("macaque")` resolves without error
2. All 6 macaque tests pass (not skipped)
3. Macaque annotation produces valid AIRR output

**Files to modify:**
- `src/sadie/germlines/igblast/Ig/internal_data/macaque/` — Create database
- `src/sadie/germlines/igblast/aux_db/macaque_gl.aux` — Create aux file
- `tests/unit/airr/test_airr.py` — Remove `@skip_no_macaque` markers
- `tests/unit/airr/test_methods.py` — Remove `@skip_no_macaque` markers

---

## Phase 26: Add AIRR Package Dependency

**Goal:** Add airr package to dependencies so AIRR validation test runs

**Depends on:** Phase 25

**Status:** ✓ Complete

**Context:** 1 test skipped with "airr package not installed":
- `test_write_and_check_airr`

**Requirements:**
- AIRR-01: Add `airr` package to pyproject.toml dependencies
- AIRR-02: Verify airr package installs correctly
- AIRR-03: Remove importorskip from test

**Success Criteria:**
1. ✓ `pip install sadie` includes airr package
2. ✓ `test_write_and_check_airr` passes (not skipped)
3. ✓ AIRR table validation works with official airr package

**Files modified:**
- `pyproject.toml` — Moved airr from dev to main dependencies
- `tests/unit/airr/test_airr.py` — Added `import airr`, removed `pytest.importorskip("airr")`

---

## Phase 27: Remove Deprecated G3 Tests

**Goal:** Remove or migrate tests that depend on deprecated G3 API

**Depends on:** Phase 26

**Status:** ✓ Complete

**Context:** 2 tests skipped with "G3 API deprecated, will be removed after 2026-06-01":
- `test_v_gene_dir_attribute_exists`
- `test_aux_path_attribute_exists`

**Requirements:**
- G3-01: Review what these tests are validating
- G3-02: Determine if equivalent germlines module tests exist
- G3-03: Either migrate tests to germlines or remove if redundant
- G3-04: Update deprecation timeline documentation

**Success Criteria:**
1. ✓ No tests reference deprecated G3 API
2. ✓ All test functionality covered by germlines module tests
3. ✓ Clear documentation of G3 deprecation timeline

**Files modified:**
- `tests/unit/germlines/test_germline_data_legacy.py` — Removed `TestGermlineDataPaths` class
- `docs/G3-Deprecation.md` — Created deprecation documentation (new file)

---

## Phase 28: Fix Germline Priority Order

**Goal:** Reorder germline provider priority to ['vdjbase', 'ogrdb', 'imgt', 'custom'] for optimal data quality

**Depends on:** Phase 27

**Status:** ✓ Complete

**Context:** Current priority order doesn't reflect data quality hierarchy:
- VDJbase: Best for human and macaque (curated, validated alleles)
- OGRDB: Good fallback for human, excellent for mouse (community-curated)
- IMGT: Good for species diversity (comprehensive reference)
- Custom: Fill gaps from internal lab data

**Requirements:**
- PRIO-01: Update default provider priority in GermlineManager
- PRIO-02: Document priority rationale in code comments
- PRIO-03: Verify priority order used in fallback resolution
- PRIO-04: Test priority order with multi-source queries

**Success Criteria:**
1. ✓ `GermlineManager` default priority is `['vdjbase', 'ogrdb', 'imgt', 'custom']`
2. ✓ Human/macaque queries prefer VDJbase alleles when available
3. ✓ Mouse queries get OGRDB alleles (VDJbase has limited mouse data)
4. ✓ Fallback chain works correctly when preferred source lacks data

**Files modified:**
- `src/sadie/germlines/manager.py` — Updated DEFAULT_PROVIDERS order + docstrings
- `src/sadie/germlines/pipeline.py` — Updated hardcoded source iteration
- `tests/unit/germlines/test_compliance.py` — Updated test assertions
- `tests/unit/germlines/test_manager.py` — Updated test assertions
- `src/sadie/germlines/__init__.py` — Updated module docstrings
- Documentation files (7 files) — Updated priority references

---

*Milestone v1.3 created: 2026-01-25*
*Milestone v1.3 completed: 2026-01-25*

---

## Phase 29: Add Germline Source Tracking to AIRR Output

**Goal:** Add columns to AIRR output showing which germline database each gene call came from

**Depends on:** Phase 28

**Status:** Not started

**Context:** Users cannot currently tell which germline source (imgt, vdjbase, ogrdb, custom) each matched gene came from in the AIRR output. The `source` field exists in `GermlineGene` objects at the GermlineManager level, but is not propagated to the final annotation results.

**Use Case:**
- Researchers need to cite the correct germline database for publications
- Quality control: verify expected sources are being used
- Debugging: understand why specific alleles were matched
- Reproducibility: document exact data sources used

**Requirements:**
- SRC-01: Add `v_call_source`, `d_call_source`, `j_call_source`, `c_call_source` columns to AIRR output
- SRC-02: Populate source columns during IgBLAST result parsing
- SRC-03: Handle cases where source lookup fails (return "unknown")
- SRC-04: Include sources in LinkedAirrTable with appropriate suffixes

**Success Criteria:**
1. AIRR output contains `v_call_source`, `d_call_source`, `j_call_source`, `c_call_source` columns
2. Source values match germline provider names: "imgt", "vdjbase", "ogrdb", "custom"
3. Source columns populated for all matched genes
4. Unmatched genes (NaN calls) have NaN sources
5. LinkedAirrTable has `_heavy`/`_light` suffixed source columns

**Implementation Notes:**
- Create gene name → source lookup table during GermlineData initialization
- Cache lookup table for performance (avoid repeated source lookups)
- Consider adding source info to aux file or separate metadata file

**Files to modify:**
- `src/sadie/airr/igblast/germline.py` — Build gene→source lookup table
- `src/sadie/airr/igblast/igblast.py` — Add source columns during result parsing
- `src/sadie/airr/airrtable/airrtable.py` — Handle source columns in table operations
- `tests/unit/airr/test_airr.py` — Add tests for source columns

---

*Milestone v1.3 created: 2026-01-25*
