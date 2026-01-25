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

**Phases:** 19-23 (continuing from v1.1)
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

**Note:** `--use-germlines` flag has gap — germlines adapter missing IMGT region fields

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

*Milestone v1.2 created: 2026-01-23*
*Milestone v1.2 completed: 2026-01-23*
