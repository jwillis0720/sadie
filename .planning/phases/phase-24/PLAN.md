# Phase 24: v1.2 Gap Closure (Validation + IMGT Regions)

## Goal
Close the remaining v1.2 gaps so `use_germlines` flows work end-to-end and the milestone can be verified.

## Context Analysis

### Gaps from v1.2 Audit
1. **SRC-02 missing**: No species-level source validation when `use_germlines=True`.
2. **IMGT region fields missing**: `GermlineToG3Adapter` does not populate IMGT region positions, so `make_airr_database()` fails.
3. **VDJbase ungapped sequences**: No IMGT gapped sequences, blocking region position derivation.
4. **Verification artifacts missing**: Phase 19 lacks SUMMARY; phases 19-23 lack VERIFICATION.

### Target Files
- `src/sadie/reference/reference.py`
- `src/sadie/germlines/manager.py`
- `src/sadie/germlines/g3_adapter.py`
- `src/sadie/germlines/providers/vdjbase.py`
- `src/sadie/germlines/builders/gapper.py` (reuse)
- `tests/unit/reference/test_reference.py`
- `tests/unit/germlines/test_reference_integration.py`
- `tests/unit/germlines/test_vdjbase_provider.py`
- `docs/reference-workflow.md`
- `.planning/phases/phase-19/SUMMARY.md`
- `.planning/phases/phase-19/VERIFICATION.md`
- `.planning/phases/phase-20/VERIFICATION.md`
- `.planning/phases/phase-21/VERIFICATION.md`
- `.planning/phases/phase-22/VERIFICATION.md`
- `.planning/phases/phase-23/VERIFICATION.md`

---

## Task 1: Implement SRC-02 (Source Availability Validation)

**Requirement:** SRC-02 — Validate source exists in germlines before processing

**Files:** `src/sadie/reference/reference.py`, `src/sadie/germlines/manager.py`

### Changes
1. Add a helper to validate provider availability for a species:
   - Option A: new `GermlineManager.get_provider(name)` and `validate_species(name, species)`.
   - Option B: new `Reference._validate_source_species(source, species)` using provider `is_available()`.
2. Call validation in `add_gene()` and `add_genes()` (germlines path) before fetching genes.
3. Error message should include provider, species, and available species list.

### Verification
- Unit test that `Reference(use_germlines=True)` raises a clear error when a source has no data for a species.

---

## Task 2: Add IMGT Region Positions for V Genes

**Goal:** Populate `imgt.fwr1_start`...`imgt.fwr3_end` for V genes so internal data can be built.

**Files:** `src/sadie/germlines/g3_adapter.py`, new helper in `src/sadie/germlines/builders/`

### Changes
1. Extract IMGT position logic from `germlines/scripts/build_internal_data.py` into a reusable helper:
   - Compute ungapped FR/CDR positions from IMGT-gapped sequence.
   - Convert to **0-based** positions to match G3 output.
2. In `GermlineToG3Adapter.to_g3_format()`:
   - If `gene.segment == "V"` and `gene.sequence_gapped` exists, derive and populate:
     - `imgt.fwr1_start`, `imgt.fwr1_end`
     - `imgt.cdr1_start`, `imgt.cdr1_end`
     - `imgt.fwr2_start`, `imgt.fwr2_end`
     - `imgt.cdr2_start`, `imgt.cdr2_end`
     - `imgt.fwr3_start`, `imgt.fwr3_end`
   - Optionally fill region sequences (`imgt.fwr1`, `imgt.cdr1`, ...) for parity.

### Verification
- Unit test that IMGT positions exist for an IMGT V gene and are integers.

---

## Task 3: Gap VDJbase V/J Sequences When Needed

**Goal:** Ensure VDJbase genes can derive IMGT positions.

**Files:** `src/sadie/germlines/providers/vdjbase.py`

### Changes
1. Use `GapperService` (as in `CustomProvider`) when `sequence_gapped` is missing for V/J genes.
2. Use IMGT templates if available; fall back to human templates when species-specific templates are missing.
3. Preserve existing behavior for D genes (no gapping).

### Verification
- Unit test that VDJbase V genes gain `sequence_gapped` when templates exist.

---

## Task 4: Fail Fast on Missing IMGT Fields During Build

**Goal:** Prevent silent partial builds when IMGT positions are missing.

**Files:** `src/sadie/reference/reference.py`, `src/sadie/app.py`

### Changes
1. In `_make_internal_annotaion_file()`, detect missing IMGT positions in V genes and raise `ValueError` with a list of offending genes.
2. Update CLI `sadie reference build` to surface the error clearly (no partial output).

### Verification
- Integration test: `References.from_yaml(..., use_germlines=True)` builds a database for `short_reference.yml`.
- Negative test: missing IMGT fields triggers a clear error.

---

## Task 5: Documentation + Verification Artifacts

**Files:** planning docs and `docs/reference-workflow.md`

### Changes
1. Create `.planning/phases/phase-19/SUMMARY.md` (SRC-01 complete, SRC-02 moved to Phase 24).
2. Create VERIFICATION.md files for phases 19-23 using evidence/tests.
3. Update `docs/reference-workflow.md` to remove `--use-germlines` gap note after fix.

---

## Execution Order

1. Task 1 — Source validation (SRC-02)
2. Task 2 — IMGT region positions for V genes
3. Task 3 — VDJbase gapping
4. Task 4 — Fail-fast on missing IMGT fields
5. Task 5 — Documentation and verification

---

## Test Strategy

- `pytest tests/unit/reference/test_reference.py -k "germlines or source"`
- `pytest tests/unit/germlines/test_reference_integration.py -k "adapter"`
- `pytest tests/unit/germlines/test_vdjbase_provider.py`
- `pytest tests/integration/reference/test_reference_integration.py -k "reference build"`

---

## Success Criteria

- `sadie reference build <yaml> --output <path> --use-germlines` completes with full database structure.
- `References.from_yaml(use_germlines=True)` builds internal data without missing IMGT positions.
- SRC-02 errors are explicit with provider/species availability details.
- Phase 19 SUMMARY and phases 19-23 VERIFICATION.md files exist.
