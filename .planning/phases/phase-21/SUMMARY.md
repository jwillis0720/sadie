# Phase 21: Build CLI — SUMMARY

**Status:** Complete
**Date:** 2026-01-23

## Tasks Executed

### Task 1: Add `build` command to reference group ✓
**File:** `src/sadie/app.py`
**Commit:** `5469a16a`

- Added `@reference.command("build")` with Click decorators
- Positional `yaml_path` argument
- Required `--output/-o` option
- Optional `--use-germlines` flag
- Progress output: "Loading YAML...", "Fetching genes...", "Building databases...", "Complete"
- Error handling with exit code 1 on failure

### Task 2: Phase 20 dependency ✓
Already complete - `use_germlines` parameter added in Phase 20.

### Task 3: Verify output structure ✓
Tested with `short_reference.yml` fixture:
```
/tmp/test_build_db/
├── .references_dataframe.csv.gz
├── aux_db/imgt/
└── Ig/
    ├── blastdb/
    └── internal_data/
```

## Requirements Satisfied

| Requirement | Description | Status |
|-------------|-------------|--------|
| CLI-01 | Add `sadie reference build <yaml> --output <path>` command | ✓ |
| CLI-02 | Build generates complete IgBLAST database structure | ✓ |
| CLI-03 | Progress output during build | ✓ |

## Success Criteria Verification

1. ✓ `sadie reference build reference.yml --output ./db` creates database directory
2. ✓ Output contains: `Ig/blastdb/`, `Ig/internal_data/`, `aux_db/`, `.references_dataframe.csv.gz`
3. ✓ Progress output shows: "Loading YAML...", "Fetching genes...", "Building databases...", "Complete"
4. ✓ Exit code 0 on success, non-zero with error message on failure
5. ✓ Database structure identical to `References.make_airr_database()` output

## Known Limitation

The `--use-germlines` flag has a gap: the germlines adapter doesn't provide all IMGT region fields (`imgt.fwr1_start`, `imgt.cdr1_end`, etc.) that `make_airr_database()` requires for creating internal annotation files.

**Workaround:** Use without `--use-germlines` (G3 API path) for now.

**Future work:** Enhance germlines adapter to include all IMGT region data.

## Test Results

```
sadie reference build --help  # Works ✓
sadie reference build short_reference.yml --output /tmp/db  # Works ✓ (G3 path)
sadie reference build ref.yml -o /tmp/db --use-germlines  # Gap: missing IMGT fields
```
