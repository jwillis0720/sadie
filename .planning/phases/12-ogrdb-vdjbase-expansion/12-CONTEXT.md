# Phase 12: Provider Auto-Population — Context

## Scope

Make all germline providers (IMGT, OGRDB, VDJbase) auto-populate their data programmatically via CLI command.

## Decisions

### 1. Trigger Mechanism

**Decision**: Explicit CLI command (same pattern as `sadie make`)

- Add `sadie germlines populate` CLI command
- Supports `--provider` flag to target specific providers (imgt, ogrdb, vdjbase, all)
- Supports `--species` flag for specific species (optional)
- No automatic download on import or install
- User explicitly triggers population

**Example usage**:
```bash
sadie germlines populate --provider imgt
sadie germlines populate --provider ogrdb --species human mouse
sadie germlines populate --provider all
```

### 2. Species Scope

**Decision**: All species available from each provider's API

- IMGT: Download all species in SPECIES_MAP (33 mapped species)
- OGRDB: Query API for all available species and download all
- VDJbase: Query API for all available species and download all
- No filtering or subsetting by default
- `--species` flag allows user to override if needed

### 3. Error Handling

**Decision**: Fail fast

- Stop on first error
- Clear error message indicating which species/provider failed
- No silent skipping of failed species
- User can re-run after fixing the issue
- Idempotent: re-running after failure continues from where it left off (via checkpoint)

### 4. Update Strategy

**Decision**: Version-check and update if newer

- Check provider metadata/version before downloading
- Skip download if local version matches or is newer
- `--force` flag to override and re-download regardless
- Store version metadata in provider source directories

**Version tracking**:
- IMGT: Track release date from IMGT website
- OGRDB: Track Zenodo archive version
- VDJbase: Track API response version/date

## Implementation Notes

### IMGT Provider Changes

The IMGT provider currently has `download()` raising `NotImplementedError`. Must implement:
1. Move logic from `scripts/download_imgt.py` into `IMGTProvider.download()`
2. Support version checking against IMGT release
3. Integrate with CLI command

### CLI Integration

Follow existing `sadie make` pattern:
- Located in `src/sadie/app.py` or similar
- Uses Click/Typer for CLI framework (check existing)
- Progress output during download
- Summary at completion

### Checkpoint/Resume

For fail-fast with resume capability:
- Track completed species in checkpoint file
- On re-run, skip already-completed species
- Clear checkpoint on successful full completion

### 5. Progress Reporting

**Decision**: Progress bar using rich package

- Use `rich.progress` for download progress bars
- Show species-level progress (e.g., "Downloading mouse [3/33]")
- Show file-level progress within each species
- Display ETA and download speed
- Clean output suitable for CI/CD logs

### 6. Parallel vs Sequential

**Decision**: Sequential downloads (one species at a time)

- Avoid rate limiting / IP bans from providers
- Simpler error handling and checkpointing
- More predictable progress reporting
- No `--parallel` flag needed

### 7. Post-Download Actions

**Decision**: Complete build pipeline after download

The `sadie germlines populate` command performs full pipeline:
1. Download source FASTA files
2. Build BLAST databases (`makeblastdb`)
3. Generate auxiliary files (*.aux)
4. Build internal_data directories (for IgBLAST)
5. Update organism.yaml

This reuses existing scripts:
- `build_aux_files.py`
- `build_internal_data.py`
- `update_organism_yaml.py`

### 8. Output Location

**Decision**: Always use default location

- Sources: `src/sadie/germlines/sources/{provider}/`
- Databases: `src/sadie/germlines/igblast/database/`
- No `--output-dir` flag
- No environment variable override
- Keeps package self-contained

## Out of Scope

- Automatic download on import (rejected)
- Post-install hooks (rejected)
- Partial failure tolerance (rejected - fail fast chosen)
- Subset filtering by default (all species downloaded)
- Parallel downloads (rejected - avoid rate limiting)
- Custom output directory (rejected - use defaults)

---
*Updated: 2026-01-22*
