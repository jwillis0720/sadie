# Phase 12: Provider Auto-Population — Execution Plan

## Overview

**Goal**: Make all germline providers (IMGT, OGRDB, VDJbase) auto-populate via CLI command

**Decisions** (from 12-CONTEXT.md):
- Trigger: `sadie germlines populate` CLI command
- Scope: All species available from each provider's API
- Errors: Fail fast with checkpoint for resume
- Updates: Version-check and update if newer
- Progress: Rich progress bars with species/file-level tracking
- Sequential: One species at a time (avoid rate limiting)
- Post-download: Complete build (BLAST DBs, aux files, internal_data)
- Location: Always use default paths (no custom output)

---

## Task Breakdown

### T080: Implement `sadie germlines populate` CLI command

**File**: `src/sadie/app.py`

**Changes**:
1. Add new Click group `germlines` under `sadie`
2. Add `populate` subcommand with options

```python
@sadie.group()
def germlines():
    """Germline database management commands."""
    pass

@germlines.command("populate")
@click.option("--provider", "-p", type=click.Choice(["imgt", "ogrdb", "vdjbase", "all"]), default="all")
@click.option("--species", "-s", multiple=True, help="Specific species to download (optional)")
@click.option("--force", "-f", is_flag=True, help="Force re-download even if up-to-date")
@click.option("--dry-run", is_flag=True, help="Show what would be downloaded without downloading")
def populate(provider: str, species: tuple, force: bool, dry_run: bool):
    """Download germline data from providers."""
    from sadie.germlines.cli import populate_germlines
    populate_germlines(provider, list(species) if species else None, force, dry_run)
```

**Acceptance Criteria**:
- `sadie germlines populate --help` shows usage
- `--provider` accepts imgt, ogrdb, vdjbase, all
- `--species` accepts multiple species names
- `--force` flag available
- `--dry-run` flag available

**Dependencies**: Add `rich` to project dependencies (likely already present)

---

### T081: Implement `IMGTProvider.download()` from existing script logic

**File**: `src/sadie/germlines/providers/imgt.py`

**Changes**:
1. Move download logic from `scripts/download_imgt.py` into provider class
2. Replace `NotImplementedError` with working implementation

```python
def download(self, species: List[str]) -> None:
    """
    Download IMGT data for species.

    Parameters
    ----------
    species : List[str]
        Species to download. If empty, downloads all SPECIES_MAP species.
    """
    from sadie.germlines.scripts.download_imgt import IMGTDownloader

    downloader = IMGTDownloader(output_dir=self.data_dir)

    if not species:
        species = list(SPECIES_MAP.keys())

    for sp in species:
        logger.info(f"Downloading IMGT data for {sp}...")
        downloader.download_species(sp)
```

**Refactoring needed in `scripts/download_imgt.py`**:
- Extract `IMGTDownloader` class with `download_species()` method
- Keep CLI as thin wrapper around class

**Acceptance Criteria**:
- `IMGTProvider().download(["human"])` works
- `IMGTProvider().download([])` downloads all species
- Progress logged during download

---

### T082: Add version tracking for IMGT releases

**File**: `src/sadie/germlines/providers/imgt.py` and `sources/imgt/VERSION.json`

**Changes**:
1. Create `VERSION.json` file in each provider's source directory
2. Check version before downloading
3. Skip if local version matches remote

```python
def _get_local_version(self) -> Optional[str]:
    version_file = self.data_dir / "VERSION.json"
    if version_file.exists():
        return json.loads(version_file.read_text()).get("version")
    return None

def _save_version(self, version: str):
    version_file = self.data_dir / "VERSION.json"
    version_file.write_text(json.dumps({
        "version": version,
        "downloaded_at": datetime.now().isoformat(),
        "provider": "imgt"
    }, indent=2))

def _get_remote_version(self) -> str:
    # Parse IMGT website for release version
    # or use release date as version
    return f"release-{datetime.now().strftime('%Y%m')}"
```

**Acceptance Criteria**:
- `VERSION.json` created after download
- Skip download if version matches (unless --force)
- Version includes download timestamp

---

### T083: Audit and download all OGRDB available species

**File**: `src/sadie/germlines/providers/ogrdb.py`

**Changes**:
1. Add method to query OGRDB API for available species
2. Ensure `download()` method works for all species

```python
def get_available_species_from_api(self) -> List[str]:
    """Query OGRDB API for all available species."""
    # OGRDB uses Zenodo archive - check what species are available
    # Return list of species names
    pass
```

**Acceptance Criteria**:
- `OGRDBProvider().download([])` downloads all available species
- Downloaded species count logged
- Data integrity verified

---

### T084: Audit and download all VDJbase available species

**File**: `src/sadie/germlines/providers/vdjbase.py`

**Changes**:
1. Add method to query VDJbase API for available species
2. Ensure `download()` method works for all species

```python
def get_available_species_from_api(self) -> List[str]:
    """Query VDJbase API for all available species."""
    # VDJbase API endpoint for species list
    pass
```

**Acceptance Criteria**:
- `VDJbaseProvider().download([])` downloads all available species
- Downloaded species count logged
- Data integrity verified

---

### T085: Add `--force` flag for re-download

**File**: `src/sadie/germlines/cli.py` (new file)

**Changes**:
1. Create CLI module with `populate_germlines()` function
2. Pass force flag to provider download methods

```python
def populate_germlines(
    provider: str,
    species: Optional[List[str]],
    force: bool,
    dry_run: bool
):
    """Main entry point for germlines populate command."""
    providers_to_run = []

    if provider == "all":
        providers_to_run = ["imgt", "ogrdb", "vdjbase"]
    else:
        providers_to_run = [provider]

    for prov_name in providers_to_run:
        prov = get_provider(prov_name)

        if not force and prov.is_up_to_date():
            click.echo(f"{prov_name}: Already up-to-date, skipping")
            continue

        if dry_run:
            click.echo(f"{prov_name}: Would download {len(species or prov.get_all_species())} species")
            continue

        prov.download(species or [], force=force)
```

**Acceptance Criteria**:
- `--force` bypasses version check
- Without `--force`, up-to-date providers skipped
- Dry run shows what would happen

---

### T086: Add checkpoint/resume for fail-fast recovery

**File**: `src/sadie/germlines/providers/base.py` and provider implementations

**Changes**:
1. Add checkpoint methods to base provider class
2. Track completed species in `.checkpoint.json`
3. Resume from checkpoint on re-run

```python
class GermlineProvider(ABC):
    def _load_checkpoint(self) -> Set[str]:
        checkpoint_file = self.data_dir / ".checkpoint.json"
        if checkpoint_file.exists():
            return set(json.loads(checkpoint_file.read_text()).get("completed", []))
        return set()

    def _save_checkpoint(self, completed: Set[str]):
        checkpoint_file = self.data_dir / ".checkpoint.json"
        checkpoint_file.write_text(json.dumps({
            "completed": list(completed),
            "updated_at": datetime.now().isoformat()
        }, indent=2))

    def _clear_checkpoint(self):
        checkpoint_file = self.data_dir / ".checkpoint.json"
        if checkpoint_file.exists():
            checkpoint_file.unlink()
```

**Acceptance Criteria**:
- Checkpoint created during download
- Re-run skips completed species
- Checkpoint cleared on full success

---

### T087: Add rich progress bars for download tracking

**File**: `src/sadie/germlines/cli.py`

**Changes**:
1. Add rich progress bars for species-level and file-level tracking
2. Show ETA and download speed

```python
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

def populate_germlines(...):
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
    ) as progress:
        species_task = progress.add_task(f"[cyan]{provider}...", total=len(species_list))

        for sp in species_list:
            progress.update(species_task, description=f"[cyan]Downloading {sp}...")
            prov.download_species(sp)
            progress.advance(species_task)
```

**Acceptance Criteria**:
- Progress bar shows current species
- Progress bar shows X/Y species completed
- ETA displayed during long downloads

---

### T088: Integrate post-download build pipeline

**File**: `src/sadie/germlines/cli.py`

**Changes**:
1. After download, run build_aux_files for downloaded species
2. Run build_internal_data for downloaded species
3. Update organism.yaml

```python
def _run_post_download_build(species_list: List[str]):
    """Run complete build pipeline after download."""
    from sadie.germlines.scripts.build_aux_files import build_aux_file
    from sadie.germlines.scripts.build_internal_data import build_internal_data
    from sadie.germlines.scripts.update_organism_yaml import update_organism_yaml

    console.print("[bold]Building BLAST databases...[/bold]")
    # makeblastdb is called during download

    console.print("[bold]Generating auxiliary files...[/bold]")
    for species in species_list:
        build_aux_file(species)

    console.print("[bold]Building internal_data...[/bold]")
    for species in species_list:
        build_internal_data(species)

    console.print("[bold]Updating organism.yaml...[/bold]")
    update_organism_yaml()
```

**Acceptance Criteria**:
- BLAST databases built after download
- Aux files generated for all species
- internal_data directories created
- organism.yaml updated with new species

---

### T089: Test CLI command with all providers

**File**: `tests/unit/germlines/test_cli.py` (new file)

**Tests**:
1. `test_populate_help()` - CLI help works
2. `test_populate_imgt_single_species()` - IMGT single species download
3. `test_populate_dry_run()` - Dry run doesn't download
4. `test_populate_force_redownload()` - Force bypasses version check
5. `test_populate_checkpoint_resume()` - Checkpoint resumes correctly

```python
def test_populate_help(cli_runner):
    result = cli_runner.invoke(sadie, ["germlines", "populate", "--help"])
    assert result.exit_code == 0
    assert "--provider" in result.output

def test_populate_dry_run(cli_runner, monkeypatch):
    result = cli_runner.invoke(sadie, ["germlines", "populate", "--dry-run"])
    assert result.exit_code == 0
    assert "Would download" in result.output
```

**Acceptance Criteria**:
- All CLI options tested
- Checkpoint behavior tested
- Force flag tested

---

### T090: Verify downloaded data integrity

**File**: `src/sadie/germlines/scripts/validate.py` (update)

**Changes**:
1. Add validation after each provider download
2. Verify FASTA files are valid
3. Verify expected gene counts

```python
def validate_provider_data(provider_name: str) -> bool:
    """Validate downloaded provider data."""
    provider = get_provider(provider_name)

    for species in provider.get_available_species():
        for segment in ["V", "D", "J"]:
            for chain in ["H", "K", "L"]:
                genes = provider.fetch_genes(species, segment, chain)
                # Verify genes have required fields
                for gene in genes:
                    assert gene.name
                    assert gene.sequence

    return True
```

**Acceptance Criteria**:
- Validation runs after download
- Invalid data causes fail-fast error
- Validation errors are clear

---

## Execution Order

```
T080 (CLI command)
    ↓
T081 (IMGT download) ←──┐
    ↓                   │
T082 (IMGT versioning)  │
    ↓                   │
T083 (OGRDB download)   │ Can run in parallel
    ↓                   │
T084 (VDJbase download) │
    ↓                   │
T085 (force flag) ──────┘
    ↓
T086 (checkpoint/resume)
    ↓
T087 (rich progress)
    ↓
T088 (post-download build)
    ↓
T089 (tests)
    ↓
T090 (validation)
```

## Estimated Effort

| Task | Effort | Risk |
|------|--------|------|
| T080 | 15 min | Low |
| T081 | 45 min | Medium (refactoring) |
| T082 | 20 min | Low |
| T083 | 30 min | Medium (API discovery) |
| T084 | 30 min | Medium (API discovery) |
| T085 | 15 min | Low |
| T086 | 30 min | Medium |
| T087 | 20 min | Low |
| T088 | 30 min | Low |
| T089 | 30 min | Low |
| T090 | 20 min | Low |

**Total**: ~5 hours

---

## Files to Create/Modify

### New Files
- `src/sadie/germlines/cli.py` — CLI logic for populate command

### Modified Files
- `src/sadie/app.py` — Add germlines group and populate command
- `src/sadie/germlines/providers/imgt.py` — Implement download()
- `src/sadie/germlines/providers/ogrdb.py` — Add version tracking
- `src/sadie/germlines/providers/vdjbase.py` — Add version tracking
- `src/sadie/germlines/providers/base.py` — Add checkpoint methods
- `tests/unit/germlines/test_cli.py` — New test file

---
*Created: 2026-01-22*
