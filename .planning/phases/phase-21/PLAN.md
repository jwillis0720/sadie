# Phase 21: Build CLI

**Goal:** Add CLI command to build IgBLAST database from reference.yml

**Depends on:** Phase 20 (use_germlines param in References.from_yaml)

**Requirements:**
- CLI-01: Add `sadie reference build <yaml> --output <path>` command
- CLI-02: Build generates complete IgBLAST database structure
- CLI-03: Progress output during build

## Pre-Implementation Analysis

### Existing CLI Structure
- `src/sadie/app.py` contains all CLI commands using Click
- `@sadie.group() def reference()` already exists (line 231)
- `@reference.command("make")` already exists with similar functionality
- Pattern: Commands are added directly to app.py, not separate cli.py files

### Key Differences from Existing `sadie reference make`:
| Aspect | `make` (existing) | `build` (new) |
|--------|-------------------|---------------|
| YAML path | `--reference/-d` option | Positional `<yaml>` argument |
| Output path | `--outpath/-o` option | `--output` option |
| Progress | Basic echo statements | Structured: "Loading YAML...", "Fetching genes...", "Building databases...", "Complete" |
| Use germlines | Implicit (G3 API) | Explicit `--use-germlines` flag (Phase 20 integration) |

### Output Structure (from make_airr_database):
```
<output>/
├── Ig/
│   ├── blastdb/          # BLAST databases per species
│   │   └── {name}/
│   │       ├── {name}_V.*
│   │       ├── {name}_D.*
│   │       └── {name}_J.*
│   └── internal_data/    # Annotation files
│       └── {name}/
│           └── {name}.ndm.imgt
├── aux_db/
│   └── imgt/
│       └── {name}_gl.aux
└── .references_dataframe.csv.gz
```

---

## Tasks

### Task 1: Add `build` command to reference group
**File:** `src/sadie/reference/app.py`
**Location:** After existing `make_igblast_reference` command (~line 280)

**Code to add:**

```python
@reference.command("build")
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=4,
    help="Verbosity level, ex. -vvvvv for debug level logging",
)
@click.option(
    "--output",
    "-o",
    required=True,
    help="Output path for IgBLAST database structure",
    type=click.Path(resolve_path=True, dir_okay=True, writable=True),
)
@click.option(
    "--use-germlines",
    is_flag=True,
    default=False,
    help="Use local germlines module instead of G3 API",
)
@click.argument(
    "yaml_path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
)
def build_reference(verbose: int, output: str, use_germlines: bool, yaml_path: str) -> None:
    """Build IgBLAST database from reference.yml configuration.

    Creates complete database structure including blast databases, internal
    annotation files, and auxiliary files needed for AIRR annotation.

    \b
    Examples:
        sadie reference build reference.yml --output ./db
        sadie reference build my_refs.yml -o /data/germlines --use-germlines
    """
    import sys
    from pathlib import Path

    # Set logging level
    numeric_level = getVerbosityLevel(verbose)
    logging.basicConfig(level=numeric_level)

    output_path = Path(output)

    try:
        # Progress: Loading YAML
        click.echo("Loading YAML...")
        reference_object = References.from_yaml(Path(yaml_path), use_germlines=use_germlines)

        # Progress: Fetching genes
        click.echo("Fetching genes...")
        # get_dataframe() triggers the actual gene fetching from G3/germlines
        _ = reference_object.get_dataframe()

        # Ensure output directory exists
        output_path.mkdir(parents=True, exist_ok=True)

        # Progress: Building databases
        click.echo("Building databases...")
        germline_path = reference_object.make_airr_database(output_path)

        # Progress: Complete
        click.echo("Complete")
        click.echo(f"Database written to: {germline_path}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)
```

### Task 2: Update from_yaml signature for use_germlines
**File:** `src/sadie/reference/reference.py`
**Note:** This task depends on Phase 20 completing INT-01. If Phase 20 is not complete, implement minimal version.

**Check if Phase 20 complete:** Look for `use_germlines` parameter in `from_yaml()` method.

**If not present (Phase 20 not complete), add minimal support:**
```python
@staticmethod
def from_yaml(yaml_path: Optional[Path] = None, use_germlines: bool = False) -> "References":
    """Parse a yaml file into a references file object

    Parameters
    ----------
    yaml_path : Path to yaml file
    use_germlines : bool, optional
        If True, use local germlines module instead of G3 API. Defaults to False.

    Returns
    -------
    Reference - Reference Object
    """
    yaml_ref_object = YamlRef(yaml_path)
    yaml_ref = yaml_ref_object.yaml
    references_object = References()

    for name in yaml_ref:
        reference_object = Reference(use_germlines=use_germlines)  # Pass flag
        for source in yaml_ref.get(name):
            for species in yaml_ref.get(name).get(source):
                logger.info(f"Adding {species} from {source} to {name}")
                list_of_genes: List[str] = yaml_ref[name][source][species]
                reference_object.add_genes(species, source, list_of_genes)
        references_object.add_reference(name, reference_object)
    return references_object
```

### Task 3: Verify output structure
**Manual verification steps:**
1. Run `sadie reference build reference.yml --output ./test_db`
2. Check directory structure:
   - `./test_db/Ig/blastdb/` exists with species subdirs
   - `./test_db/Ig/internal_data/` exists with `.ndm.imgt` files
   - `./test_db/aux_db/imgt/` exists with `*_gl.aux` files
   - `./test_db/.references_dataframe.csv.gz` exists
3. Verify exit code 0
4. Test error case (bad YAML path) returns non-zero exit code

### Task 4: Add integration test
**File:** `tests/integration/reference/test_cli.py` (create if needed)

```python
"""Integration tests for reference CLI build command."""

import subprocess
import tempfile
from pathlib import Path

import pytest


class TestReferenceBuildCLI:
    """Tests for `sadie reference build` command."""

    def test_build_creates_database_structure(self):
        """Verify build creates expected directory structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yaml_path = Path(__file__).parent.parent.parent.parent / "src/sadie/reference/data/reference.yml"
            
            result = subprocess.run(
                ["sadie", "reference", "build", str(yaml_path), "--output", tmpdir],
                capture_output=True,
                text=True,
            )
            
            # Check exit code
            assert result.returncode == 0, f"Build failed: {result.stderr}"
            
            # Check progress output
            assert "Loading YAML..." in result.stdout
            assert "Fetching genes..." in result.stdout
            assert "Building databases..." in result.stdout
            assert "Complete" in result.stdout
            
            # Check directory structure
            output_path = Path(tmpdir)
            assert (output_path / "Ig" / "blastdb").is_dir()
            assert (output_path / "Ig" / "internal_data").is_dir()
            assert (output_path / "aux_db").is_dir()
            assert (output_path / ".references_dataframe.csv.gz").is_file()

    def test_build_invalid_yaml_returns_error(self):
        """Verify build returns non-zero exit code for invalid YAML."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                ["sadie", "reference", "build", "/nonexistent/file.yml", "--output", tmpdir],
                capture_output=True,
                text=True,
            )
            
            # Should fail with non-zero exit code
            assert result.returncode != 0

    def test_build_help_shows_usage(self):
        """Verify --help shows command documentation."""
        result = subprocess.run(
            ["sadie", "reference", "build", "--help"],
            capture_output=True,
            text=True,
        )
        
        assert result.returncode == 0
        assert "Build IgBLAST database" in result.stdout
        assert "--output" in result.stdout
        assert "yaml_path" in result.stdout.lower()
```

---

## Success Criteria Verification

| Criterion | How to Verify |
|-----------|---------------|
| 1. `sadie reference build reference.yml --output ./db` creates database directory | Run command, check exit code 0 |
| 2. Output contains `Ig/blastdb/`, `Ig/internal_data/`, `aux_db/`, `.references_dataframe.csv.gz` | `ls -R ./db` |
| 3. Progress output shows structured messages | Check stdout for exact strings |
| 4. Exit code 0 on success, non-zero with error message on failure | Test both paths |
| 5. Database structure identical to `References.make_airr_database()` | Compare to existing `sadie reference make` output |

---

## Dependency Notes

- **Phase 20 dependency:** Task 2 requires `use_germlines` param in `from_yaml()`. If Phase 20 is not complete, implement minimal version as shown.
- **No pyproject.toml changes needed:** The `reference` group already exists in the CLI, so no new entry point is required.
- **No src/sadie/reference/cli.py needed:** Following existing pattern of adding commands directly to `app.py`.

---

## Execution Order

1. **Task 1** — Add `build` command (main implementation)
2. **Task 2** — Ensure `from_yaml()` supports `use_germlines` (check Phase 20 status first)
3. **Task 3** — Manual verification of output structure
4. **Task 4** — Add integration tests (optional, can be done in parallel with Task 3)

**Estimated complexity:** 2 tasks (Task 1 + Task 2), verification is lightweight

---
*Created: 2026-01-23*
