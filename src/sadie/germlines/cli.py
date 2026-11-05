"""
Germline CLI Module
===================

Command-line interface logic for germline database management.

This module provides the core logic for the `sadie germlines populate` command,
including provider downloads, version checking, progress tracking, and
post-download build pipeline.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Set

import click

from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

logger = logging.getLogger(__name__)
console = Console()


CHECKPOINT_FILE = ".populate_checkpoint.json"
VERSION_FILE = "VERSION.json"


def get_provider(provider_name: str):
    """Get provider instance by name."""
    from .providers.imgt import IMGTProvider
    from .providers.ogrdb import OGRDBProvider
    from .providers.vdjbase import VDJbaseProvider

    providers = {
        "imgt": IMGTProvider,
        "ogrdb": OGRDBProvider,
        "vdjbase": VDJbaseProvider,
    }

    if provider_name not in providers:
        raise ValueError(f"Unknown provider: {provider_name}")

    return providers[provider_name]()


def get_provider_data_dir(provider_name: str) -> Path:
    """Get the data directory for a provider."""
    base_dir = Path(__file__).parent / "sources"
    return base_dir / provider_name


def load_checkpoint(provider_name: str) -> Set[str]:
    """Load download checkpoint for a provider."""
    checkpoint_path = get_provider_data_dir(provider_name) / CHECKPOINT_FILE
    if checkpoint_path.exists():
        try:
            data = json.loads(checkpoint_path.read_text())
            return set(data.get("completed_species", []))
        except Exception:
            pass
    return set()


def save_checkpoint(provider_name: str, completed: Set[str], total: int):
    """Save download checkpoint."""
    checkpoint_path = get_provider_data_dir(provider_name) / CHECKPOINT_FILE
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_path.write_text(
        json.dumps(
            {
                "completed_species": list(completed),
                "total_species": total,
                "updated_at": datetime.now().isoformat(),
            },
            indent=2,
        )
    )


def clear_checkpoint(provider_name: str):
    """Clear download checkpoint after successful completion."""
    checkpoint_path = get_provider_data_dir(provider_name) / CHECKPOINT_FILE
    if checkpoint_path.exists():
        checkpoint_path.unlink()


def get_local_version(provider_name: str) -> Optional[dict]:
    """Get local version info for a provider."""
    version_path = get_provider_data_dir(provider_name) / VERSION_FILE
    if version_path.exists():
        try:
            return json.loads(version_path.read_text())
        except Exception:
            pass
    return None


def save_version(provider_name: str, version: str, species_count: int):
    """Save version info after download."""
    version_path = get_provider_data_dir(provider_name) / VERSION_FILE
    version_path.parent.mkdir(parents=True, exist_ok=True)
    version_path.write_text(
        json.dumps(
            {
                "version": version,
                "downloaded_at": datetime.now().isoformat(),
                "provider": provider_name,
                "species_count": species_count,
            },
            indent=2,
        )
    )


def get_current_version_string() -> str:
    """Get current version string based on year-month."""
    return f"release-{datetime.now().strftime('%Y%m')}"


def get_upstream_version_string(provider_name: str) -> str:
    """Return the current upstream version string for a provider.

    For OGRDB this is derived from the Zenodo publication date so the check
    reflects actual data releases rather than the calendar month.
    """
    if provider_name == "ogrdb":
        try:
            from .scripts.download_ogrdb import OGRDBDownloader
            return OGRDBDownloader().get_latest_version_string()
        except Exception:
            pass
    return get_current_version_string()


def is_up_to_date(provider_name: str) -> bool:
    """Check if provider data is up-to-date."""
    local = get_local_version(provider_name)
    if not local:
        return False

    current = get_upstream_version_string(provider_name)
    return local.get("version") == current


def get_all_provider_species(provider_name: str) -> List[str]:
    """Get all available species for a provider."""
    if provider_name == "imgt":
        from .scripts.download_imgt import SPECIES_MAP

        return list(SPECIES_MAP.keys())
    elif provider_name == "ogrdb":
        return ["human", "mouse"]
    elif provider_name == "vdjbase":
        provider = get_provider(provider_name)
        return provider.get_available_species() or ["human", "rhesus_macaque"]
    return []


def populate_provider(
    provider_name: str,
    species: Optional[List[str]],
    force: bool,
    dry_run: bool,
    progress: Progress,
) -> dict:
    """
    Populate data for a single provider.

    Parameters
    ----------
    provider_name : str
        Provider name (imgt, ogrdb, vdjbase)
    species : List[str], optional
        Specific species to download, or None for all
    force : bool
        Force re-download even if up-to-date
    dry_run : bool
        Show what would be downloaded without downloading
    progress : Progress
        Rich progress instance

    Returns
    -------
    dict
        Results with species counts and status
    """
    all_species = get_all_provider_species(provider_name)
    target_species = species if species else all_species

    if dry_run:
        console.print(f"[cyan]{provider_name}[/cyan]: Would download {len(target_species)} species")
        for sp in target_species:
            console.print(f"  - {sp}")
        return {"status": "dry_run", "species_count": len(target_species)}

    if not force and is_up_to_date(provider_name):
        console.print(f"[yellow]{provider_name}[/yellow]: Already up-to-date, skipping")
        return {"status": "skipped", "reason": "up-to-date"}

    completed = load_checkpoint(provider_name) if not force else set()
    remaining = [sp for sp in target_species if sp not in completed]

    if not remaining:
        console.print(f"[green]{provider_name}[/green]: All species already downloaded")
        clear_checkpoint(provider_name)
        return {"status": "complete", "species_count": len(target_species)}

    task = progress.add_task(f"[cyan]{provider_name}[/cyan]", total=len(target_species))

    if completed:
        progress.update(task, completed=len(completed))
        console.print(f"[yellow]Resuming from checkpoint: {len(completed)}/{len(target_species)} complete[/yellow]")

    results = {"status": "success", "downloaded": [], "failed": [], "skipped": list(completed)}

    provider = get_provider(provider_name)

    for sp in remaining:
        progress.update(task, description=f"[cyan]{provider_name}[/cyan]: {sp}")

        try:
            if provider_name == "imgt":
                from .scripts.download_imgt import IMGTDownloader

                downloader = IMGTDownloader(output_dir=provider.data_dir)
                downloader.download([sp], force=force)
            else:
                provider.download([sp])

            # Verify files actually landed — 404-only runs produce empty dirs
            species_dir = provider.data_dir / sp
            fasta_files = list(species_dir.glob("*.fasta")) if species_dir.exists() else []
            if not fasta_files:
                if provider_name == "imgt":
                    from .scripts.download_imgt import SPECIES_MAP, IMGTDownloader as _D
                    imgt_name = SPECIES_MAP.get(sp, sp)
                    reason = _D(output_dir=provider.data_dir).check_ig_unavailable_reason(imgt_name, sp)
                else:
                    reason = f"{sp}: no sequences downloaded"
                logger.warning(f"{provider_name}: {reason}")
                results["failed"].append({"species": sp, "error": reason})
            else:
                completed.add(sp)
                results["downloaded"].append(sp)
                save_checkpoint(provider_name, completed, len(target_species))

        except Exception as e:
            logger.error(f"Failed to download {sp} from {provider_name}: {e}")
            results["failed"].append({"species": sp, "error": str(e)})

        progress.advance(task)

    # Only stamp the version when at least one species produced real data
    if results["downloaded"] and not results["failed"]:
        clear_checkpoint(provider_name)
        version = get_upstream_version_string(provider_name)
        save_version(provider_name, version, len(results["downloaded"]))

    results["species_count"] = len(results["downloaded"]) + len(results["skipped"])
    return results


def run_post_download_build(provider_name: str, species_list: List[str], progress: Progress):
    """
    Run post-download build pipeline.

    Builds BLAST databases, auxiliary files, and internal_data for downloaded species.

    Parameters
    ----------
    provider_name : str
        Provider name
    species_list : List[str]
        Species to build for
    progress : Progress
        Rich progress instance
    """
    if provider_name != "imgt":
        console.print(f"[dim]Post-download build only needed for IMGT[/dim]")
        return

    console.print("\n[bold]Running post-download build pipeline...[/bold]")

    task = progress.add_task("[cyan]Building databases...", total=4)

    try:
        from .scripts.build_aux_files import build_aux_file_for_species
        from .scripts.build_internal_data import build_internal_data_for_species
    except ImportError:
        console.print("[yellow]Build scripts not found, skipping post-download build[/yellow]")
        return

    progress.update(task, description="[cyan]Building auxiliary files...")
    for sp in species_list:
        try:
            build_aux_file_for_species(sp)
        except Exception as e:
            logger.warning(f"Failed to build aux file for {sp}: {e}")
    progress.advance(task)

    progress.update(task, description="[cyan]Building internal_data...")
    for sp in species_list:
        try:
            build_internal_data_for_species(sp)
        except Exception as e:
            logger.warning(f"Failed to build internal_data for {sp}: {e}")
    progress.advance(task)

    progress.update(task, description="[cyan]Building HMMs for renumbering...")
    try:
        from .renumbering_integration import LocalHMMBuilder

        hmm_builder = LocalHMMBuilder()
        hmm_count = 0
        for sp in species_list:
            for chain in ["H", "K", "L"]:
                try:
                    hmm_builder.get_hmm(sp, chain)
                    hmm_count += 1
                except Exception as e:
                    logger.debug(f"Could not build HMM for {sp} {chain}: {e}")
        if hmm_count > 0:
            console.print(f"[dim]Built {hmm_count} HMM models[/dim]")
    except ImportError:
        console.print("[yellow]HMM builder not available, skipping[/yellow]")
    progress.advance(task)

    progress.update(task, description="[cyan]Build complete")
    progress.advance(task)

    console.print("[green]Post-download build complete[/green]")


def validate_provider_data(provider_name: str) -> bool:
    """
    Validate downloaded provider data.

    Parameters
    ----------
    provider_name : str
        Provider name

    Returns
    -------
    bool
        True if validation passes
    """
    console.print(f"\n[bold]Validating {provider_name} data...[/bold]")

    provider = get_provider(provider_name)
    metadata = provider.get_metadata()

    if not metadata.species_available:
        console.print(f"[red]No species data found for {provider_name}[/red]")
        return False

    valid = True
    for species in metadata.species_available:
        if not provider.is_available(species):
            console.print(f"[red]Missing data for {species}[/red]")
            valid = False
            continue

        for segment in ["V", "D", "J"]:
            for chain in ["H", "K", "L"]:
                genes = provider.fetch_genes(species, segment, chain)
                if segment == "D" and chain != "H":
                    continue
                for gene in genes:
                    if not gene.name or not gene.sequence:
                        console.print(f"[red]Invalid gene in {species}/{segment}/{chain}[/red]")
                        valid = False
                        break

    if valid:
        console.print(f"[green]{provider_name} data validation passed[/green]")
    else:
        console.print(f"[red]{provider_name} data validation failed[/red]")

    return valid


def populate_germlines(
    provider: str,
    species: Optional[List[str]],
    force: bool,
    dry_run: bool,
):
    """
    Main entry point for germlines populate command.

    Parameters
    ----------
    provider : str
        Provider name or "all"
    species : List[str], optional
        Specific species to download
    force : bool
        Force re-download
    dry_run : bool
        Show what would happen without doing it
    """
    providers_to_run = []

    if provider == "all":
        providers_to_run = ["imgt", "ogrdb", "vdjbase"]
    else:
        providers_to_run = [provider]

    console.print("\n[bold]SADIE Germline Database Population[/bold]")
    console.print("=" * 50)

    if dry_run:
        console.print("[yellow]DRY RUN - no changes will be made[/yellow]\n")

    table = Table(title="Provider Status")
    table.add_column("Provider", style="cyan")
    table.add_column("Status")
    table.add_column("Species")

    for prov in providers_to_run:
        local_version = get_local_version(prov)
        status = "Up-to-date" if is_up_to_date(prov) else "Needs update"
        species_count = local_version.get("species_count", 0) if local_version else 0
        table.add_row(prov, status, str(species_count) if species_count else "Not downloaded")

    console.print(table)
    console.print()

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        all_results = {}

        for prov_name in providers_to_run:
            console.print(f"\n[bold]Processing {prov_name}...[/bold]")

            results = populate_provider(prov_name, species, force, dry_run, progress)
            all_results[prov_name] = results

            if not dry_run and results.get("status") == "success":
                downloaded = results.get("downloaded", [])
                if downloaded:
                    run_post_download_build(prov_name, downloaded, progress)

                if not validate_provider_data(prov_name):
                    console.print(f"[yellow]Warning: {prov_name} validation had issues[/yellow]")

    console.print("\n[bold]Summary[/bold]")
    console.print("=" * 50)

    summary_table = Table()
    summary_table.add_column("Provider", style="cyan")
    summary_table.add_column("Status")
    summary_table.add_column("Downloaded")
    summary_table.add_column("Failed")

    for prov_name, results in all_results.items():
        status = results.get("status", "unknown")
        downloaded = len(results.get("downloaded", []))
        failed = len(results.get("failed", []))

        status_style = "green" if status == "success" else "yellow" if status in ["skipped", "dry_run"] else "red"
        summary_table.add_row(
            prov_name,
            f"[{status_style}]{status}[/{status_style}]",
            str(downloaded),
            str(failed),
        )

    console.print(summary_table)

    total_failed = sum(len(r.get("failed", [])) for r in all_results.values())
    if total_failed > 0:
        console.print(f"\n[red]Warning: {total_failed} species failed to download[/red]")
        console.print("Re-run the command to resume from checkpoint")


def get_all_downloaded_species() -> List[str]:
    """
    Discover all species with downloaded source data across all providers.

    Returns
    -------
    List[str]
        Sorted unique species names that have data in any provider
    """
    species = set()
    for provider_name in ["imgt", "ogrdb", "vdjbase", "custom"]:
        data_dir = get_provider_data_dir(provider_name)
        if not data_dir.exists():
            continue
        for d in data_dir.iterdir():
            if d.is_dir() and not d.name.startswith("."):
                species.add(d.name)
    return sorted(species)


def rebuild_germlines(
    species: Optional[List[str]],
):
    """
    Rebuild germline databases from already-downloaded data (no downloading).

    Uses the Reference module to build IgBLAST databases:
    1. Auto-generates reference YAML config if missing
    2. Loads reference config via References.from_yaml()
    3. Builds IgBLAST databases via References.make_airr_database()
    4. Places output in the germlines/igblast/ directory

    Parameters
    ----------
    species : List[str], optional
        Specific species to rebuild, or None for all downloaded
    """
    import tempfile

    import yaml

    from sadie.reference.generate import generate_reference_yaml
    from sadie.reference.reference import References

    germlines_root = Path(__file__).parent
    igblast_dir = germlines_root / "igblast"

    all_downloaded = get_all_downloaded_species()

    if species:
        target = [sp for sp in species if sp in all_downloaded]
        missing = [sp for sp in species if sp not in all_downloaded]
        if missing:
            console.print(f"[yellow]Species not found locally: {', '.join(missing)}[/yellow]")
    else:
        target = all_downloaded

    if not target:
        console.print("[yellow]No downloaded germline data found. Run 'sadie germlines populate' first.[/yellow]")
        return

    console.print("\n[bold]SADIE Germline Database Rebuild[/bold]")
    console.print("=" * 50)
    console.print(f"[dim]Rebuilding {len(target)} species via Reference module (no downloads)[/dim]\n")

    # Step 1: Get or generate reference YAML config
    yaml_path = Path(__file__).parent.parent / "reference" / "data" / "reference.yml"

    if not yaml_path.exists():
        console.print("[cyan]Reference config not found, auto-generating...[/cyan]")
        generate_reference_yaml(output_path=yaml_path)
        console.print(f"[green]Generated reference config at {yaml_path}[/green]")

    # Step 2: Read the full config to filter by target species
    with open(yaml_path) as f:
        full_config = yaml.safe_load(f)

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("[cyan]Rebuilding germlines...", total=len(target))

        for sp in target:
            progress.update(task, description=f"[cyan]Rebuilding {sp} via Reference module...")
            try:
                # Find the reference name for this species in the config
                ref_name = _find_reference_name_for_species(full_config, sp)
                if ref_name is None:
                    console.print(f"  [yellow]{sp}[/yellow]: not found in reference config, skipping")
                    progress.advance(task)
                    continue

                # Build a filtered config for just this species/reference
                filtered_config = {ref_name: full_config[ref_name]}
                filtered_yaml_content = yaml.dump(filtered_config, default_flow_style=False, sort_keys=False)

                # Write to a temporary YAML file and build
                with tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False) as tmp_yaml:
                    tmp_yaml.write(filtered_yaml_content)
                    tmp_yaml_path = Path(tmp_yaml.name)

                try:
                    refs = References.from_yaml(tmp_yaml_path, use_germlines=True)
                    refs.make_airr_database(igblast_dir)
                    console.print(f"  [green]{sp}[/green]: rebuilt successfully via Reference module")
                finally:
                    tmp_yaml_path.unlink(missing_ok=True)

            except Exception as e:
                console.print(f"  [red]{sp}[/red]: failed - {e}")
                logger.error(f"Failed to rebuild {sp}: {e}")
            progress.advance(task)

    console.print("\n[green]Rebuild complete[/green]")


def _find_reference_name_for_species(config: dict, species: str) -> Optional[str]:
    """Find the reference name that contains data for a given species.

    Searches through the reference config to find a reference name that
    has the species as a sub-key under any provider.

    Parameters
    ----------
    config : dict
        Full reference YAML config.
    species : str
        Species to find.

    Returns
    -------
    str or None
        Reference name if found, None otherwise.
    """
    # First check if species name is directly a reference name
    if species in config:
        return species

    # Search through all reference names for the species
    for ref_name in config:
        for source in config[ref_name]:
            if species in config[ref_name][source]:
                return ref_name

    return None


# ---------------------------------------------------------------------------
# Click command definitions
# ---------------------------------------------------------------------------


@click.group(name="germlines", help="Germline database management commands.")
def germlines() -> None:
    pass


@germlines.command("populate")
@click.option(
    "--provider",
    "-p",
    default="all",
    type=click.Choice(["all", "imgt", "ogrdb", "vdjbase"], case_sensitive=False),
    show_default=True,
    help="Provider to populate.",
)
@click.option(
    "--species",
    "-s",
    multiple=True,
    help="Species to download (may be repeated). Defaults to all species for the provider.",
)
@click.option("--force", "-f", is_flag=True, default=False, help="Force re-download even if up-to-date.")
@click.option("--dry-run", is_flag=True, default=False, help="Show what would be downloaded without doing it.")
def populate_cmd(provider: str, species: tuple, force: bool, dry_run: bool) -> None:
    """Download and update germline sequences from a provider."""
    populate_germlines(
        provider=provider,
        species=list(species) if species else None,
        force=force,
        dry_run=dry_run,
    )


@germlines.command("status")
def status_cmd() -> None:
    """Show download status for all germline providers."""
    table = Table(title="Germline Provider Status")
    table.add_column("Provider", style="cyan")
    table.add_column("Version")
    table.add_column("Species")
    table.add_column("Up-to-date")

    for prov in ["imgt", "ogrdb", "vdjbase", "custom"]:
        local = get_local_version(prov)
        if local:
            version_str = local.get("version", "unknown")
            species_count = str(local.get("species_count", "?"))
            current = is_up_to_date(prov)
            status_str = "[green]yes[/green]" if current else "[yellow]no[/yellow]"
        else:
            version_str = "not downloaded"
            species_count = "-"
            status_str = "[red]no[/red]"
        table.add_row(prov, version_str, species_count, status_str)

    console.print(table)


@germlines.command("rebuild")
@click.option(
    "--species",
    "-s",
    multiple=True,
    help="Species to rebuild (may be repeated). Defaults to all downloaded species.",
)
def rebuild_cmd(species: tuple) -> None:
    """Rebuild IgBLAST databases from already-downloaded germline data via the Reference module.

    Runs References.from_yaml() -> make_airr_database() to regenerate BLAST
    databases in germlines/igblast/ without re-downloading source FASTA files.
    """
    rebuild_germlines(species=list(species) if species else None)
