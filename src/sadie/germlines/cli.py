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

from rich.console import Console
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
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


def is_up_to_date(provider_name: str) -> bool:
    """Check if provider data is up-to-date."""
    local = get_local_version(provider_name)
    if not local:
        return False

    current = get_current_version_string()
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
    if not force and is_up_to_date(provider_name):
        console.print(f"[yellow]{provider_name}[/yellow]: Already up-to-date, skipping")
        return {"status": "skipped", "reason": "up-to-date"}

    all_species = get_all_provider_species(provider_name)
    target_species = species if species else all_species

    if dry_run:
        console.print(
            f"[cyan]{provider_name}[/cyan]: Would download {len(target_species)} species"
        )
        for sp in target_species:
            console.print(f"  - {sp}")
        return {"status": "dry_run", "species_count": len(target_species)}

    completed = load_checkpoint(provider_name) if not force else set()
    remaining = [sp for sp in target_species if sp not in completed]

    if not remaining:
        console.print(f"[green]{provider_name}[/green]: All species already downloaded")
        clear_checkpoint(provider_name)
        return {"status": "complete", "species_count": len(target_species)}

    task = progress.add_task(
        f"[cyan]{provider_name}[/cyan]", total=len(target_species)
    )

    if completed:
        progress.update(task, completed=len(completed))
        console.print(
            f"[yellow]Resuming from checkpoint: {len(completed)}/{len(target_species)} complete[/yellow]"
        )

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

            completed.add(sp)
            results["downloaded"].append(sp)
            save_checkpoint(provider_name, completed, len(target_species))

        except Exception as e:
            logger.error(f"Failed to download {sp} from {provider_name}: {e}")
            results["failed"].append({"species": sp, "error": str(e)})

        progress.advance(task)

    if not results["failed"]:
        clear_checkpoint(provider_name)
        version = get_current_version_string()
        save_version(provider_name, version, len(target_species))

    results["species_count"] = len(results["downloaded"]) + len(results["skipped"])
    return results


def run_post_download_build(
    provider_name: str, species_list: List[str], progress: Progress
):
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

    task = progress.add_task("[cyan]Building databases...", total=3)

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
                        console.print(
                            f"[red]Invalid gene in {species}/{segment}/{chain}[/red]"
                        )
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
        table.add_row(
            prov, status, str(species_count) if species_count else "Not downloaded"
        )

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
                    console.print(
                        f"[yellow]Warning: {prov_name} validation had issues[/yellow]"
                    )

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

        status_style = (
            "green"
            if status == "success"
            else "yellow"
            if status in ["skipped", "dry_run"]
            else "red"
        )
        summary_table.add_row(
            prov_name,
            f"[{status_style}]{status}[/{status_style}]",
            str(downloaded),
            str(failed),
        )

    console.print(summary_table)

    total_failed = sum(len(r.get("failed", [])) for r in all_results.values())
    if total_failed > 0:
        console.print(
            f"\n[red]Warning: {total_failed} species failed to download[/red]"
        )
        console.print("Re-run the command to resume from checkpoint")
