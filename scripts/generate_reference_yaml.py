#!/usr/bin/env python3
"""
Generate reference.yml file for SADIE reference databases.

This script queries the Germline Gene Gateway (G3) API to help create
or update reference.yml files used by the 'sadie make-all' command.

Usage:
    python scripts/generate_reference_yaml.py
    python scripts/generate_reference_yaml.py --species human --output my_reference.yml
    python scripts/generate_reference_yaml.py --dry-run --verbose
    python scripts/generate_reference_yaml.py --all-functional --species human
    python scripts/generate_reference_yaml.py --interactive
"""

import argparse
import json
import shutil
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import requests
import yaml

try:
    from rich import print as rprint
    from rich.console import Console
    from rich.panel import Panel
    from rich.progress import Progress, SpinnerColumn, TextColumn
    from rich.prompt import Confirm, Prompt
    from rich.table import Table
except ImportError as e:
    print(f"Error: Missing required dependencies. Please install: {e}")
    print("Run: pip install rich requests pyyaml")
    sys.exit(1)

console = Console()

# G3 API Configuration
G3_API_BASE = "https://g3.jordanrwillis.com/api/v1"
DEFAULT_OUTPUT = Path("reference.yml")


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate reference.yml file for SADIE reference databases",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --name human --species human --all-functional
  %(prog)s --interactive
  %(prog)s --update --output existing.yml --species mouse
  %(prog)s --name custom --species human --segments V J --source imgt
  %(prog)s --dry-run --verbose
        """
    )

    # Basic options
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output YAML file path (default: {DEFAULT_OUTPUT})"
    )

    parser.add_argument(
        "--name",
        help="Reference name (e.g., human, mouse, custom)"
    )

    parser.add_argument(
        "--species",
        help="Species to query (e.g., human, mouse, rat, rabbit, dog, macaque)"
    )

    parser.add_argument(
        "--source",
        choices=["imgt", "custom", "all"],
        default="imgt",
        help="Gene source database (default: imgt)"
    )

    parser.add_argument(
        "--segments",
        nargs="+",
        choices=["V", "D", "J", "C"],
        help="Gene segments to include (default: V D J)"
    )

    # Convenience options
    parser.add_argument(
        "--all-functional",
        action="store_true",
        help="Include all functional genes for the species"
    )

    parser.add_argument(
        "--common-genes",
        action="store_true",
        help="Include only commonly used genes (recommended subset)"
    )

    # Mode options
    parser.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Interactive gene selection mode"
    )

    parser.add_argument(
        "--update",
        action="store_true",
        help="Update existing reference.yml file"
    )

    # Standard options
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output"
    )

    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Skip backing up existing reference file"
    )

    parser.add_argument(
        "--limit",
        type=int,
        default=1000,
        help="Maximum number of genes to query from G3 API (default: 1000)"
    )

    return parser.parse_args()


def query_g3_api(
    species: Optional[str] = None,
    source: Optional[str] = None,
    segment: Optional[str] = None,
    functional: Optional[bool] = None,
    limit: int = 1000,
    verbose: bool = False
) -> List[Dict[str, Any]]:
    """Query the G3 API for available genes."""
    params = {"limit": limit}

    if species:
        params["common"] = species
    if source and source != "all":
        params["source"] = source
    if segment:
        params["segment"] = segment
    if functional is not None:
        params["functional"] = str(functional).lower()

    url = f"{G3_API_BASE}/genes"

    if verbose:
        console.print(f"[dim]Querying: {url}")
        console.print(f"[dim]Params: {params}")

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        console.print(f"[red]Error querying G3 API: {e}[/red]")
        return []


def get_available_species(verbose: bool = False) -> List[str]:
    """Get list of available species from G3 API."""
    # Query with limit=1 for each common species to check availability
    common_species = ["human", "mouse", "rat", "rabbit", "dog", "macaque"]
    available = []

    for species in common_species:
        genes = query_g3_api(species=species, limit=1, verbose=verbose)
        if genes:
            available.append(species)

    return available


def organize_genes_by_type(genes: List[Dict[str, Any]]) -> Dict[str, Dict[str, List[str]]]:
    """Organize genes by source and species."""
    organized = defaultdict(lambda: defaultdict(list))

    for gene in genes:
        source = gene.get("source", "unknown")
        species = gene.get("common", "unknown")
        gene_name = gene.get("gene", "")

        if gene_name:
            organized[source][species].append(gene_name)

    # Convert to regular dict and sort gene lists
    result = {}
    for source in organized:
        result[source] = {}
        for species in organized[source]:
            result[source][species] = sorted(organized[source][species])

    return result


def filter_common_genes(genes: List[str]) -> List[str]:
    """Filter to include only commonly used genes based on patterns."""
    common_patterns = {
        "V": ["IGHV3-", "IGHV4-", "IGHV1-", "IGKV1-", "IGKV3-", "IGLV1-", "IGLV2-"],
        "D": ["IGHD3-", "IGHD2-", "IGHD1-"],
        "J": ["IGHJ4", "IGHJ6", "IGKJ", "IGLJ"],
    }

    filtered = []
    for gene in genes:
        for patterns in common_patterns.values():
            if any(pattern in gene for pattern in patterns):
                filtered.append(gene)
                break

    return filtered


def create_reference_structure(
    name: str,
    genes_by_source: Dict[str, Dict[str, List[str]]]
) -> Dict[str, Any]:
    """Create the reference YAML structure."""
    return {name: genes_by_source}


def load_existing_reference(file_path: Path) -> Dict[str, Any]:
    """Load and parse existing reference.yml file."""
    if not file_path.exists():
        return {}

    try:
        with open(file_path, 'r') as f:
            return yaml.safe_load(f) or {}
    except Exception as e:
        console.print(f"[red]Error loading existing reference: {e}[/red]")
        return {}


def merge_references(
    existing: Dict[str, Any],
    new: Dict[str, Any],
    verbose: bool = False
) -> Tuple[Dict[str, Any], Dict[str, int]]:
    """Merge new reference with existing, returning merged dict and stats."""
    merged = existing.copy()
    stats = {"added": 0, "existing": 0, "total": 0}

    for name, sources in new.items():
        if name not in merged:
            merged[name] = {}

        for source, species_genes in sources.items():
            if source not in merged[name]:
                merged[name][source] = {}

            for species, genes in species_genes.items():
                if species not in merged[name][source]:
                    merged[name][source][species] = []

                existing_genes = set(merged[name][source][species])

                for gene in genes:
                    stats["total"] += 1
                    if gene not in existing_genes:
                        merged[name][source][species].append(gene)
                        stats["added"] += 1
                        if verbose:
                            console.print(f"[green]+ {gene}[/green]")
                    else:
                        stats["existing"] += 1
                        if verbose:
                            console.print(f"[dim]= {gene} (already exists)[/dim]")

                # Sort genes
                merged[name][source][species].sort()

    return merged, stats


def display_gene_summary(genes: List[Dict[str, Any]]) -> None:
    """Display summary of genes in a table."""
    # Count by segment and functional status
    segment_counts = defaultdict(lambda: {"functional": 0, "non_functional": 0})

    for gene in genes:
        segment = gene.get("segment", "Unknown")
        # Check different possible keys for functional status
        is_functional = gene.get("functional", gene.get("is_functional", True))

        if is_functional:
            segment_counts[segment]["functional"] += 1
        else:
            segment_counts[segment]["non_functional"] += 1

    table = Table(title="Gene Summary", show_header=True)
    table.add_column("Segment", style="cyan")
    table.add_column("Functional", style="green", justify="right")
    table.add_column("Non-Functional", style="yellow", justify="right")
    table.add_column("Total", style="magenta", justify="right")

    total_functional = 0
    total_non_functional = 0

    for segment in ["V", "D", "J", "C"]:
        if segment in segment_counts:
            counts = segment_counts[segment]
            functional = counts["functional"]
            non_functional = counts["non_functional"]
            total = functional + non_functional

            total_functional += functional
            total_non_functional += non_functional

            table.add_row(
                segment,
                str(functional),
                str(non_functional),
                str(total)
            )

    # Add total row
    table.add_row(
        "[bold]Total[/bold]",
        f"[bold]{total_functional}[/bold]",
        f"[bold]{total_non_functional}[/bold]",
        f"[bold]{total_functional + total_non_functional}[/bold]",
        style="bold"
    )

    console.print(table)


def interactive_mode() -> Optional[Dict[str, Any]]:
    """Interactive mode for selecting genes."""
    console.print(Panel.fit(
        "[bold blue]Interactive Reference Generation[/bold blue]\n"
        "Follow the prompts to select genes for your reference",
        border_style="blue"
    ))

    # Get reference name
    name = Prompt.ask("\n[cyan]Reference name[/cyan]", default="custom")

    # Get available species
    console.print("\n[cyan]Checking available species...[/cyan]")
    available_species = get_available_species()

    if not available_species:
        console.print("[red]No species available from G3 API[/red]")
        return None

    console.print(f"Available species: {', '.join(available_species)}")
    species = Prompt.ask(
        "[cyan]Select species[/cyan]",
        choices=available_species,
        default=available_species[0]
    )

    # Get source
    source = Prompt.ask(
        "\n[cyan]Gene source[/cyan]",
        choices=["imgt", "custom", "all"],
        default="imgt"
    )

    # Get segments
    console.print("\n[cyan]Select gene segments (space-separated)[/cyan]")
    segments_input = Prompt.ask(
        "Segments",
        default="V D J"
    )
    segments = segments_input.upper().split()

    # Functional only?
    functional_only = Confirm.ask(
        "\n[cyan]Include functional genes only?[/cyan]",
        default=True
    )

    # Query genes
    all_genes = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True
    ) as progress:

        for segment in segments:
            task = progress.add_task(f"Querying {segment} genes...", total=None)

            genes = query_g3_api(
                species=species,
                source=source if source != "all" else None,
                segment=segment,
                functional=functional_only if functional_only else None,
                limit=1000
            )

            all_genes.extend(genes)
            progress.update(task, completed=True)

    if not all_genes:
        console.print("[red]No genes found matching criteria[/red]")
        return None

    # Display summary
    console.print(f"\n[green]Found {len(all_genes)} genes[/green]")
    display_gene_summary(all_genes)

    # Common genes filter?
    if len(all_genes) > 50:
        use_common = Confirm.ask(
            "\n[cyan]Filter to commonly used genes only?[/cyan]",
            default=False
        )

        if use_common:
            gene_names = [g["gene"] for g in all_genes]
            common_genes = filter_common_genes(gene_names)
            console.print(f"[yellow]Filtered to {len(common_genes)} common genes[/yellow]")

            # Filter all_genes to match
            all_genes = [g for g in all_genes if g["gene"] in common_genes]

    # Final confirmation
    if not Confirm.ask("\n[cyan]Create reference with these genes?[/cyan]", default=True):
        return None

    # Organize and return
    organized = organize_genes_by_type(all_genes)
    return create_reference_structure(name, organized)


def backup_file(file_path: Path) -> Optional[Path]:
    """Create a timestamped backup of a file if it exists."""
    if not file_path.exists():
        return None

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = file_path.with_suffix(f".bak.{timestamp}.yml")

    shutil.copy2(file_path, backup_path)
    return backup_path


def generate_reference_yaml(
    args: argparse.Namespace
) -> None:
    """Main function to generate reference YAML."""

    # Interactive mode
    if args.interactive:
        reference_data = interactive_mode()
        if not reference_data:
            console.print("[yellow]Reference generation cancelled[/yellow]")
            return

    # Command line mode
    else:
        # Validate required arguments
        if not args.all_functional and not args.update:
            if not args.name or not args.species:
                console.print("[red]Error: --name and --species required (or use --interactive)[/red]")
                sys.exit(1)

        # Set defaults
        if args.all_functional and not args.name:
            args.name = args.species

        segments = args.segments or ["V", "D", "J"]

        # Query genes
        all_genes = []

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
            transient=True
        ) as progress:

            task = progress.add_task(
                f"Querying G3 API for {args.species} genes...",
                total=None
            )

            if args.all_functional:
                # Query all segments
                for segment in segments:
                    genes = query_g3_api(
                        species=args.species,
                        source=args.source if args.source != "all" else None,
                        segment=segment,
                        functional=True,
                        limit=args.limit,
                        verbose=args.verbose
                    )
                    all_genes.extend(genes)
            else:
                # Query with specified parameters
                for segment in segments:
                    genes = query_g3_api(
                        species=args.species,
                        source=args.source if args.source != "all" else None,
                        segment=segment,
                        functional=None,
                        limit=args.limit,
                        verbose=args.verbose
                    )
                    all_genes.extend(genes)

            progress.update(task, completed=True)

        if not all_genes:
            console.print(f"[red]No genes found for {args.species}[/red]")
            return

        console.print(f"\n[green]Found {len(all_genes)} genes[/green]")

        # Apply common genes filter if requested
        if args.common_genes:
            gene_names = [g["gene"] for g in all_genes]
            common_genes = filter_common_genes(gene_names)
            all_genes = [g for g in all_genes if g["gene"] in common_genes]
            console.print(f"[yellow]Filtered to {len(common_genes)} common genes[/yellow]")

        # Display summary
        display_gene_summary(all_genes)

        # Organize genes
        organized = organize_genes_by_type(all_genes)
        reference_data = create_reference_structure(args.name, organized)

    # Handle update mode
    if args.update and args.output.exists():
        existing = load_existing_reference(args.output)
        reference_data, stats = merge_references(
            existing,
            reference_data,
            verbose=args.verbose
        )

        console.print(f"\n[cyan]Update Summary:[/cyan]")
        console.print(f"  Added: {stats['added']}")
        console.print(f"  Existing: {stats['existing']}")
        console.print(f"  Total: {stats['total']}")

    # Dry run
    if args.dry_run:
        console.print("\n[yellow]DRY RUN - No files will be written[/yellow]")
        console.print("\n[cyan]Would write:[/cyan]")
        console.print(yaml.dump(reference_data, default_flow_style=False, sort_keys=False))
        return

    # Backup existing file
    if args.output.exists() and not args.no_backup:
        backup_path = backup_file(args.output)
        if backup_path:
            console.print(f"\n✓ Backed up existing file to {backup_path.name}")

    # Write output
    try:
        with open(args.output, 'w') as f:
            yaml.dump(reference_data, f, default_flow_style=False, sort_keys=False)

        console.print(f"\n[green]✓ Generated {args.output}[/green]")

        # Count total genes
        total_genes = sum(
            len(genes)
            for sources in reference_data.values()
            for species_dict in sources.values()
            for genes in species_dict.values()
        )

        console.print(f"[green]Total genes in reference: {total_genes}[/green]")

        # Remind about next step
        console.print("\n[dim]Next step: Run 'sadie make-all' to build the reference database[/dim]")

    except Exception as e:
        console.print(f"[red]Error writing file: {e}[/red]")
        sys.exit(1)


def main():
    """Main execution function."""
    args = parse_arguments()

    # Display header
    console.print(Panel.fit(
        "[bold blue]Reference YAML Generation Tool[/bold blue]\n"
        "Create reference.yml for SADIE databases",
        border_style="blue"
    ))

    try:
        generate_reference_yaml(args)
    except KeyboardInterrupt:
        console.print("\n[yellow]Operation cancelled by user[/yellow]")
        sys.exit(1)
    except Exception as e:
        console.print(f"\n[red]Error: {e}[/red]")
        if args.verbose:
            import traceback
            console.print(f"[dim]{traceback.format_exc()}[/dim]")
        sys.exit(1)


if __name__ == "__main__":
    main()
