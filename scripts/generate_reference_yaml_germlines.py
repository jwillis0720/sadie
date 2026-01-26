#!/usr/bin/env python3
"""
Generate reference.yml file using the SADIE germlines module.

This script uses the local germlines database (IMGT, OGRDB, VDJbase) to generate
reference.yml files used by 'sadie make-all' and 'sadie reference build' commands.

Unlike the legacy G3-based script, this version:
- Works offline (no network required after initial germlines population)
- Includes alleles from VDJbase and OGRDB (not just IMGT)
- Supports the full set of species available in the germlines module

Usage:
    python scripts/generate_reference_yaml_germlines.py --generate-all
    python scripts/generate_reference_yaml_germlines.py --species human --output my_reference.yml
    python scripts/generate_reference_yaml_germlines.py --dry-run --verbose
    python scripts/generate_reference_yaml_germlines.py --interactive
"""

import argparse
import shutil
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

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
    print("Run: pip install rich pyyaml")
    sys.exit(1)

# Import SADIE germlines module
try:
    from sadie.germlines import GermlineManager
    from sadie.germlines.models import GermlineGene
except ImportError as e:
    print(f"Error: SADIE germlines module not available: {e}")
    print("Make sure SADIE is installed: pip install -e .")
    sys.exit(1)

console = Console()

DEFAULT_OUTPUT = Path("reference.yml")

# Reference configuration - defines how each reference should be generated
# This encodes the patterns found in the manually curated reference.yml
REFERENCE_CONFIG = {
    # Standard references: name=species, include ALL functional genes from all providers
    "dog": {
        "type": "standard",
        "species": "dog",
        "providers": ["vdjbase", "ogrdb", "imgt"],
    },
    "human": {
        "type": "standard",
        "species": "human",
        "providers": ["vdjbase", "ogrdb", "imgt"],
    },
    "mouse": {
        "type": "standard",
        "species": "mouse",
        "providers": ["vdjbase", "ogrdb", "imgt"],
    },
    "rabbit": {
        "type": "standard",
        "species": "rabbit",
        "providers": ["vdjbase", "ogrdb", "imgt"],
    },
    "rat": {
        "type": "standard",
        "species": "rat",
        "providers": ["vdjbase", "ogrdb", "imgt"],
    },
    "macaque": {
        "type": "standard",
        "species": "macaque",
        "providers": ["vdjbase", "ogrdb", "imgt"],
    },
    # Chimeric references: combine curated human subset with full mouse
    "clk": {
        "type": "chimeric",
        "components": {
            "human": {
                "providers": ["imgt"],
                "genes": [
                    "IGHV1-2*02",
                    "IGKV1-33*01",
                    "IGKJ3*01",
                    "IGKJ4*01",
                    "IGKJ4*02",
                    "IGHJ2*01",
                    "IGHJ3*02",
                    "IGHJ5*02",
                    "IGHD3-10*01",
                    "IGHD3-16*02",
                    "IGHD6-13*01",
                    "IGKV1-5*03",
                    "IGHJ4*02",
                    "IGHD3-9*01",
                    "IGLV2-11*01",
                    "IGLJ1*01",
                ],
            },
            "mouse": {
                "providers": ["vdjbase", "ogrdb", "imgt"],
                "genes": "ALL_FUNCTIONAL",
            },
        },
    },
    "se09": {
        "type": "chimeric",
        "components": {
            "human": {
                "providers": ["imgt"],
                "genes": [
                    "IGHV1-2*02",
                    "IGKV1-33*01",
                    "IGHJ2*01",
                ],
            },
            "mouse": {
                "providers": ["vdjbase", "ogrdb", "imgt"],
                "genes": "ALL_FUNCTIONAL",
            },
        },
    },
}

# Segments and chains to query
SEGMENTS = ["V", "D", "J"]
CHAINS = ["H", "K", "L"]

# Provider priority for deduplication (higher priority sources keep the gene)
PROVIDER_PRIORITY = ["vdjbase", "ogrdb", "imgt", "custom"]

# Maximum gene name length (IgBLAST local ID limit is 50)
MAX_GENE_NAME_LENGTH = 50

# Minimum V gene sequence length (full IGHV is ~296bp, truncated genes cause issues)
MIN_V_GENE_LENGTH = 280


def deduplicate_reference(reference: Dict[str, Any], verbose: bool = False) -> Dict[str, Any]:
    """
    Remove duplicate genes within a reference, keeping only the highest-priority source.

    Deduplication is done PER SPECIES - the same gene name in different species
    (e.g., IGHJ3*02 in human vs mouse) are NOT duplicates since they have different sequences.

    Args:
        reference: Dict of {source: {species: [genes]}}
        verbose: Print deduplication info

    Returns:
        Deduplicated reference with same structure
    """
    # First, collect all genes by species across all sources
    # Structure: {species: {gene_name: source}}
    gene_to_source: Dict[str, Dict[str, str]] = defaultdict(dict)

    # Process sources in priority order
    for source in PROVIDER_PRIORITY:
        if source not in reference:
            continue
        for species, genes in reference[source].items():
            if genes is None:
                continue
            for gene in genes:
                # Only set if not already seen (higher priority source wins)
                if gene not in gene_to_source[species]:
                    gene_to_source[species][gene] = source

    # Now rebuild the reference with only the winning genes per source
    result: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    duplicates_removed: Dict[str, int] = defaultdict(int)

    for source in PROVIDER_PRIORITY:
        if source not in reference:
            continue
        for species, genes in reference[source].items():
            if genes is None:
                continue
            for gene in genes:
                if gene_to_source[species].get(gene) == source:
                    result[source][species].append(gene)
                else:
                    duplicates_removed[f"{source}/{species}"] += 1

    # Report duplicates removed
    if verbose:
        for key, count in duplicates_removed.items():
            if count > 0:
                console.print(f"[dim]  Removed {count} duplicates from {key}[/dim]")

    # Convert to regular dicts and sort gene lists
    return {
        source: {species: sorted(genes) for species, genes in species_data.items() if genes}
        for source, species_data in result.items()
        if any(species_data.values())
    }


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate reference.yml file using SADIE germlines module",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --generate-all                    # Generate full reference.yml with all references
  %(prog)s --generate-all --output ref.yml   # Generate to custom path
  %(prog)s --species human --name human      # Generate single species
  %(prog)s --species macaque --name macaque  # Generate macaque with VDJbase alleles
  %(prog)s --dry-run --verbose               # Show what would be generated
  %(prog)s --interactive                     # Interactive mode
  %(prog)s --list-species                    # List available species
        """,
    )

    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output YAML file path (default: {DEFAULT_OUTPUT})",
    )

    parser.add_argument(
        "--generate-all",
        action="store_true",
        help="Generate complete reference.yml with all configured references",
    )

    parser.add_argument(
        "--name",
        help="Reference name (e.g., human, mouse, clk)",
    )

    parser.add_argument(
        "--species",
        help="Species to query (e.g., human, mouse, macaque)",
    )

    parser.add_argument(
        "--providers",
        nargs="+",
        choices=["vdjbase", "ogrdb", "imgt", "custom"],
        default=["vdjbase", "ogrdb", "imgt"],
        help="Germline providers to include (default: vdjbase ogrdb imgt)",
    )

    parser.add_argument(
        "--segments",
        nargs="+",
        choices=["V", "D", "J", "C"],
        default=["V", "D", "J", "C"],
        help="Gene segments to include (default: V D J C)",
    )

    parser.add_argument(
        "--functional-only",
        action="store_true",
        default=True,
        help="Include only functional genes (default: True)",
    )

    parser.add_argument(
        "--include-non-functional",
        action="store_true",
        help="Include non-functional genes as well",
    )

    parser.add_argument(
        "--interactive",
        "-i",
        action="store_true",
        help="Interactive mode for gene selection",
    )

    parser.add_argument(
        "--update",
        action="store_true",
        help="Update existing reference.yml file (merge with existing)",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Show detailed output",
    )

    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Skip backing up existing reference file",
    )

    parser.add_argument(
        "--list-species",
        action="store_true",
        help="List available species from germlines module",
    )

    parser.add_argument(
        "--compare",
        type=Path,
        help="Compare generated output with existing reference.yml",
    )

    return parser.parse_args()


def get_available_species(providers: List[str], verbose: bool = False) -> Dict[str, List[str]]:
    """Get available species from germlines module for each provider."""
    # The normalized directory is organized by species, not provider
    # We need to check each species to see which providers have data
    from sadie.germlines import get_germlines_base_dir

    normalized_dir = get_germlines_base_dir() / "normalized"

    if not normalized_dir.exists():
        if verbose:
            console.print("[yellow]Normalized directory not found[/yellow]")
        return {p: [] for p in providers}

    # Get all species
    all_species = [d.name for d in normalized_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]

    # For simplicity, report all species as available from all providers
    # The actual provider filtering happens when we query genes
    species_by_provider: Dict[str, List[str]] = {}
    for provider in providers:
        species_by_provider[provider] = sorted(all_species)
        if verbose:
            console.print(f"[dim]{provider}: {len(all_species)} species available[/dim]")

    return species_by_provider


def get_all_species() -> List[str]:
    """Get list of all available species from normalized directory."""
    from sadie.germlines import get_germlines_base_dir

    normalized_dir = get_germlines_base_dir() / "normalized"

    if not normalized_dir.exists():
        return []

    return sorted([d.name for d in normalized_dir.iterdir() if d.is_dir() and not d.name.startswith('.')])


def is_valid_gene(gene: "GermlineGene", segment: str) -> bool:
    """
    Check if a gene is valid for inclusion in reference.yml.

    V genes require gapped sequences for IgBLAST internal_data.
    D and J genes don't require gapping.
    """
    # Gene name length check
    if len(gene.name) > MAX_GENE_NAME_LENGTH:
        return False

    # V genes need gapped sequences (attribute is sequence_gapped)
    if segment == "V":
        gapped = gene.sequence_gapped
        if not gapped or len(gapped.replace(".", "")) < 50:
            return False

    return True


def get_all_genes_for_species(
    species: str,
    providers: List[str],
    segments: List[str] = SEGMENTS,
    functional_only: bool = True,
    verbose: bool = False,
) -> Dict[str, List[str]]:
    """
    Get all genes for a species from specified providers.

    Returns dict organized by source for YAML output:
    {
        "imgt": ["IGHV1-2*01", "IGHV1-2*02", ...],
        "vdjbase": ["IGHV1-2*03", ...],
        ...
    }

    Filtering applied:
    - Gene names longer than MAX_GENE_NAME_LENGTH are excluded (IgBLAST limit)
    - V genes without gapped sequences are excluded (required for internal_data)
    """
    genes_by_source: Dict[str, Set[str]] = defaultdict(set)
    skipped_long_names = 0
    skipped_no_gapped = 0

    manager = GermlineManager(providers=providers)

    for segment in segments:
        for chain in CHAINS:
            try:
                genes = manager.get_genes(species, segment, chain, functional_only=functional_only)
                for gene in genes:
                    if len(gene.name) > MAX_GENE_NAME_LENGTH:
                        skipped_long_names += 1
                    elif segment == "V" and len(gene.sequence) < MIN_V_GENE_LENGTH:
                        # Skip truncated V genes (missing FR1/CDR1)
                        skipped_no_gapped += 1
                    elif segment == "V" and (not gene.sequence_gapped or gene.sequence_gapped.startswith("." * 20)):
                        # Skip V genes without proper gapped sequences (too many leading gaps)
                        skipped_no_gapped += 1
                    else:
                        genes_by_source[gene.source].add(gene.name)
                if verbose and genes:
                    console.print(f"[dim]  {species}/{segment}/{chain}: {len(genes)} genes[/dim]")
            except Exception as e:
                if verbose:
                    console.print(f"[dim]  {species}/{segment}/{chain}: no data ({e})[/dim]")

    # Also get C genes (constant regions) - they use different naming
    try:
        # C genes are typically fetched differently - check if available
        for chain in CHAINS:
            c_genes = manager.get_genes(species, "C", chain, functional_only=functional_only)
            for gene in c_genes:
                if len(gene.name) <= MAX_GENE_NAME_LENGTH:
                    genes_by_source[gene.source].add(gene.name)
                else:
                    skipped_long_names += 1
    except Exception:
        pass  # C genes may not be available for all species

    if verbose and skipped_long_names > 0:
        console.print(f"[yellow]  Skipped {skipped_long_names} genes with names > {MAX_GENE_NAME_LENGTH} chars[/yellow]")
    if verbose and skipped_no_gapped > 0:
        console.print(f"[yellow]  Skipped {skipped_no_gapped} V genes without gapped sequences[/yellow]")

    # Convert sets to sorted lists
    return {source: sorted(list(genes)) for source, genes in genes_by_source.items() if genes}


def get_specific_genes(
    species: str,
    gene_names: List[str],
    providers: List[str],
    verbose: bool = False,
) -> Dict[str, List[str]]:
    """Get specific genes by name, organized by source."""
    genes_by_source: Dict[str, List[str]] = defaultdict(list)

    manager = GermlineManager(providers=providers)

    for gene_name in gene_names:
        try:
            gene = manager.get_gene_by_name(gene_name, species)
            if gene:
                genes_by_source[gene.source].append(gene.name)
            elif verbose:
                console.print(f"[yellow]Warning: Gene {gene_name} not found for {species}[/yellow]")
        except Exception as e:
            if verbose:
                console.print(f"[yellow]Warning: Could not fetch {gene_name}: {e}[/yellow]")

    return dict(genes_by_source)


def generate_standard_reference(
    name: str,
    species: str,
    providers: List[str],
    functional_only: bool = True,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Generate a standard reference (name=species, all functional genes)."""
    if verbose:
        console.print(f"\n[cyan]Generating standard reference: {name} ({species})[/cyan]")

    genes_by_source = get_all_genes_for_species(
        species=species,
        providers=providers,
        functional_only=functional_only,
        verbose=verbose,
    )

    if not genes_by_source:
        console.print(f"[yellow]Warning: No genes found for {species}[/yellow]")
        return {}

    # Build reference structure: {source: {species: [genes]}}
    reference = {}
    for source, genes in genes_by_source.items():
        if genes:
            reference[source] = {species: genes}

    return reference


def generate_chimeric_reference(
    name: str,
    components: Dict[str, Dict[str, Any]],
    functional_only: bool = True,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Generate a chimeric reference (multiple species combined)."""
    if verbose:
        console.print(f"\n[cyan]Generating chimeric reference: {name}[/cyan]")

    reference: Dict[str, Dict[str, List[str]]] = defaultdict(dict)

    for species, config in components.items():
        providers = config.get("providers", ["imgt"])
        genes_config = config.get("genes", "ALL_FUNCTIONAL")

        if verbose:
            console.print(f"  [dim]Component: {species}[/dim]")

        if genes_config == "ALL_FUNCTIONAL":
            genes_by_source = get_all_genes_for_species(
                species=species,
                providers=providers,
                functional_only=functional_only,
                verbose=verbose,
            )
        else:
            # Specific gene list
            genes_by_source = get_specific_genes(
                species=species,
                gene_names=genes_config,
                providers=providers,
                verbose=verbose,
            )

        # Merge into reference structure
        for source, genes in genes_by_source.items():
            if genes:
                reference[source][species] = genes

    return dict(reference)


def generate_reference(
    name: str,
    config: Dict[str, Any],
    functional_only: bool = True,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Generate a reference based on its configuration."""
    ref_type = config.get("type", "standard")

    if ref_type == "standard":
        reference = generate_standard_reference(
            name=name,
            species=config["species"],
            providers=config.get("providers", ["vdjbase", "ogrdb", "imgt"]),
            functional_only=functional_only,
            verbose=verbose,
        )
    elif ref_type == "chimeric":
        reference = generate_chimeric_reference(
            name=name,
            components=config["components"],
            functional_only=functional_only,
            verbose=verbose,
        )
    else:
        console.print(f"[red]Unknown reference type: {ref_type}[/red]")
        return {}

    # Deduplicate: remove same gene appearing in multiple sources for same species
    if reference:
        reference = deduplicate_reference(reference, verbose=verbose)

    return reference


def generate_all_references(
    functional_only: bool = True,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Generate complete reference.yml with all configured references."""
    all_references = {}

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Generating references...", total=len(REFERENCE_CONFIG))

        for name, config in REFERENCE_CONFIG.items():
            progress.update(task, description=f"Generating {name}...")

            reference = generate_reference(
                name=name,
                config=config,
                functional_only=functional_only,
                verbose=verbose,
            )

            if reference:
                all_references[name] = reference

            progress.advance(task)

    return all_references


def display_reference_summary(references: Dict[str, Any]) -> None:
    """Display summary of generated references."""
    table = Table(title="Generated References Summary", show_header=True)
    table.add_column("Reference", style="cyan")
    table.add_column("Source", style="green")
    table.add_column("Species", style="yellow")
    table.add_column("Genes", style="magenta", justify="right")

    for ref_name, ref_data in references.items():
        first_row = True
        for source, species_data in ref_data.items():
            for species, genes in species_data.items():
                table.add_row(
                    ref_name if first_row else "",
                    source,
                    species,
                    str(len(genes)),
                )
                first_row = False

    console.print(table)


def compare_references(
    generated: Dict[str, Any],
    existing_path: Path,
    verbose: bool = False,
) -> None:
    """Compare generated references with existing file."""
    if not existing_path.exists():
        console.print(f"[red]Comparison file not found: {existing_path}[/red]")
        return

    with open(existing_path) as f:
        existing = yaml.safe_load(f)

    console.print("\n[bold]Comparison with existing reference.yml[/bold]")

    table = Table(show_header=True)
    table.add_column("Reference", style="cyan")
    table.add_column("Existing", style="yellow", justify="right")
    table.add_column("Generated", style="green", justify="right")
    table.add_column("Diff", style="magenta", justify="right")

    for ref_name in set(list(existing.keys()) + list(generated.keys())):
        existing_count = 0
        generated_count = 0

        if ref_name in existing:
            for source in existing[ref_name]:
                for species in existing[ref_name][source]:
                    existing_count += len(existing[ref_name][source][species])

        if ref_name in generated:
            for source in generated[ref_name]:
                for species in generated[ref_name][source]:
                    generated_count += len(generated[ref_name][source][species])

        diff = generated_count - existing_count
        diff_str = f"+{diff}" if diff > 0 else str(diff)

        table.add_row(ref_name, str(existing_count), str(generated_count), diff_str)

    console.print(table)


def backup_file(file_path: Path) -> Optional[Path]:
    """Create timestamped backup of existing file."""
    if not file_path.exists():
        return None

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = file_path.with_suffix(f".bak.{timestamp}.yml")
    shutil.copy2(file_path, backup_path)
    return backup_path


def interactive_mode() -> Optional[Dict[str, Any]]:
    """Interactive mode for reference generation."""
    console.print(Panel.fit("[bold cyan]SADIE Reference Generator - Interactive Mode[/bold cyan]"))

    # List available species
    console.print("\n[cyan]Checking available species...[/cyan]")
    species_by_provider = get_available_species(["vdjbase", "ogrdb", "imgt"], verbose=True)

    # Merge all species
    all_species = set()
    for species_list in species_by_provider.values():
        all_species.update(species_list)

    console.print(f"\nAvailable species: {', '.join(sorted(all_species))}")

    # Get reference name
    name = Prompt.ask("[cyan]Reference name[/cyan]", default="custom")

    # Get species
    species = Prompt.ask("[cyan]Species[/cyan]", default="human")

    # Get providers
    providers_input = Prompt.ask(
        "[cyan]Providers (space-separated)[/cyan]",
        default="vdjbase ogrdb imgt",
    )
    providers = providers_input.split()

    # Functional only?
    functional_only = Confirm.ask("[cyan]Include only functional genes?[/cyan]", default=True)

    # Generate
    console.print(f"\n[green]Generating reference '{name}' for {species}...[/green]")

    reference = generate_standard_reference(
        name=name,
        species=species,
        providers=providers,
        functional_only=functional_only,
        verbose=True,
    )

    if not reference:
        console.print("[red]No genes found![/red]")
        return None

    # Show summary
    total_genes = sum(len(genes) for source_data in reference.values() for genes in source_data.values())
    console.print(f"\n[green]Found {total_genes} genes[/green]")

    # Confirm
    if not Confirm.ask("\n[cyan]Create reference with these genes?[/cyan]", default=True):
        return None

    return {name: reference}


def main() -> None:
    """Main entry point."""
    args = parse_arguments()

    # List species mode
    if args.list_species:
        console.print(Panel.fit("[bold cyan]Available Species by Provider[/bold cyan]"))
        species_by_provider = get_available_species(args.providers, verbose=True)

        for provider, species_list in species_by_provider.items():
            console.print(f"\n[bold]{provider.upper()}[/bold]: {len(species_list)} species")
            for sp in species_list:
                console.print(f"  - {sp}")
        return

    # Interactive mode
    if args.interactive:
        references = interactive_mode()
        if not references:
            console.print("[yellow]Reference generation cancelled[/yellow]")
            return
    # Generate all mode
    elif args.generate_all:
        console.print(Panel.fit("[bold cyan]Generating Complete reference.yml[/bold cyan]"))
        functional_only = not args.include_non_functional
        references = generate_all_references(
            functional_only=functional_only,
            verbose=args.verbose,
        )
    # Single species mode
    elif args.species:
        if not args.name:
            args.name = args.species

        console.print(f"[cyan]Generating reference for {args.species}...[/cyan]")

        functional_only = not args.include_non_functional
        reference = generate_standard_reference(
            name=args.name,
            species=args.species,
            providers=args.providers,
            functional_only=functional_only,
            verbose=args.verbose,
        )

        if not reference:
            console.print(f"[red]No genes found for {args.species}[/red]")
            return

        references = {args.name: reference}
    else:
        console.print("[red]Error: Specify --generate-all, --species, or --interactive[/red]")
        return

    # Display summary
    display_reference_summary(references)

    # Compare if requested
    if args.compare:
        compare_references(references, args.compare, verbose=args.verbose)

    # Dry run - don't write
    if args.dry_run:
        console.print("\n[yellow]DRY RUN - No files written[/yellow]")
        if args.verbose:
            console.print("\n[dim]Generated YAML:[/dim]")
            console.print(yaml.dump(references, default_flow_style=False, sort_keys=False))
        return

    # Backup existing file
    if args.output.exists() and not args.no_backup:
        backup_path = backup_file(args.output)
        if backup_path:
            console.print(f"[dim]Backed up existing file to {backup_path}[/dim]")

    # Update mode - merge with existing
    if args.update and args.output.exists():
        with open(args.output) as f:
            existing = yaml.safe_load(f) or {}
        existing.update(references)
        references = existing

    # Write output
    with open(args.output, "w") as f:
        yaml.dump(references, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

    console.print(f"\n[green]Successfully wrote {args.output}[/green]")

    # Show total gene count
    total_genes = 0
    for ref_data in references.values():
        for source_data in ref_data.values():
            for genes in source_data.values():
                total_genes += len(genes)

    console.print(f"[dim]Total: {len(references)} references, {total_genes} genes[/dim]")


if __name__ == "__main__":
    main()
