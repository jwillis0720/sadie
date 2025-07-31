#!/usr/bin/env python3
"""
Regenerate CATNAP reference files for integration tests.

This script should be run whenever the IgBLAST databases are updated
to ensure integration tests use the latest reference data.

Usage:
    python scripts/regenerate_catnap_references.py
    python scripts/regenerate_catnap_references.py --no-backup --verbose
    python scripts/regenerate_catnap_references.py --dry-run
    python scripts/regenerate_catnap_references.py --run-tests
"""

import argparse
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

try:
    import pandas as pd
    from rich import print as rprint
    from rich.console import Console
    from rich.panel import Panel
    from rich.progress import (
        BarColumn,
        Progress,
        SpinnerColumn,
        TextColumn,
        TimeElapsedColumn,
    )
    from rich.table import Table
except ImportError as e:
    print(f"Error: Missing required dependencies. Please install: {e}")
    print("Run: pip install pandas rich")
    sys.exit(1)

try:
    from sadie.airr import Airr
except ImportError:
    print("Error: SADIE is not installed or not in PYTHONPATH")
    print("Please ensure SADIE is properly installed")
    sys.exit(1)

console = Console()

# Default paths
DEFAULT_INPUT_DIR = Path("tests/data/fixtures/fasta_inputs")
DEFAULT_OUTPUT_DIR = Path("tests/data/fixtures/airr_tables")
INPUT_FILES = {
    "heavy": "catnap_nt_heavy.fasta",
    "light": "catnap_nt_light.fasta"
}
OUTPUT_FILES = {
    "heavy": "catnap_heavy_airrtable.feather",
    "light": "catnap_light_airrtable.feather"
}


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Regenerate CATNAP reference files for integration tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                    # Basic usage with defaults
  %(prog)s --no-backup        # Skip backing up existing files
  %(prog)s --dry-run          # Show what would be done without changes
  %(prog)s --run-tests        # Run integration tests after regeneration
  %(prog)s --verbose          # Show detailed output
        """
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Input directory containing FASTA files (default: {DEFAULT_INPUT_DIR})"
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory for reference files (default: {DEFAULT_OUTPUT_DIR})"
    )

    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Skip backing up existing reference files"
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes"
    )

    parser.add_argument(
        "--run-tests",
        action="store_true",
        help="Run integration tests after regeneration"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output"
    )

    return parser.parse_args()


def backup_file(file_path: Path) -> Optional[Path]:
    """Create a timestamped backup of a file if it exists."""
    if not file_path.exists():
        return None

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = file_path.with_suffix(f".bak.{timestamp}{file_path.suffix}")

    shutil.copy2(file_path, backup_path)
    return backup_path


def compare_airrtables(old_path: Path, new_df: pd.DataFrame) -> Dict[str, Any]:
    """Compare old and new AirrTable dataframes."""
    comparison = {
        "exists": False,
        "row_diff": 0,
        "column_diff": set(),
        "shape_old": None,
        "shape_new": new_df.shape,
        "corrected_changes": {"old": 0, "new": 0}
    }

    if not old_path.exists():
        return comparison

    try:
        old_df = pd.read_feather(old_path)
        comparison["exists"] = True
        comparison["shape_old"] = old_df.shape
        comparison["row_diff"] = len(new_df) - len(old_df)

        # Column differences
        old_cols = set(old_df.columns)
        new_cols = set(new_df.columns)
        comparison["column_diff"] = {
            "added": new_cols - old_cols,
            "removed": old_cols - new_cols
        }

        # Check v_germline_alignment_aa_corrected changes
        if "v_germline_alignment_aa_corrected" in old_df.columns and "v_germline_alignment_aa_corrected" in new_df.columns:
            old_corrected = old_df["v_germline_alignment_aa_corrected"].sum() if "v_germline_alignment_aa_corrected" in old_df.columns else 0
            new_corrected = new_df["v_germline_alignment_aa_corrected"].sum() if "v_germline_alignment_aa_corrected" in new_df.columns else 0
            comparison["corrected_changes"] = {
                "old": int(old_corrected),
                "new": int(new_corrected)
            }

    except Exception as e:
        comparison["error"] = str(e)

    return comparison


def regenerate_reference(
    chain_type: str,
    input_path: Path,
    output_path: Path,
    dry_run: bool = False,
    verbose: bool = False
) -> Tuple[Optional[pd.DataFrame], float]:
    """Regenerate a single reference file."""
    start_time = time.time()

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if dry_run:
        console.print(f"[yellow]DRY RUN:[/yellow] Would process {input_path}")
        return None, 0.0

    # Run Airr analysis
    airr_api = Airr("human", adaptable=True)
    result_df = airr_api.run_fasta(str(input_path))

    # Save to feather format
    result_df.to_feather(output_path)

    elapsed_time = time.time() - start_time
    return result_df, elapsed_time


def display_comparison_results(comparisons: Dict[str, Dict[str, Any]]) -> None:
    """Display comparison results in a formatted table."""
    table = Table(title="Reference File Comparison", show_header=True)
    table.add_column("Chain", style="cyan")
    table.add_column("Old Shape", style="yellow")
    table.add_column("New Shape", style="green")
    table.add_column("Row Diff", style="magenta")
    table.add_column("Corrected (Old→New)", style="red")

    for chain_type, comp in comparisons.items():
        old_shape = f"{comp['shape_old']}" if comp['shape_old'] else "N/A"
        new_shape = f"{comp['shape_new']}"
        row_diff = f"{comp['row_diff']:+d}" if comp['exists'] else "N/A"
        corrected = f"{comp['corrected_changes']['old']}→{comp['corrected_changes']['new']}"

        table.add_row(
            chain_type.capitalize(),
            old_shape,
            new_shape,
            row_diff,
            corrected
        )

    console.print(table)

    # Show column changes if any
    for chain_type, comp in comparisons.items():
        if comp.get("column_diff"):
            if comp["column_diff"]["added"]:
                console.print(f"\n[green]Added columns in {chain_type}:[/green] {', '.join(comp['column_diff']['added'])}")
            if comp["column_diff"]["removed"]:
                console.print(f"[red]Removed columns in {chain_type}:[/red] {', '.join(comp['column_diff']['removed'])}")


def run_integration_tests() -> bool:
    """Run the CATNAP integration test."""
    console.print("\n[cyan]Running integration tests...[/cyan]")

    cmd = [
        "pytest", "-xvs",
        "tests/integration/airr/test_airr_intergration.py::test_catnap_integration"
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            console.print("[green]✓ Integration tests passed![/green]")
            return True
        else:
            console.print("[red]✗ Integration tests failed![/red]")
            console.print(f"[dim]{result.stdout}[/dim]")
            if result.stderr:
                console.print(f"[red]{result.stderr}[/red]")
            return False

    except Exception as e:
        console.print(f"[red]Error running tests: {e}[/red]")
        return False


def main():
    """Main execution function."""
    args = parse_arguments()

    # Display header
    console.print(Panel.fit(
        "[bold blue]CATNAP Reference Regeneration Tool[/bold blue]\n"
        "Regenerates reference files for integration tests",
        border_style="blue"
    ))

    # Validate directories
    if not args.input_dir.exists():
        console.print(f"[red]Error: Input directory not found: {args.input_dir}[/red]")
        sys.exit(1)

    if not args.dry_run and not args.output_dir.exists():
        console.print(f"[yellow]Creating output directory: {args.output_dir}[/yellow]")
        args.output_dir.mkdir(parents=True, exist_ok=True)

    # Prepare file paths
    files_to_process = {}
    for chain_type, filename in INPUT_FILES.items():
        input_path = args.input_dir / filename
        output_path = args.output_dir / OUTPUT_FILES[chain_type]
        files_to_process[chain_type] = {
            "input": input_path,
            "output": output_path
        }

    # Backup existing files
    backups = {}
    if not args.no_backup and not args.dry_run:
        console.print("\n[cyan]Creating backups...[/cyan]")
        for chain_type, paths in files_to_process.items():
            backup_path = backup_file(paths["output"])
            if backup_path:
                backups[chain_type] = backup_path
                console.print(f"  ✓ Backed up {paths['output'].name} → {backup_path.name}")

    # Process files
    results = {}
    comparisons = {}

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True
    ) as progress:

        for chain_type, paths in files_to_process.items():
            task_desc = f"Processing {chain_type} chain"

            if args.dry_run:
                console.print(f"[yellow]DRY RUN:[/yellow] Would process {chain_type} chain")
                # Still do comparison for dry run
                dummy_df = pd.DataFrame()
                comparisons[chain_type] = compare_airrtables(paths["output"], dummy_df)
            else:
                task = progress.add_task(task_desc, total=None)

                try:
                    # Compare before regeneration
                    result_df, elapsed = regenerate_reference(
                        chain_type,
                        paths["input"],
                        paths["output"],
                        dry_run=args.dry_run,
                        verbose=args.verbose
                    )

                    if result_df is not None:
                        comparisons[chain_type] = compare_airrtables(
                            backups.get(chain_type, paths["output"]),
                            result_df
                        )
                        results[chain_type] = {
                            "success": True,
                            "elapsed": elapsed,
                            "shape": result_df.shape
                        }

                    progress.update(task, completed=True)

                except Exception as e:
                    console.print(f"[red]Error processing {chain_type}: {e}[/red]")
                    results[chain_type] = {"success": False, "error": str(e)}
                    progress.update(task, completed=True)

    # Display results
    if not args.dry_run:
        console.print("\n[green]Processing complete![/green]")
        summary_table = Table(title="Processing Summary", show_header=True)
        summary_table.add_column("Chain", style="cyan")
        summary_table.add_column("Status", style="green")
        summary_table.add_column("Time", style="yellow")
        summary_table.add_column("Output Shape", style="magenta")

        for chain_type, result in results.items():
            if result["success"]:
                status = "✓ Success"
                time_str = f"{result['elapsed']:.2f}s"
                shape_str = str(result["shape"])
            else:
                status = "✗ Failed"
                time_str = "N/A"
                shape_str = "N/A"

            summary_table.add_row(chain_type.capitalize(), status, time_str, shape_str)

        console.print(summary_table)

    # Display comparison results
    if comparisons:
        console.print()
        display_comparison_results(comparisons)

    # Run tests if requested
    if args.run_tests and not args.dry_run:
        test_passed = run_integration_tests()
        if not test_passed:
            sys.exit(1)

    # Final message
    if args.dry_run:
        console.print("\n[yellow]DRY RUN COMPLETE - No files were modified[/yellow]")
    else:
        console.print("\n[green]✓ Reference files regenerated successfully![/green]")
        console.print("[dim]Run integration tests with: pytest -xvs tests/integration/airr/test_airr_intergration.py::test_catnap_integration[/dim]")


if __name__ == "__main__":
    main()
