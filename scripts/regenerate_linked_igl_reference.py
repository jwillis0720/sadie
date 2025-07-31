#!/usr/bin/env python3
"""
Regenerate linked IGL reference file for test_hard_igl_seqs_linked.

This script regenerates the expected output file (bum_link_solution.feather) used in
the test_hard_igl_seqs_linked test case. It should be run whenever the IgBLAST
databases or termini buffer logic is updated.

Usage:
    python scripts/regenerate_linked_igl_reference.py
    python scripts/regenerate_linked_igl_reference.py --no-backup --verbose
    python scripts/regenerate_linked_igl_reference.py --dry-run
    python scripts/regenerate_linked_igl_reference.py --run-test
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
    from sadie.airr import LinkedAirrTable
    from sadie.airr import methods as airr_methods
except ImportError:
    print("Error: SADIE is not installed or not in PYTHONPATH")
    print("Please ensure SADIE is properly installed")
    sys.exit(1)

console = Console()

# Default paths
INPUT_FILE = Path("tests/data/fixtures/airr_tables/bum_link_input.feather")
OUTPUT_FILE = Path("tests/data/fixtures/airr_tables/bum_link_solution.feather")


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Regenerate linked IGL reference file for test_hard_igl_seqs_linked",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                    # Basic usage with defaults
  %(prog)s --no-backup        # Skip backing up existing files
  %(prog)s --dry-run          # Show what would be done without changes
  %(prog)s --run-test         # Run test after regeneration
  %(prog)s --verbose          # Show detailed output
        """
    )

    parser.add_argument(
        "--input-file",
        type=Path,
        default=INPUT_FILE,
        help=f"Input feather file (default: {INPUT_FILE})"
    )

    parser.add_argument(
        "--output-file",
        type=Path,
        default=OUTPUT_FILE,
        help=f"Output feather file (default: {OUTPUT_FILE})"
    )

    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Skip backing up existing reference file"
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes"
    )

    parser.add_argument(
        "--run-test",
        action="store_true",
        help="Run test_hard_igl_seqs_linked after regeneration"
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


def compare_airrtables(old_path: Path, new_df: pd.DataFrame, verbose: bool = False) -> Dict[str, Any]:
    """Compare old and new LinkedAirrTable dataframes, focusing on sequence changes."""
    comparison = {
        "exists": False,
        "row_diff": 0,
        "column_diff": set(),
        "shape_old": None,
        "shape_new": new_df.shape,
        "sequence_changes": {
            "heavy": {"count": 0, "percentage": 0.0, "indices": []},
            "light": {"count": 0, "percentage": 0.0, "indices": []},
            "total": {"count": 0, "percentage": 0.0}
        }
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

        # Check sequence column changes for both heavy and light chains
        if len(old_df) == len(new_df):
            # Check for heavy chain sequence changes
            if "sequence_heavy" in old_df.columns and "sequence_heavy" in new_df.columns:
                heavy_diff = old_df["sequence_heavy"] != new_df["sequence_heavy"]
                heavy_changed = heavy_diff[heavy_diff].index.tolist()

                comparison["sequence_changes"]["heavy"] = {
                    "count": len(heavy_changed),
                    "percentage": (len(heavy_changed) / len(old_df)) * 100 if len(old_df) > 0 else 0,
                    "indices": heavy_changed[:10] if verbose else []
                }

            # Check for light chain sequence changes
            if "sequence_light" in old_df.columns and "sequence_light" in new_df.columns:
                light_diff = old_df["sequence_light"] != new_df["sequence_light"]
                light_changed = light_diff[light_diff].index.tolist()

                comparison["sequence_changes"]["light"] = {
                    "count": len(light_changed),
                    "percentage": (len(light_changed) / len(old_df)) * 100 if len(old_df) > 0 else 0,
                    "indices": light_changed[:10] if verbose else []
                }

            # Calculate total changes
            total_changed = len(set(heavy_changed + light_changed)) if 'heavy_changed' in locals() and 'light_changed' in locals() else 0
            comparison["sequence_changes"]["total"] = {
                "count": total_changed,
                "percentage": (total_changed / len(old_df)) * 100 if len(old_df) > 0 else 0
            }

    except Exception as e:
        comparison["error"] = str(e)

    return comparison


def regenerate_linked_igl_reference(
    input_path: Path,
    output_path: Path,
    dry_run: bool = False,
    verbose: bool = False
) -> Tuple[Optional[pd.DataFrame], float]:
    """Regenerate the linked IGL reference file using the same process as the test."""
    start_time = time.time()

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if dry_run:
        console.print(f"[yellow]DRY RUN:[/yellow] Would process {input_path}")
        # Still load and process to show what would happen
        lat = LinkedAirrTable(pd.read_feather(input_path))
        return lat, 0.0

    # Load input data
    if verbose:
        console.print(f"Loading input file: {input_path}")
    lat = LinkedAirrTable(pd.read_feather(input_path))

    # Apply the same transformations as the test
    if verbose:
        console.print("Running termini buffers...")
    out_lat = airr_methods.run_termini_buffers(lat)

    if verbose:
        console.print("Running IGL assignment...")
    igl_df = airr_methods.run_igl_assignment(out_lat)

    # Save output
    if not dry_run:
        igl_df.to_feather(output_path)
        if verbose:
            console.print(f"Saved output to: {output_path}")

    elapsed_time = time.time() - start_time
    return igl_df, elapsed_time


def display_comparison_results(comparison: Dict[str, Any]) -> None:
    """Display comparison results in a formatted table."""
    table = Table(title="Reference File Comparison", show_header=True)
    table.add_column("Metric", style="cyan")
    table.add_column("Old", style="yellow")
    table.add_column("New", style="green")
    table.add_column("Difference", style="magenta")

    # Shape comparison
    old_shape = f"{comparison['shape_old']}" if comparison['shape_old'] else "N/A"
    new_shape = f"{comparison['shape_new']}"
    table.add_row("Shape", old_shape, new_shape, f"{comparison['row_diff']:+d} rows" if comparison['exists'] else "N/A")

    # Sequence changes for heavy chain
    if comparison['exists'] and comparison['sequence_changes']['heavy']['count'] > 0:
        heavy_changes = comparison['sequence_changes']['heavy']
        table.add_row(
            "Heavy Sequences Modified",
            "-",
            f"{heavy_changes['count']}",
            f"{heavy_changes['percentage']:.2f}%"
        )

    # Sequence changes for light chain
    if comparison['exists'] and comparison['sequence_changes']['light']['count'] > 0:
        light_changes = comparison['sequence_changes']['light']
        table.add_row(
            "Light Sequences Modified",
            "-",
            f"{light_changes['count']}",
            f"{light_changes['percentage']:.2f}%"
        )

    # Total sequence changes
    if comparison['exists'] and comparison['sequence_changes']['total']['count'] > 0:
        total_changes = comparison['sequence_changes']['total']
        table.add_row(
            "Total Sequences Modified",
            "-",
            f"{total_changes['count']}",
            f"{total_changes['percentage']:.2f}%"
        )

    console.print(table)

    # Show column changes if any
    if comparison.get("column_diff"):
        if comparison["column_diff"]["added"]:
            console.print(f"\n[green]Added columns:[/green] {', '.join(comparison['column_diff']['added'])}")
        if comparison["column_diff"]["removed"]:
            console.print(f"[red]Removed columns:[/red] {', '.join(comparison['column_diff']['removed'])}")

    # Show modified sequence indices if verbose
    if comparison.get("sequence_changes", {}).get("heavy", {}).get("indices"):
        indices = comparison["sequence_changes"]["heavy"]["indices"]
        console.print(f"\n[yellow]Modified heavy chain sequence indices:[/yellow] {', '.join(map(str, indices))}")

    if comparison.get("sequence_changes", {}).get("light", {}).get("indices"):
        indices = comparison["sequence_changes"]["light"]["indices"]
        console.print(f"[yellow]Modified light chain sequence indices:[/yellow] {', '.join(map(str, indices))}")


def run_test() -> bool:
    """Run the test_hard_igl_seqs_linked test."""
    console.print("\n[cyan]Running test_hard_igl_seqs_linked...[/cyan]")

    cmd = [
        "pytest", "-xsv",
        "tests/unit/airr/test_airr.py::test_hard_igl_seqs_linked"
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            console.print("[green]✓ Test passed![/green]")
            return True
        else:
            console.print("[red]✗ Test failed![/red]")
            # Show only the assertion error, not the full output
            if "AssertionError" in result.stdout:
                lines = result.stdout.split('\n')
                for i, line in enumerate(lines):
                    if "AssertionError" in line:
                        # Show a few lines around the assertion
                        start = max(0, i - 2)
                        end = min(len(lines), i + 5)
                        console.print("[dim]" + '\n'.join(lines[start:end]) + "[/dim]")
                        break
            return False

    except Exception as e:
        console.print(f"[red]Error running test: {e}[/red]")
        return False


def main():
    """Main execution function."""
    args = parse_arguments()

    # Display header
    console.print(Panel.fit(
        "[bold blue]Linked IGL Reference Regeneration Tool[/bold blue]\n"
        "Regenerates reference file for test_hard_igl_seqs_linked",
        border_style="blue"
    ))

    # Validate input file
    if not args.input_file.exists():
        console.print(f"[red]Error: Input file not found: {args.input_file}[/red]")
        sys.exit(1)

    # Create output directory if needed
    if not args.dry_run and not args.output_file.parent.exists():
        console.print(f"[yellow]Creating output directory: {args.output_file.parent}[/yellow]")
        args.output_file.parent.mkdir(parents=True, exist_ok=True)

    # Backup existing file
    backup_path = None
    if not args.no_backup and not args.dry_run and args.output_file.exists():
        console.print("\n[cyan]Creating backup...[/cyan]")
        backup_path = backup_file(args.output_file)
        if backup_path:
            console.print(f"  ✓ Backed up {args.output_file.name} → {backup_path.name}")

    # Process file
    result_df = None
    comparison = {}

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True
    ) as progress:

        task = progress.add_task("Processing linked IGL reference file", total=None)

        try:
            # Regenerate reference
            result_df, elapsed = regenerate_linked_igl_reference(
                args.input_file,
                args.output_file,
                dry_run=args.dry_run,
                verbose=args.verbose
            )

            if result_df is not None:
                # Compare with old version
                comparison = compare_airrtables(
                    backup_path if backup_path else args.output_file,
                    result_df,
                    verbose=args.verbose
                )

            progress.update(task, completed=True)
            success = True

        except Exception as e:
            console.print(f"[red]Error processing file: {e}[/red]")
            if args.verbose:
                import traceback
                console.print(f"[dim]{traceback.format_exc()}[/dim]")
            progress.update(task, completed=True)
            success = False
            sys.exit(1)

    # Display results
    if success and not args.dry_run:
        console.print(f"\n[green]Processing complete![/green] Time: {elapsed:.2f}s")

        summary_table = Table(title="Processing Summary", show_header=True)
        summary_table.add_column("File", style="cyan")
        summary_table.add_column("Status", style="green")
        summary_table.add_column("Time", style="yellow")
        summary_table.add_column("Output Shape", style="magenta")

        summary_table.add_row(
            "linked_igl",
            "✓ Success",
            f"{elapsed:.2f}s",
            str(result_df.shape) if result_df is not None else "N/A"
        )

        console.print(summary_table)

    # Display comparison results
    if comparison:
        console.print()
        display_comparison_results(comparison)

    # Run test if requested
    if args.run_test and not args.dry_run:
        test_passed = run_test()
        if not test_passed:
            console.print("\n[yellow]Note: Test may still fail if other aspects changed.[/yellow]")
            console.print("[yellow]Review the changes and update test expectations if needed.[/yellow]")

    # Final message
    if args.dry_run:
        console.print("\n[yellow]DRY RUN COMPLETE - No files were modified[/yellow]")
    else:
        console.print("\n[green]✓ Reference file regenerated successfully![/green]")
        if not args.run_test:
            console.print("[dim]Run test with: pytest -xsv tests/unit/airr/test_airr.py::test_hard_igl_seqs_linked[/dim]")


if __name__ == "__main__":
    main()
