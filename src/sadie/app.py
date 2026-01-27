"""This is our main entry point"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, List, Optional, Union

import click

try:
    from importlib.metadata import version
except ImportError:
    # Fallback for Python < 3.8
    from importlib_metadata import version

# airr
from sadie.airr import Airr
from sadie.airr.methods import run_igl_assignment, run_mutational_analysis

# reference
from sadie.reference.reference import References

# Renumbering
from sadie.renumbering import Renumbering

# utility
from sadie.utility import SadieInputDir, SadieInputFile, SadieOutput
from sadie.utility.util import get_project_root, getVerbosityLevel

try:
    __version__ = version("sadie-antibody")
except Exception:
    try:
        __version__ = version("sadie-test")
    except Exception:
        __version__ = "unknown"


@click.group()
@click.version_option(version=__version__)
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=3,
    help="Vebosity level, ex. -vvvvv for debug level logging",
)
def sadie(verbose: int) -> None:
    numeric_level = getVerbosityLevel(verbose)
    logging.basicConfig(level=numeric_level)


@sadie.command("airr")
@click.option(
    "--name", "-n", type=click.Choice(Airr.get_available_species()), help="Species to annotate", default="human"
)
@click.option("--skip-igl", is_flag=True, help="Skip the igl assignment")
@click.option("--skip-mutation", is_flag=True, help="Skip the somewhat time instansive mutational analysis")
@click.argument(
    "input_path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True),
)
@click.argument(
    "output_path", required=False, type=click.Path(file_okay=True, dir_okay=False, writable=True, resolve_path=True)
)
def airr(name: str, skip_igl: bool, skip_mutation: bool, input_path: Path, output_path: Union[None, Path, str]) -> None:
    """Run the AIRR annotation pipeline from the command line on a single file or a directory of abi files.


    if you give a directory of abi files, it will combine all the records into a single file and annotate that.

    if you do not provide an output, the default is airr .tsv file
    """
    airr = Airr(name)
    input_path = Path(input_path)

    # if we gave a directory of abi files
    if input_path.is_dir():
        input_object = SadieInputDir(input_path)
        # get all the records inside of the directory
        _records = input_object.get_combined_seq_records()

    # else we listed a single file
    else:
        input_object = SadieInputFile(input_path)
        _records = input_object.get_seq_records()

    # they did not give an output path, defaults to airr
    if output_path is None:
        output_path = input_path.parent / Path(input_path.stem + ".tsv.gz")
    output_object = SadieOutput(output_path)
    airr_table = airr.run_records(_records)

    # run the mutational analysis
    if not skip_mutation:
        airr_table = run_mutational_analysis(airr_table, "kabat")

    # run the iGL assignment
    if not skip_igl:
        airr_table = run_igl_assignment(airr_table)

    # output handler
    if output_object.output_format in ["airr", "tsv"]:
        airr_table.to_airr(str(output_object.output_path))
    elif output_object.output_format == "csv":
        airr_table.to_csv(output_object.output_path)
    elif output_object.output_format == "gb":
        airr_table.to_genbank(str(output_object.output_path))
    elif output_object.output_format == "feather":
        airr_table.to_feather(output_object.output_path)
    else:
        raise ValueError(f"Output format {output_object.output_format} not recognized")
    click.echo(f"File written to: {str(output_object.output_path)}")


def _validate_numbering_objects(ctx: click.Context, param: Any, value: str) -> List[str]:
    """Private method for click context to evaluate comma seperated lists and make sure each field is okay"""
    columns = [c.strip() for c in value.split(",")]
    param_name = param.human_readable_name
    if param_name == "allowed_species":
        avail_columns = Renumbering.get_allowed_species()
    elif param_name == "allowed_chains":
        avail_columns = Renumbering.get_allowed_chains()
    else:
        raise ValueError(f"{param.human_readable_name} not recognized as a valid param")
    for c in columns:
        if c not in avail_columns:
            raise click.BadOptionUsage(param, f"{c} is not available. Only have, {avail_columns}")
    return columns


@sadie.command("renumbering")
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=5,
    help="Vebosity level, ex. -vvvvv for debug level logging",
)
@click.option(
    "--query",
    "-q",
    required=False,
    type=click.Path(exists=True, file_okay=True, readable=True, resolve_path=True),
    help="""The input file can be compressed or uncompressed file of fasta""",
)
@click.option(
    "--seq",
    "-i",
    required=False,
    type=str,
    help="""The input seq""",
)
@click.option(
    "--scheme",
    "-s",
    is_flag=False,
    default="imgt",
    type=click.Choice(Renumbering.get_available_numbering_schemes()),
    show_default=True,
    help="The numbering scheme to use.",
)
@click.option(
    "--region",
    "-r",
    is_flag=False,
    default="imgt",
    type=click.Choice(Renumbering.get_available_region_definitions()),
    show_default=True,
    help="The framework and cdr defition to use",
)
@click.option(
    "--allowed-species",
    "-a",
    is_flag=False,
    default=",".join(Renumbering.get_allowed_species()),
    show_default=True,
    callback=_validate_numbering_objects,
    help="A comma seperated list of species to align against",
)
@click.option(
    "--allowed-chains",
    "-c",
    is_flag=False,
    default=",".join(Renumbering.get_allowed_chains()),
    show_default=True,
    callback=_validate_numbering_objects,
    help="A comma seperated list of species to align against",
)
@click.option(
    "--out",
    "-o",
    type=click.Path(writable=True, resolve_path=True),
    help="""The output file, type is inferred from extensions""",
)
@click.option(
    "--compress",
    "-z",
    type=click.Choice(["gz", "bz2"]),
    help="opitonal file compression on output",
)
@click.option(
    "--file-format",
    "-f",
    type=click.Choice(["json", "csv", "feather"]),
    help="output file type format",
    default="csv",
)
def renumbering(
    verbose: bool,
    query: Path,
    seq: str,
    scheme: str,
    region: str,
    allowed_species: List[str],
    allowed_chains: List[str],
    out: Path,
    compress: str,
    file_format: str,
) -> None:
    numeric_level = getVerbosityLevel(verbose)
    logging.basicConfig(level=numeric_level)
    logger = logging.getLogger("NUMBERING")

    # No reason to use click echo over print except to show e can
    click.echo(f"Logging with level=>{logging.getLevelName(logger.getEffectiveLevel())}")
    logger.info(f"Running Renumbering on renumbering: {query}")
    logger.info(f"Allowed-species {allowed_species}")
    logger.info(f"Allowed-chains: {allowed_chains}")
    logger.info(f"Numbering: {scheme}")
    logger.info(f"Region Def: {region}")

    # setup object
    renumbering_api = Renumbering(
        scheme=scheme,
        region_assign=region,
        allowed_chain=allowed_chains,
        allowed_species=allowed_species,
        use_numbering_hmms=True,
    )

    # # run file on query
    if query is not None:
        renumbering_results = renumbering_api.run_file(query)
    else:
        renumbering_results = renumbering_api.run_single(seq_id="0", seq=seq)

    # deal with output
    # if no output file, name after input
    if out:
        out = Path(out)
        segment_out = str(out.stem) + "_renumbering_results" + str(out.suffix)
        align_out = str(out.stem) + "_numbering_alignment" + str(out.suffix)
    else:
        input_path = Path(query) if query else Path("query")
        if compress and file_format.lower() != "feather":
            # feather can't be compressed
            compress = "." + compress
        else:
            compress = ""
        segment_out = input_path.stem + f"_numbering_segment.{file_format.lower()}{compress}"
        align_out = input_path.stem + f"_numbering_alignment.{file_format.lower()}{compress}"

    # deal with file format
    # csv
    if file_format.lower() == "csv":
        renumbering_results.to_csv(segment_out)
        renumbering_results.get_alignment_table().to_csv(align_out)

    # json
    elif file_format.lower() == "json":
        renumbering_results.to_json(segment_out, orient="records")
        renumbering_results.get_alignment_table().to_json(align_out, orient="records")

    # feather
    elif file_format.lower() == "feather":
        renumbering_results.reset_index(drop=True).to_feather(segment_out)
        renumbering_results.get_alignment_table().reset_index().to_feather(align_out)

    # shouldn't get here, but if they specify a file format that is not recognized using an invoke method
    else:
        raise ValueError(f"{file_format} not recognized")
    logger.info(f"Output: {segment_out}")
    logger.info(f"Output: {align_out}")


@sadie.group()
def reference() -> None:
    pass


@reference.command("make")
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=4,
    help="Vebosity level, ex. -vvvvv for debug level logging",
)
@click.option(
    "--outpath",
    "-o",
    help="Output path to generate blast_db, internal_data and auxillary files",
    type=click.Path(resolve_path=True, dir_okay=True, writable=True),
    default=os.path.join(get_project_root(), "airr/data/germlines"),
    show_default=True,
)
@click.option(
    "--reference",
    "-d",
    help="Path to reference.yml",
    type=click.Path(resolve_path=True, dir_okay=True, exists=True),
    default=os.path.join(os.path.dirname(__file__), "reference/data/reference.yml"),
    show_default=True,
)
def make_igblast_reference(verbose: int, outpath: Path, reference: Path) -> None:
    """make the igblast reference files

    DEPRECATED: This command is deprecated. Use 'sadie reference build' instead.

    This script will make the imgt reference files used by igblast or airr, including internal data, the blast
    the blast database, and the auxillary files. It uses the reference.yml to configure select genes and species.
    If you update the reference.yml file, run this again.
    """
    # Show deprecation warning
    import warnings

    warnings.warn(
        "The 'sadie reference make' command is deprecated and will be removed in a future release.\n"
        "Please use 'sadie reference build' instead:\n"
        "  sadie reference build reference.yml -o /path/to/output",
        DeprecationWarning,
        stacklevel=2,
    )

    click.echo("⚠️  WARNING: 'sadie reference make' is deprecated. Use 'sadie reference build' instead.\n")

    # Set the root logger in the console script
    # Get back a numeric level associated with number of clicks
    numeric_level = getVerbosityLevel(verbose)
    logging.basicConfig(level=numeric_level)
    click.echo(f"reference path {reference}")

    if not outpath:
        outpath = Path(__file__).parent.joinpath("reference/data/germlines").resolve()
        click.echo(f"No outpath specified, using {outpath}")

    # make reference
    click.echo(f"Getting G3 genes from {reference}")

    # read in yaml file with all statric reference data
    reference_object = References.from_yaml(reference)
    germline_path = reference_object.make_airr_database(outpath)
    click.echo(f"Wrote germline database to {germline_path}")
    click.echo("Done!")


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
    "--use-germlines/--no-use-germlines",
    default=True,
    help="Use local germlines module instead of G3 API (default: germlines module)",
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

        # Progress: Building databases...
        click.echo("Building databases...")
        germline_path = reference_object.make_airr_database(output_path)

        # Progress: Complete
        click.echo("Complete")
        click.echo(f"Database written to: {germline_path}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)


@reference.command("generate")
@click.option(
    "--output",
    "-o",
    type=click.Path(resolve_path=True),
    default="reference.yml",
    help="Output YAML file path (default: reference.yml)",
    show_default=True,
)
@click.option(
    "--generate-all",
    is_flag=True,
    help="Generate complete reference.yml with all configured references",
)
@click.option(
    "--name",
    help="Reference name (e.g., human, mouse, clk)",
)
@click.option(
    "--species",
    help="Species to query (e.g., human, mouse, macaque)",
)
@click.option(
    "--providers",
    multiple=True,
    type=click.Choice(["vdjbase", "ogrdb", "imgt", "custom"]),
    default=["vdjbase", "ogrdb", "imgt"],
    help="Germline providers to include (default: vdjbase ogrdb imgt)",
)
@click.option(
    "--functional-only/--include-non-functional",
    default=True,
    help="Include only functional genes (default: functional-only)",
)
@click.option(
    "--interactive",
    "-i",
    is_flag=True,
    help="Interactive mode for gene selection",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show what would be done without making changes",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Show detailed output",
)
@click.option(
    "--list-species",
    is_flag=True,
    help="List available species from germlines module",
)
def generate_reference(
    output: str,
    generate_all: bool,
    name: Optional[str],
    species: Optional[str],
    providers: tuple,
    functional_only: bool,
    interactive: bool,
    dry_run: bool,
    verbose: bool,
    list_species: bool,
) -> None:
    """Generate reference.yml file using SADIE germlines module.

    This command uses the local germlines database (IMGT, OGRDB, VDJbase) to generate
    reference.yml files used by 'sadie reference build' command.

    \b
    Examples:
        sadie reference generate --generate-all
        sadie reference generate --species human --output human.yml
        sadie reference generate --interactive
        sadie reference generate --list-species
    """
    import sys
    from pathlib import Path

    # Import the script's main functionality
    script_path = Path(__file__).parent.parent.parent / "scripts" / "generate_reference_yaml_germlines.py"

    if not script_path.exists():
        click.echo("Error: generate_reference_yaml_germlines.py script not found", err=True)
        sys.exit(1)

    # Build command arguments
    cmd = [sys.executable, str(script_path)]

    if list_species:
        cmd.append("--list-species")
    elif interactive:
        cmd.append("--interactive")
    elif generate_all:
        cmd.append("--generate-all")
    elif species:
        cmd.extend(["--species", species])
        if name:
            cmd.extend(["--name", name])
    else:
        click.echo("Error: Specify --generate-all, --species, --interactive, or --list-species", err=True)
        sys.exit(1)

    cmd.extend(["--output", output])

    if providers:
        cmd.append("--providers")
        cmd.extend(providers)

    if not functional_only:
        cmd.append("--include-non-functional")

    if dry_run:
        cmd.append("--dry-run")

    if verbose:
        cmd.append("--verbose")

    # Execute the script
    result = subprocess.run(cmd, capture_output=False, text=True)
    sys.exit(result.returncode)


@sadie.command()
@click.option(
    "-o",
    "--outpath",
    help="Output path for all generated files",
    type=click.Path(resolve_path=True, dir_okay=True, exists=False, file_okay=False),
    required=False,
)
@click.option(
    "-r",
    "--reference",
    help="Path to reference.yml file",
    type=click.Path(resolve_path=True, exists=False),
    default=os.path.join(os.path.dirname(__file__), "reference/data/reference.yml"),
    show_default=True,
)
@click.option(
    "--species",
    "-s",
    help="Comma-separated list of species to build databases for (default: all)",
    type=str,
    default="all",
    show_default=True,
)
@click.option("--regenerate-catnap", is_flag=True, help="Regenerate CATNAP references from FASTA files", default=False)
@click.option("--regenerate-igl", is_flag=True, help="Regenerate IGL references from FASTA files", default=False)
@click.option(
    "--generate-reference",
    is_flag=True,
    help="Generate reference.yml file interactively or regenerate if missing",
    default=False,
)
@click.option("--skip-blast", is_flag=True, help="Skip BLAST database generation", default=False)
@click.option("--skip-aux", is_flag=True, help="Skip auxiliary file generation", default=False)
@click.option("--skip-internal", is_flag=True, help="Skip internal annotation file generation", default=False)
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=3,
    help="Vebosity level, ex. -vvvvv for debug level logging",
)
def make_all(
    verbose: int,
    outpath: Optional[Path],
    reference: Path,
    species: str,
    regenerate_catnap: bool,
    regenerate_igl: bool,
    generate_reference: bool,
    skip_blast: bool,
    skip_aux: bool,
    skip_internal: bool,
) -> None:
    """Comprehensive database generation for SADIE AIRR analysis

    DEPRECATED: This command is deprecated and will be removed in a future release.
    Use the new workflow instead:
      1. sadie germlines populate              # Download germline data
      2. sadie reference generate --generate-all  # Generate reference.yml
      3. sadie reference build reference.yml -o ./db  # Build database

    This command performs all necessary steps to set up a complete AIRR analysis environment:

    1. Generates or updates reference.yml configuration (if needed)
    2. Downloads IMGT germline sequences (if needed)
    3. Generates IgBLAST reference databases
    4. Creates auxiliary files for IgBLAST
    5. Generates internal annotation files (ndm.imgt)
    6. Optionally regenerates CATNAP and IGL references
    7. Creates reference index files

    The process creates the following directory structure:

    \b
    outpath/
    ├── reference.yml            # Reference configuration (if generated)
    ├── Ig/
    │   ├── blastdb/         # BLAST databases
    │   │   ├── human/
    │   │   ├── mouse/
    │   │   └── ...
    │   └── internal_data/   # Internal annotation files
    │       ├── human/
    │       ├── mouse/
    │       └── ...
    ├── aux_db/              # Auxiliary files
    │   └── imgt/
    │       ├── human_gl.aux
    │       ├── mouse_gl.aux
    │       └── ...
    └── .references_dataframe.csv.gz  # Reference index
    """
    # Show deprecation warning
    import warnings

    warnings.warn(
        "The 'sadie make-all' command is deprecated and will be removed in a future release.\n"
        "Please use the new workflow:\n"
        "  1. sadie germlines populate              # Download germline data\n"
        "  2. sadie reference generate --generate-all  # Generate reference.yml\n"
        "  3. sadie reference build reference.yml -o ./db  # Build database",
        DeprecationWarning,
        stacklevel=2,
    )

    click.echo("⚠️  WARNING: 'sadie make-all' is deprecated. See documentation for new workflow.\n")

    # Set logging level
    numeric_level = getVerbosityLevel(verbose)
    logging.basicConfig(level=numeric_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # Default output path
    if not outpath:
        outpath = Path(__file__).parent.joinpath("airr/data/germlines").resolve()
        click.echo(f"No outpath specified, using {outpath}")

    # Ensure output directory exists
    Path(outpath).mkdir(parents=True, exist_ok=True)

    # Step 1: Generate reference.yml if requested or missing
    if generate_reference or not Path(reference).exists():
        click.echo("\n" + "=" * 60)
        click.echo("STEP 1: Generating reference.yml configuration")
        click.echo("=" * 60)

        if not Path(reference).exists():
            click.echo(f"Reference file not found at {reference}")

        script_path = Path(__file__).parent.parent.parent / "scripts" / "generate_reference_yaml.py"
        if script_path.exists():
            # Determine if we should use interactive mode
            cmd = [sys.executable, str(script_path)]

            if generate_reference and Path(reference).exists():
                # User explicitly requested regeneration
                click.echo("Regenerating reference.yml file...")
                cmd.extend(["--update", "--output", str(reference)])
            elif not Path(reference).exists():
                # File doesn't exist, create it
                click.echo("Creating new reference.yml file...")
                cmd.extend(["--interactive", "--output", str(reference)])

            if species != "all":
                # Use specified species
                for sp in species.split(","):
                    cmd.extend(["--species", sp.strip(), "--all-functional"])

            result = subprocess.run(cmd, capture_output=False, text=True)
            if result.returncode != 0:
                click.echo("Error: Reference generation failed", err=True)
                raise click.ClickException("Failed to generate reference.yml")
            else:
                click.echo("Reference.yml generated successfully")
        else:
            click.echo(f"Warning: Reference generation script not found at {script_path}", err=True)
            if not Path(reference).exists():
                raise click.ClickException(f"Reference file not found and cannot generate: {reference}")

    # Step 2: Load reference configuration
    click.echo("\n" + "=" * 60)
    click.echo("STEP 2: Loading reference configuration")
    click.echo("=" * 60)
    reference_object = References.from_yaml(reference)

    # Step 3: Regenerate CATNAP references if requested
    if regenerate_catnap:
        click.echo("\n" + "=" * 60)
        click.echo("STEP 3: Regenerating CATNAP references")
        click.echo("=" * 60)

        script_path = Path(__file__).parent.parent.parent / "scripts" / "regenerate_catnap_references.py"
        if script_path.exists():
            cmd = [sys.executable, str(script_path), "--verbose"]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                click.echo(f"Warning: CATNAP regeneration failed: {result.stderr}", err=True)
            else:
                click.echo("CATNAP references regenerated successfully")
        else:
            click.echo(f"Warning: CATNAP regeneration script not found at {script_path}", err=True)

    # Step 4: Regenerate IGL references if requested
    if regenerate_igl:
        click.echo("\n" + "=" * 60)
        click.echo("STEP 4: Regenerating IGL references")
        click.echo("=" * 60)
        script_path = Path(__file__).parent.parent.parent / "scripts" / "regenerate_igl_reference.py"
        if script_path.exists():
            cmd = [sys.executable, str(script_path), "--verbose"]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                click.echo(f"Warning: IGL regeneration failed: {result.stderr}", err=True)
            else:
                click.echo("IGL references regenerated successfully")
        else:
            click.echo(f"Warning: IGL regeneration script not found at {script_path}", err=True)

    # Step 5: Generate all database files
    click.echo("\n" + "=" * 60)
    click.echo("STEP 5: Generating database files")
    click.echo("=" * 60)

    try:
        # The make_airr_database method handles all three components
        reference_object.make_airr_database(Path(outpath) if outpath else Path.cwd())

        click.echo("\n✓ Database generation completed successfully!")
        click.echo(f"  - BLAST databases: {outpath}/Ig/blastdb/")
        click.echo(f"  - Internal annotation files: {outpath}/Ig/internal_data/")
        click.echo(f"  - Auxiliary files: {outpath}/aux_db/")
        click.echo(f"  - Reference index: {outpath}/.references_dataframe.csv.gz")

    except Exception as e:
        click.echo(f"\n✗ Error during database generation: {str(e)}", err=True)
        raise click.ClickException(str(e))

    # Step 6: Verify generated files
    click.echo("\n" + "=" * 60)
    click.echo("STEP 6: Verifying generated files")
    click.echo("=" * 60)

    verification_passed = True

    # Check BLAST databases
    blast_dir = Path(outpath) / "Ig" / "blastdb" if outpath else Path.cwd() / "Ig" / "blastdb"
    if blast_dir.exists():
        species_dirs = list(blast_dir.glob("*/"))
        click.echo(f"✓ Found {len(species_dirs)} species BLAST databases")
        for sp_dir in species_dirs:
            v_files = list(sp_dir.glob("*_V.*"))
            d_files = list(sp_dir.glob("*_D.*"))
            j_files = list(sp_dir.glob("*_J.*"))
            click.echo(f"  - {sp_dir.name}: V={len(v_files)>0}, D={len(d_files)>0}, J={len(j_files)>0}")
    else:
        click.echo("✗ BLAST database directory not found", err=True)
        verification_passed = False

    # Check auxiliary files
    aux_dir = Path(outpath) / "aux_db" / "imgt" if outpath else Path.cwd() / "aux_db" / "imgt"
    if aux_dir.exists():
        aux_files = list(aux_dir.glob("*.aux"))
        click.echo(f"✓ Found {len(aux_files)} auxiliary files")
    else:
        click.echo("✗ Auxiliary file directory not found", err=True)
        verification_passed = False

    # Check internal annotation files
    internal_dir = Path(outpath) / "Ig" / "internal_data" if outpath else Path.cwd() / "Ig" / "internal_data"
    if internal_dir.exists():
        ndm_files = list(internal_dir.glob("*/*.ndm.imgt"))
        click.echo(f"✓ Found {len(ndm_files)} internal annotation files")
    else:
        click.echo("✗ Internal annotation directory not found", err=True)
        verification_passed = False

    # Check reference index
    ref_index = (
        Path(outpath) / ".references_dataframe.csv.gz" if outpath else Path.cwd() / ".references_dataframe.csv.gz"
    )
    if ref_index.exists():
        click.echo(f"✓ Reference index file created: {ref_index}")
    else:
        click.echo("✗ Reference index file not found", err=True)
        verification_passed = False

    if verification_passed:
        click.echo("\n" + "=" * 60)
        click.echo("✅ All database files generated successfully!")
        click.echo("=" * 60)
        click.echo(f"\nDatabase location: {outpath}")
        click.echo("\nTo use these databases with SADIE AIRR:")
        click.echo(f"  export IGDATA={outpath}/Ig")
        click.echo("  sadie airr <your_sequence.fasta>")
    else:
        click.echo("\n" + "=" * 60)
        click.echo("⚠️  Some files were not generated correctly", err=True)
        click.echo("=" * 60)


@sadie.group()
def germlines() -> None:
    """Germline database management commands."""
    pass


@germlines.command("populate")
@click.option(
    "--provider",
    "-p",
    type=click.Choice(["imgt", "ogrdb", "vdjbase", "all"]),
    default="all",
    show_default=True,
    help="Provider to download data from",
)
@click.option(
    "--species",
    "-s",
    multiple=True,
    help="Specific species to download (can be repeated)",
)
@click.option(
    "--force",
    "-f",
    is_flag=True,
    help="Force re-download even if up-to-date",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show what would be downloaded without downloading",
)
def germlines_populate(
    provider: str,
    species: tuple,
    force: bool,
    dry_run: bool,
) -> None:
    """Download germline data from providers.

    Downloads germline sequences from IMGT, OGRDB, and/or VDJbase providers
    and builds the necessary databases for AIRR annotation and renumbering.

    \b
    Examples:
        sadie germlines populate                    # Download all providers
        sadie germlines populate -p imgt            # Download IMGT only
        sadie germlines populate -p imgt -s human   # Download human from IMGT
        sadie germlines populate --dry-run          # Show what would be downloaded
        sadie germlines populate --force            # Force re-download
    """
    from sadie.germlines.cli import populate_germlines

    populate_germlines(provider, list(species) if species else None, force, dry_run)


@germlines.command("status")
def germlines_status() -> None:
    """Show status of germline databases."""
    from rich.console import Console
    from rich.table import Table

    from sadie.germlines.cli import get_local_version, is_up_to_date

    console = Console()

    table = Table(title="Germline Database Status")
    table.add_column("Provider", style="cyan")
    table.add_column("Version")
    table.add_column("Downloaded At")
    table.add_column("Species")
    table.add_column("Status")

    for prov_name in ["imgt", "ogrdb", "vdjbase"]:
        version_info = get_local_version(prov_name)

        if version_info:
            version = version_info.get("version", "unknown")
            downloaded_at = version_info.get("downloaded_at", "unknown")[:10]
            species_count = version_info.get("species_count", 0)
            status = "[green]Up-to-date[/green]" if is_up_to_date(prov_name) else "[yellow]Update available[/yellow]"
        else:
            version = "-"
            downloaded_at = "-"
            species_count = 0
            status = "[red]Not downloaded[/red]"

        table.add_row(prov_name, version, downloaded_at, str(species_count), status)

    console.print(table)


@sadie.group("regenerate-tests")
def regenerate_tests() -> None:
    """Regenerate test fixtures for unit and integration tests.

    These commands regenerate expected output files used by tests.
    Run after making changes to germline databases or AIRR annotation logic.
    """
    pass


@regenerate_tests.command("igl")
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Show detailed output",
)
@click.option(
    "--run-test",
    is_flag=True,
    help="Run associated tests after regeneration",
)
@click.option(
    "--no-backup",
    is_flag=True,
    help="Skip backing up existing files",
)
def regenerate_igl(verbose: bool, run_test: bool, no_backup: bool) -> None:
    """Regenerate IGL reference fixtures for macaque tests.

    Regenerates:
      - tests/data/fixtures/airr_tables/igl_out.feather
      - tests/data/fixtures/airr_tables/bum_link_solution.feather

    These fixtures are used by test_hard_igl_seqs and test_hard_igl_seqs_linked.

    \b
    Examples:
        sadie regenerate-tests igl                    # Regenerate fixtures
        sadie regenerate-tests igl --run-test         # Regenerate and run tests
        sadie regenerate-tests igl --verbose          # Show detailed output
    """
    import shutil
    from datetime import datetime

    import pandas as pd

    from sadie.airr import AirrTable, LinkedAirrTable
    from sadie.airr import methods as airr_methods

    # get_project_root() returns src/sadie, so go up two levels for actual project root
    project_root = get_project_root().parent.parent
    fixtures_dir = project_root / "tests" / "data" / "fixtures" / "airr_tables"

    click.echo("=" * 60)
    click.echo("Regenerating IGL Reference Fixtures")
    click.echo("=" * 60)

    # Define fixtures to regenerate
    fixtures = [
        {
            "name": "IGL Single",
            "input": fixtures_dir / "bum_igl_assignment_macaque.feather",
            "output": fixtures_dir / "igl_out.feather",
            "table_type": "AirrTable",
        },
        {
            "name": "IGL Linked",
            "input": fixtures_dir / "bum_link_input.feather",
            "output": fixtures_dir / "bum_link_solution.feather",
            "table_type": "LinkedAirrTable",
        },
    ]

    results = []

    for fixture in fixtures:
        click.echo(f"\n[{fixture['name']}]")

        # Check input exists
        if not fixture["input"].exists():
            click.echo(f"  ERROR: Input file not found: {fixture['input']}", err=True)
            results.append({"name": fixture["name"], "status": "error", "reason": "input not found"})
            continue

        # Backup if exists and not skipped
        if fixture["output"].exists() and not no_backup:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_path = fixture["output"].with_suffix(f".bak.{timestamp}.feather")
            shutil.copy2(fixture["output"], backup_path)
            if verbose:
                click.echo(f"  Backed up to: {backup_path.name}")

        try:
            # Load input
            if verbose:
                click.echo(f"  Loading: {fixture['input'].name}")

            if fixture["table_type"] == "AirrTable":
                table = AirrTable(pd.read_feather(fixture["input"]))
            else:
                table = LinkedAirrTable(pd.read_feather(fixture["input"]))

            # Process
            if verbose:
                click.echo("  Running termini buffers...")
            processed = airr_methods.run_termini_buffers(table)

            if verbose:
                click.echo("  Running IGL assignment...")
            output_df = airr_methods.run_igl_assignment(processed)

            # Save
            output_df.to_feather(fixture["output"])
            click.echo(f"  SUCCESS: {fixture['output'].name} ({len(output_df)} rows)")
            results.append({"name": fixture["name"], "status": "success", "rows": len(output_df)})

        except Exception as e:
            click.echo(f"  ERROR: {str(e)}", err=True)
            results.append({"name": fixture["name"], "status": "error", "reason": str(e)})

    # Summary
    click.echo("\n" + "=" * 60)
    success_count = sum(1 for r in results if r["status"] == "success")
    click.echo(f"Completed: {success_count}/{len(fixtures)} fixtures regenerated")

    # Run tests if requested
    if run_test and success_count > 0:
        click.echo("\n" + "=" * 60)
        click.echo("Running Tests")
        click.echo("=" * 60)

        tests_dir = project_root / "tests" / "unit" / "airr"
        test_cmd = [
            sys.executable,
            "-m",
            "pytest",
            str(tests_dir / "test_airr.py::test_hard_igl_seqs"),
            str(tests_dir / "test_airr.py::test_hard_igl_seqs_linked"),
            "-v",
        ]

        result = subprocess.run(test_cmd, capture_output=False)

        if result.returncode == 0:
            click.echo("\n✅ All tests passed!")
        else:
            click.echo("\n❌ Some tests failed", err=True)
            sys.exit(1)


@regenerate_tests.command("catnap")
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Show detailed output",
)
@click.option(
    "--run-test",
    is_flag=True,
    help="Run associated tests after regeneration",
)
@click.option(
    "--no-backup",
    is_flag=True,
    help="Skip backing up existing files",
)
def regenerate_catnap(verbose: bool, run_test: bool, no_backup: bool) -> None:
    """Regenerate CATNAP reference fixtures for integration tests.

    Regenerates:
      - tests/data/fixtures/airr_tables/catnap_heavy_airrtable.feather
      - tests/data/fixtures/airr_tables/catnap_light_airrtable.feather

    \b
    Examples:
        sadie regenerate-tests catnap                 # Regenerate fixtures
        sadie regenerate-tests catnap --run-test     # Regenerate and run tests
    """
    import shutil
    from datetime import datetime

    # get_project_root() returns src/sadie, so go up two levels for actual project root
    project_root = get_project_root().parent.parent
    fixtures_dir = project_root / "tests" / "data" / "fixtures"

    click.echo("=" * 60)
    click.echo("Regenerating CATNAP Reference Fixtures")
    click.echo("=" * 60)

    fixtures = [
        {
            "name": "CATNAP Heavy",
            "input": fixtures_dir / "fasta_inputs" / "catnap_nt_heavy.fasta",
            "output": fixtures_dir / "airr_tables" / "catnap_heavy_airrtable.feather",
        },
        {
            "name": "CATNAP Light",
            "input": fixtures_dir / "fasta_inputs" / "catnap_nt_light.fasta",
            "output": fixtures_dir / "airr_tables" / "catnap_light_airrtable.feather",
        },
    ]

    results = []
    airr_api = Airr("human", adaptable=True)

    for fixture in fixtures:
        click.echo(f"\n[{fixture['name']}]")

        if not fixture["input"].exists():
            click.echo(f"  ERROR: Input file not found: {fixture['input']}", err=True)
            results.append({"name": fixture["name"], "status": "error"})
            continue

        # Backup
        if fixture["output"].exists() and not no_backup:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_path = fixture["output"].with_suffix(f".bak.{timestamp}.feather")
            shutil.copy2(fixture["output"], backup_path)
            if verbose:
                click.echo(f"  Backed up to: {backup_path.name}")

        try:
            if verbose:
                click.echo(f"  Processing: {fixture['input'].name}")

            result = airr_api.run_fasta(str(fixture["input"]))
            result.to_feather(fixture["output"])

            click.echo(f"  SUCCESS: {fixture['output'].name} ({len(result)} rows)")
            results.append({"name": fixture["name"], "status": "success", "rows": len(result)})

        except Exception as e:
            click.echo(f"  ERROR: {str(e)}", err=True)
            results.append({"name": fixture["name"], "status": "error"})

    # Summary
    click.echo("\n" + "=" * 60)
    success_count = sum(1 for r in results if r["status"] == "success")
    click.echo(f"Completed: {success_count}/{len(fixtures)} fixtures regenerated")


@regenerate_tests.command("all")
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Show detailed output",
)
@click.option(
    "--run-tests",
    is_flag=True,
    help="Run all tests after regeneration",
)
@click.option(
    "--no-backup",
    is_flag=True,
    help="Skip backing up existing files",
)
def regenerate_all(verbose: bool, run_tests: bool, no_backup: bool) -> None:
    """Regenerate all test fixtures.

    Runs all fixture regeneration commands in sequence:
      1. IGL fixtures (macaque tests)
      2. CATNAP fixtures (integration tests)

    \b
    Examples:
        sadie regenerate-tests all                    # Regenerate all fixtures
        sadie regenerate-tests all --run-tests       # Regenerate and run tests
    """
    from click.testing import CliRunner

    runner = CliRunner()

    click.echo("=" * 60)
    click.echo("Regenerating All Test Fixtures")
    click.echo("=" * 60)

    # Build args
    args = []
    if verbose:
        args.append("--verbose")
    if no_backup:
        args.append("--no-backup")

    # Regenerate IGL
    click.echo("\n>>> IGL Fixtures")
    result = runner.invoke(regenerate_igl, args, catch_exceptions=False)
    click.echo(result.output)

    # Regenerate CATNAP
    click.echo("\n>>> CATNAP Fixtures")
    result = runner.invoke(regenerate_catnap, args, catch_exceptions=False)
    click.echo(result.output)

    # Run tests if requested
    if run_tests:
        click.echo("\n" + "=" * 60)
        click.echo("Running Full Test Suite")
        click.echo("=" * 60)

        project_root = get_project_root().parent.parent
        test_cmd = [
            sys.executable,
            "-m",
            "pytest",
            str(project_root / "tests" / "unit" / "airr"),
            "-v",
            "--tb=short",
        ]

        result = subprocess.run(test_cmd, capture_output=False)

        if result.returncode == 0:
            click.echo("\n✅ All tests passed!")
        else:
            click.echo("\n❌ Some tests failed", err=True)
            sys.exit(1)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    sadie(verbose=3)
