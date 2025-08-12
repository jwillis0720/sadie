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

__version__ = version("sadie-antibody")


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
@click.option("--name", "-n", type=click.Choice(Airr.get_available_species()), help="Species to annotate", default="human")
@click.option("--skip-igl", is_flag=True, help="Skip the igl assignment")
@click.option("--skip-mutation", is_flag=True, help="Skip the somewhat time instansive mutational analysis")
@click.argument(
    "input_path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True),
)
@click.argument("output_path", required=False, type=click.Path(file_okay=True, dir_okay=False, writable=True, resolve_path=True))
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

    This script will make the imgt reference files used by igblast or airr, including internal data, the blast
    the blast database, and the auxillary files. It uses the reference.yml to configure select genes and species.
    If you update the reference.yml file, run this again.
    """
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
    # Set logging level
    numeric_level = getVerbosityLevel(verbose)
    logging.basicConfig(level=numeric_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logger = logging.getLogger(__name__)

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
                click.echo(f"Error: Reference generation failed", err=True)
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

    # Parse species list
    if species == "all":
        species_list = None  # Process all species in reference.yml
    else:
        species_list = [s.strip() for s in species.split(",")]

    # Step 3: Regenerate CATNAP references if requested
    if regenerate_catnap:
        click.echo("\n" + "=" * 60)
        click.echo("STEP 3: Regenerating CATNAP references")
        click.echo("=" * 60)
        from pathlib import Path

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
        germline_path = reference_object.make_airr_database(Path(outpath) if outpath else Path.cwd())

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


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    sadie(verbose=3)
