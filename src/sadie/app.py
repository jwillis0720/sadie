"""This is our main entry point"""
import logging
from pathlib import Path
import sys
import os
from typing import Union
import click

# airr
from sadie.airr import Airr
from sadie.airr.methods import run_igl_assignment, run_mutational_analysis

# utility
from sadie.utility import SadieInputDir, SadieInputFile, SadieOutput
from sadie.utility.util import get_verbosity_level, get_project_root

# reference
from sadie.reference import Reference


@click.group()
@click.pass_context
def sadie(ctx: click.Context) -> None:
    pass


@sadie.command("airr")
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=5,
    help="Vebosity level, ex. -vvvvv for debug level logging",
)
@click.option(
    "--species",
    "-s",
    type=click.Choice(Airr.get_available_species()),
    help="Species to annotate",
    default="human",
)
@click.option("--db-type", type=click.Choice(["imgt", "custom"]), default="imgt", show_default=True)
@click.option("--skip-igl", is_flag=True, help="Skip the igl assignment")
@click.option("--skip-mutation", is_flag=True, help="Skip the somewhat time instansive mutational analysis")
@click.argument(
    "input_path",
    required=True,
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
)
@click.argument(
    "output_path",
    required=False,
    type=click.Path(
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
)
def airr(
    verbose: int,
    species: str,
    db_type: str,
    skip_igl: bool,
    skip_mutation: bool,
    input_path: Path,
    output_path: Union[None, Path, str],
) -> None:

    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    airr = Airr(species=species, database=db_type)
    if not output_path:
        output_path = "stdout"
    else:
        output_path = Path(output_path)

    input_path = Path(input_path)
    if input_path.is_dir():
        input_object = SadieInputDir(input_path)
        if output_path is None:
            output_path = input_path / Path(input_path.name + ".tsv.gz")
        _records = input_object.get_combined_seq_records()

    else:  # we listed a signelf iel
        input_object = SadieInputFile(input_path)  # type: ignore[assignment]
        if output_path is None:
            output_path = input_path.parent / Path(input_path.stem + ".tsv.gz")
        _records = input_object.get_seq_records()  # type: ignore[attr-defined]
    output_object = SadieOutput(output_path)
    airr_table = airr.run_records(_records)

    if not skip_mutation:
        airr_table = run_mutational_analysis(airr_table, "kabat")

    if not skip_igl:
        airr_table = run_igl_assignment(airr_table)

    if output_object.output_format == "stdout":
        output_str = str(airr_table.to_csv(sep="\t"))
        sys.stdout.write(output_str)
    else:
        airr_table.to_output(output_object)
    # handle output


@sadie.group()
@click.pass_context
def reference(ctx: click.Context) -> None:
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
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    click.echo(f"reference path {reference}")

    if not outpath:
        outpath = Path(__file__).parent.joinpath("reference/data/germlines").resolve()
        click.echo(f"No outpath specified, using {outpath}")

    # make reference
    click.echo(f"Getting G3 genes from {reference}")

    # read in yaml file with all statric reference data
    reference_object = Reference.parse_yaml(reference)
    germline_path = reference_object.make_airr_database(outpath)
    click.echo(f"Wrote germline database to {germline_path}")
    click.echo("Done!")


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    sadie()
