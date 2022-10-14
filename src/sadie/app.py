"""This is our main entry point"""
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, List, Union

import click

# airr
from sadie.airr import Airr
from sadie.airr.methods import run_igl_assignment, run_mutational_analysis

# reference
from sadie.reference.reference import References

# Renumbering
from sadie.renumbering import Renumbering

# utility
from sadie.utility import SadieInputDir, SadieInputFile, SadieOutput
from sadie.utility.util import get_project_root, get_verbosity_level


@click.group()
@click.version_option(
    subprocess.run(["git", "describe", "--tags", "--abbrev=0"], stdout=subprocess.PIPE).stdout.decode("utf-8").strip()
)
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
    "--name",
    "-n",
    type=click.Choice(Airr.get_available_species()),
    help="Species to annotate",
    default="human",
)
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
    name: str,
    skip_igl: bool,
    skip_mutation: bool,
    input_path: Path,
    output_path: Union[None, Path, str],
) -> None:

    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    airr = Airr(name)
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
    numeric_level = get_verbosity_level(verbose)
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
    reference_object = References.from_yaml(reference)
    germline_path = reference_object.make_airr_database(outpath)
    click.echo(f"Wrote germline database to {germline_path}")
    click.echo("Done!")


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    sadie()
