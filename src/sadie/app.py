"""This is our main entry point"""
import logging
import sys
import os
import click

# airr
from sadie.airr import Airr
from sadie.airr.methods import run_igl_assignment, run_mutational_analysis

# utility
from sadie.utility import SadieIO
from sadie.utility.util import get_verbosity_level, get_project_root

# reference
from sadie.reference import make_germline_database, Reference


@click.group()
@click.pass_context
def sadie(ctx):
    pass


@sadie.command("airr")
@click.pass_context
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
@click.option(
    "--compress",
    "-z",
    type=click.Choice(["gz", "bz2", "infer"]),
    help="file compression for output, gzip, bzip2 or infer from output name",
    default="infer",
)
@click.option("--skip-mutation", is_flag=True, help="Skip the somewhat time instansive mutational analysis")
@click.option("--skip-igl", is_flag=True, help="Skip the igl assignment")
@click.option(
    "--in-format",
    type=click.Choice(["infer", "fasta", "fastq", "ab1"]),
    help="The input format you have passed. ",
    default="infer",
    show_default=True,
)
@click.option(
    "--out-format",
    type=click.Choice(["infer", "json", "csv", "tsv", "feather", "stdout"]),
    help="output file type format. Default is to infer from output argument",
    default="infer",
    show_default=True,
)
@click.argument(
    "input",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True),
)
@click.argument(
    "output", required=False, type=click.Path(file_okay=True, dir_okay=False, writable=True, resolve_path=True)
)
def airr(ctx, verbose, species, db_type, compress, skip_igl, skip_mutation, in_format, out_format, input, output):
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    airr = Airr(species=species, database=db_type)

    # force this so we don't get an error
    if not output and out_format == "infer":
        out_format = "csv"
    io = SadieIO(input, output, in_format, out_format, compress)

    # airr file handling
    # only having a fasta file uncompressed allow calling directly on run_fasta
    if io.input_file_type == "fasta" and not io.input_compressed:
        airr_table = airr.run_fasta(io.input)
    else:
        airr_table = airr.run_records(io.get_input_records())
    # now get mutation analysis
    if not skip_mutation:
        airr_table = run_mutational_analysis(airr_table, "kabat")

    if not skip_igl:
        airr_table = run_igl_assignment(airr_table)

    # handle output
    click.echo(f"Writing {len(airr_table)} to output to {io.output} file")
    if io.infered_output_format == "csv":
        airr_table.to_csv(io.output)
    elif io.infered_output_format == "tsv":
        airr_table.to_csv(io.output, sep="\t")
    elif io.infered_output_format == "feather":
        airr_table.to_feather(io.output)
    elif io.infered_output_format == "json":
        airr_table.to_json(io.output)
    else:
        sys.stdout.write(airr_table.to_csv(sep="\t"))


@sadie.group()
@click.pass_context
def reference(ctx):
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
def make_igblast_reference(verbose, outpath, reference):
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
        outpath = os.path.abspath(os.path.dirname(__file__) + "reference/data/germlines")
        click.echo(f"No outpath specified, using {outpath}")

    # make reference
    # reference_object = Reference.read_file("tmp_dataframe.csv")
    click.echo(f"Getting G3 genes from {reference}")
    reference_object = Reference.parse_yaml(reference)
    germline_path = make_germline_database(reference_object, outpath)
    click.echo(f"Wrote germline database to {germline_path}")
    click.echo("Done!")


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    # airr()
    make_igblast_reference()
