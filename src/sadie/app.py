"""This is our main entry point"""
import logging
import sys
import click
from .airr import Airr
from .utility.util import get_verbosity_level
from .utility import SadieIO


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
def airr(ctx, verbose, species, db_type, compress, skip_mutation, in_format, out_format, input, output):
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
        airr_table = Airr.run_mutational_analysis(airr_table, "kabat")

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


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    airr()
