"""This is our main entry point"""
import logging
import click
import sys

from .airr import Airr
from .utility.util import get_verbosity_level


class SadieIO:
    def __init__(self, input_object, output_object, out_format, compressed):
        self.input = input_object
        self.output = output_object
        self.compressed = compressed

    def get_input_file_name(self):
        return self.input

    def get_output_file_name(self):
        return self.output


@click.group()
@click.pass_context
def sadie(ctx):
    pass


@sadie.command()
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
@click.option("--gene-type", type=click.Choice(["all", "funtional"]), default="all", show_default=True)
@click.option(
    "--db-type",
    type=click.Choice(["imgt", "custom"]),
    default="imgt",
    show_default=True,
    help="Does the input an scfv sequence contain phagemid sequences?",
)
@click.option(
    "--compress",
    "-z",
    type=click.Choice(["gzip", "bz2", "none"]),
    help="file compression on output",
    default="none",
)
@click.option(
    "--out-format",
    type=click.Choice(["infer", "json", "csv", "feather"]),
    help="output file type format",
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
def airr(ctx, verbose, species, gene_type, db_type, compress, out_format, input, output):
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    airr = Airr(species=species, functional=gene_type, database=db_type)
    io = SadieIO(input, output, compress, out_format)

    # airr file handling
    airr_table = airr.run_file(io.get_input_file_name())
    airr_table = Airr.run_mutational_analysis(airr_table, "kabat")
    if io.output:
        airr_table.to_csv(io.get_output_file_name())
    else:
        sys.stdout.write(airr_table.to_csv())


if __name__ == "__main__":
    airr()
