import logging
import os

import click

from .util import get_verbosity_level
from .airr import Airr


@click.command()
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
    required=True,
    type=click.Path(exists=True, file_okay=True, readable=True, resolve_path=True),
    help="""The input file can be compressed or uncompressed file of fasta""",
)
@click.option(
    "--species",
    "-s",
    type=click.Choice(Airr.get_available_species()),
    help="Species to annotate",
    default="human",
)
@click.option(
    "--scfv",
    "-p",
    is_flag=True,
    help="Does the input sequence contain phagemid sequences?",
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
    type=click.Choice(["gzip", "bz2", "none"]),
    help="file compression on output",
    default="none",
)
@click.option(
    "--file_format",
    "-f",
    type=click.Choice(["json", "csv", "gb"]),
    help="output file type format",
    default="csv",
)
def run_airr(verbose, query, species, scfv, out, compress, file_format):
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    logger = logging.getLogger("Airr")

    # No reason to use click echo over print except to show e can
    click.echo("Logging with level=>{}".format(logging.getLevelName(logger.getEffectiveLevel())))
    logger.info("Running Airr on %s using %s", query, species)

    # setup object
    airr = Airr(species)

    # Run annotations on a file
    if scfv:
        scfv_airr = airr.run_file(query, scfv=True)
    else:
        airr_table = airr.run_file(query)

    if compress == "none":
        compress = None

    # Setup output if not specified
    if not out:

        def lookup(x):
            return {"gzip": "gz", "bz2": "bz2"}[x]

        out = os.path.basename(query)
        prefix = out.split(".")[0]
        if not compress:
            out = f"{prefix}.{file_format}"
            out_heavy = f"{prefix}_heavy.{file_format}"
            out_light = f"{prefix}_light.{file_format}"
            out_scfv = f"{prefix}_scfv.{file_format}"
        else:
            out = f"{prefix}.{file_format}.{lookup(compress)}"
            out_heavy = f"{prefix}_heavy.{file_format}.{lookup(compress)}"
            out_light = f"{prefix}_light.{file_format}.{lookup(compress)}"
            out_scfv = f"{prefix}_scfv.{file_format}.{lookup(compress)}"

        if scfv:
            out = out_scfv

    if file_format == "json":
        logger.info("Writing %s", out)
        if scfv:
            scfv_airr.to_json(out, compression=compress)
        else:
            airr_table.to_json(out, compression=compress)

    elif file_format == "csv":
        if scfv:
            logger.info("Writing %s and %s", out_heavy, out_light)
            scfv_airr.to_csv(out, compression=compress)
        else:
            logger.info("Writing %s", out)
            airr_table.to_csv(out, compression=compress)
    elif file_format == "gb":
        if scfv:
            # scfv_airr.to_genbank(out_scfv, compression=compress)
            scfv_airr.to_genbank(out)
        else:
            logger.info("Writing %s", out)
            airr_table.to_genbank(out, compression=compress)


if __name__ == "__main__":
    run_airr()
