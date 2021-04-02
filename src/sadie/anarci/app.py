import logging
from pathlib import Path

import click

from ..utility.util import get_verbosity_level
from .anarci import Anarci

"""This is our main command line tag"""


def _validate_anarci_objects(ctx, param, value):
    """Private method for click context to evaluate comma seperated lists and make sure each field is okay"""
    columns = [c.strip() for c in value.split(",")]
    param_name = param.human_readable_name
    if param_name == "allowed_species":
        avail_columns = Anarci.get_allowed_species()
    elif param_name == "allowed_chains":
        avail_columns = Anarci.get_allowed_chains()
    else:
        raise ValueError(f"{param.human_readable_name} not recognized as a valid param")
    for c in columns:
        if c not in avail_columns:
            raise click.BadOptionUsage(param, f"{c} is not available. Only have, {avail_columns}")
    return columns


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
    "--scheme",
    "-s",
    is_flag=False,
    default="imgt",
    type=click.Choice(Anarci.get_available_numbering_schemes()),
    show_default=True,
    help="The numbering scheme to use.",
)
@click.option(
    "--region",
    "-r",
    is_flag=False,
    default="imgt",
    type=click.Choice(Anarci.get_available_region_definitions()),
    show_default=True,
    help="The framework and cdr defition to use",
)
@click.option(
    "--allowed-species",
    "-a",
    is_flag=False,
    default=",".join(Anarci.get_allowed_species()),
    show_default=True,
    callback=_validate_anarci_objects,
    help="A comma seperated list of species to align against",
)
@click.option(
    "--allowed-chains",
    "-c",
    is_flag=False,
    default=",".join(Anarci.get_allowed_chains()),
    show_default=True,
    callback=_validate_anarci_objects,
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
def run_anarci(verbose, query, scheme, region, allowed_species, allowed_chains, out, compress, file_format):
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    logger = logging.getLogger("Anarci")

    # No reason to use click echo over print except to show e can
    click.echo(f"Logging with level=>{logging.getLevelName(logger.getEffectiveLevel())}")
    logger.info(f"Running Anarci on renumbering: {query}")
    logger.info(f"Allowed-species {allowed_species}")
    logger.info(f"Allowed-chains: {allowed_chains}")
    logger.info(f"Numbering: {scheme}")
    logger.info(f"Region Def: {region}")

    # setup object
    anarci_api = Anarci(
        scheme=scheme, region_assign=region, allowed_chain=allowed_chains, allowed_species=allowed_species
    )

    # run file on query
    anarci_results = anarci_api.run_file(query)

    # deal with output
    # if no output file, name after input
    if out:
        out = Path(out)
        segment_out = str(out.stem) + "_anarci_results" + str(out.suffix)
        align_out = str(out.stem) + "_anarci_alignment" + str(out.suffix)
    else:
        input_path = Path(query)
        if compress and file_format.lower() != "feather":
            # feather can't be compressed
            compress = "." + compress
        else:
            compress = ""
        segment_out = Path(str(input_path.stem) + f"_anarci_segment.{file_format.lower()}{compress}")
        align_out = Path(str(input_path.stem) + f"_anarci_alignment.{file_format.lower()}{compress}")

    # deal with file format
    # csv
    if file_format.lower() == "csv":
        anarci_results.to_csv(segment_out)
        anarci_results.get_alignment_table().to_csv(align_out)

    # json
    elif file_format.lower() == "json":
        anarci_results.to_json(segment_out, orient="records")
        anarci_results.get_alignment_table().to_json(align_out, orient="records")

    # feather
    elif file_format.lower() == "feather":
        anarci_results.reset_index(drop=True).to_feather(segment_out)
        anarci_results.get_alignment_table().reset_index().to_feather(align_out)

    # shouldn't get here, but if they specify a file format that is not recognized using an invoke method
    else:
        raise ValueError(f"{file_format} not recognized")
    logger.info(f"Output: {segment_out}")
    logger.info(f"Output: {align_out}")


if __name__ == "__main__":
    run_anarci()
