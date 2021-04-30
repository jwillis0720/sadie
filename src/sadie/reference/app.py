# pylint: disable=no-value-for-parameter
""" Applicaiton for generating annotation references for antibodies"""

# third party
import logging
import os
import click

# Module imports
from .internal_data import generate_internal_annotaion_file_from_db
from .igblast_ref import make_igblast_ref_database
from .aux_file import make_auxillary_file
from .genbank import generate_genbank
from ..utility.util import get_verbosity_level, get_project_root


@click.command()
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
)
@click.option(
    "--reference",
    "-d",
    help="Current database folder path. Will contain two json files",
    type=click.Path(resolve_path=True, dir_okay=True, exists=True),
    default=os.path.join(os.path.dirname(__file__), "data/reference.yml"),
)
def make_igblast_reference(verbose, outpath, reference):
    """make the igblast reference files

    This script will make the imgt reference files used by igblast or airr, including internal data, the blast
    the blast database, and the auxillary files. It uses the reference.yml to configure select genes and species.

    Parameters
    ----------
    verbose : string
        the verbosity level (e.g vvv)
    outpath : path
       the output data directory

    reference : path
        the reference.yml file with gene configurestion
    """
    # Set the root logger in the console script
    # Get back a numeric level associated with number of clicks
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    click.echo(f"reference path {reference}")

    if not outpath:
        outpath = os.path.abspath(os.path.dirname(__file__) + "/data/germlines")
        click.echo(f"No outpath specified, using {outpath}")

    generate_internal_annotaion_file_from_db(reference, outpath)
    click.echo(f"Generated Internal Data {outpath}/internal_data")
    make_igblast_ref_database(reference, outpath)
    click.echo("Successfully made blast data")

    # Send it to specialized function
    make_auxillary_file(reference, outpath)
    click.echo("Successfully made auxillary data")
    click.echo("Done!")


@click.command()
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=4,
    help="Vebosity level, ex. -vvvvv for debug level logging",
)
@click.option(
    "--database",
    type=click.Path(exists=True, dir_okay=True, readable=True, resolve_path=True),
    help="Where is the germline database, Will ",
    default=os.path.join(os.path.dirname(__file__), "../data", "new_database.json.gz"),
)
@click.option(
    "--outpath",
    "-o",
    type=str,
    help="where to dump the output genbank files?",
    default="./local_gene_reference",
)
def make_genebank_files_from_db(verbose, database, outpath):
    """Make genebank files from IMGT Fasta files (vquest),IMGT DB files, and auxilary files

    Parameters
    ----------
    verbose : str
        the verbosity level, e.g vvv
    database: path
        the database path
    outpath : path
        Output path

    Raises
    ------
    Exception
        Missing IMGT files
    Exception
        Missing AUX Path
    """
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)

    if not os.path.exists(database):
        click.echo("No  dtabase path present")
        raise Exception(f"No DB {database}")

    if not os.path.exists(outpath):
        click.echo(f"Creating {outpath}")
        os.makedirs(outpath)

    generate_genbank(database, outpath)
    click.echo("Done generating files")


if __name__ == "__main__":
    make_igblast_reference()
