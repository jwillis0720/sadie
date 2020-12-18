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
from ..utility.util import get_verbosity_level


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
    help="Output path to generate blast_db, internal_data and auxillary files",
    type=click.Path(resolve_path=True, dir_okay=True, writable=True),
)
def make_igblast_reference(verbose, outpath, only_functional=True):
    """make the igblast reference files

    This script will make the imgt reference files used by igblast or airr, including internal data, the blast
    the blast database, and the auxillary files

    Parameters
    ----------
    verbose : string
        the verbosity level (e.g vvv)
    out_path : path
       the output data directory

    Raises
    ------
    Exception
        Missing IMGT fasta db files
    """
    # Set the root logger in the console script
    # Get back a numeric level associated with number of clicks
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    db_path = os.path.abspath(os.path.dirname(__file__) + "/data/ig_database.json.gz")
    click.echo(f"Default outpath {db_path}")

    if not os.path.exists(db_path):
        click.echo("Can't find imgt database")
        raise Exception(f"Can't find  IMGT database {db_path}")

    if not outpath:
        outpath = os.path.abspath(os.path.dirname(__file__) + "/data/germlines")
        click.echo(f"No outpath specified, using {outpath}")

    generate_internal_annotaion_file_from_db(db_path, outpath, only_functional=only_functional)
    click.echo(f"Generated Internal Data {outpath}/internal_data")
    make_igblast_ref_database(db_path, outpath, only_functional=only_functional)
    click.echo("Successfully made blast data")

    # Send it to specialized function
    make_auxillary_file(db_path, outpath)
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
    help="Where is the germline database",
    default=os.path.join(os.path.dirname(__file__), "data", "ig_database.json.gz"),
)
@click.option(
    "--outpath",
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

    # No reason to use click echo over print except to show e can
    generate_genbank(database, outpath)
    click.echo("Done generating files")


@click.command()
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=4,
    help="Vebosity level, ex. -vvvvv for debug level logging",
)
@click.option(
    "--out-path",
    "-o",
    type=click.Path(file_okay=False, dir_okay=True, resolve_path=True, writable=True),
    required=True,
    help="Specify the output path to dump dataframes",
)
@click.option("--imgt_db", type=str, help="Specify the path for imgt_db")
@click.option("--blast-db", type=str, help="Specify the path for blastdb data")
@click.option("--aux-path", type=str, help="Specify the path for auxillary data ")
def make_germline_segments(verbose, out_path, imgt_db, blast_db, aux_path):
    """Generate germline gene segments file

    This script will generate the germline gene segments that are used in the backend by airr.

    Parameters
    ----------
    verbose : str
        verbosity level, e.g vvv
    out_path : path
        the out directory
    imgt_db : path
        the path to the imgt database file  (ex IMGT_DB.csv.gz). Generated with make-imgt-datbasej
    blast_db : path
        the path to the blast db director
    aux_path : path
        the path to the auxillary databasess

    Raises
    ------
    Exception
        Missing IMGT fasta files
    Exception
        Missing IMGT db file
    """

    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    # No reason to use click echo over print except to show e can
    click.echo("Logging with level =>{}".format(logging.getLevelName(numeric_level)))

    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        click.echo("No output path exists, creating")
        os.makedirs(out_path)
        click.echo("Created output path: {}".format(out_path))

    if not imgt_db:
        # grab the imgt databases
        imgt_db = os.path.join(os.path.dirname(__file__), "data/IMGT_REFERENCE.csv.gz")
    if not os.path.exists(imgt_db):
        click.echo("No IMGT DB exists, use `parse-imgt-db` to generate")
        raise Exception("IMGT DB missingsMissing")

    if not blast_db:
        blast_db = os.path.join(os.path.dirname(__file__), "data/germlines/blastdb/")
        blast_db_path = os.path.abspath(blast_db)
        if not os.path.exists(blast_db_path):
            raise Exception(
                """
                The path {} does not exist and may have not been generated.
                Please use the biopharma-reference module to generate the blast db from IMGT data,
                and use --blast-db-path to specify it's location (ex germline/blastdb)
                """.format(
                    blast_db_path
                )
            )

    if not aux_path:
        aux_path = os.path.join(os.path.dirname(__file__), "data/germlines/aux_data/")
        aux_path = os.path.abspath(aux_path)
        if not os.path.exists(aux_path):
            raise Exception(
                """
                The path {} does not exist and may have not been generated.
                Please use the biopharma-reference module to generate the  aux_datafrom IMGT data,
                and use --aux-path to specify it's location (ex germline/aux_data)
                """.format(
                    aux_path
                )
            )


if __name__ == "__main__":
    # make_imgt_database()
    # make_igblast_reference()
    make_genebank_files_from_db()
