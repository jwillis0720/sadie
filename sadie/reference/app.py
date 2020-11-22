""" Applicaiton for generating annotation references for antibodies"""

# third party
import logging
import os
import click
from sqlalchemy import create_engine

# Module imports
from .imgt import write_imgt_local, get_imgt_file_structure, make_imgt_db
from .internal_data import generate_internal_annotaion_file_from_db
from .igblast_ref import make_igblast_ref_database
from .aux_file import make_auxillary_file
from .genbank import generate_genbank
from .genesegment import generate_v_segment_data, generate_j_gene_data
from ..utility.util import get_verbosity_level


@click.command()
@click.option(
    "--outdir",
    required=True,
    type=click.Path(resolve_path=True, dir_okay=True, writable=True),
    help="Output to generate IMGT local copy",
)
def make_imgt_fasta_files(outdir):
    """Pull all the IMGT fasta files down from IMGT v quest http database

    Arguments:
        outdir {file} -- path to output
    """
    imgt_local_dir = os.path.abspath(outdir)
    click.echo(f"Writing to {imgt_local_dir}")
    write_imgt_local(outdir)
    click.echo(f"Wrote to {imgt_local_dir}")


@click.command()
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=4,
    help="Vebosity level, ex. -vvvvv for debug level logging",
)
@click.option(
    "--imgt-fasta",
    type=click.Path(exists=True, resolve_path=True, dir_okay=True),
    help="Where is the local IMGT from Vquest file path is. Can generate with `get-imgt` script",
)
@click.option(
    "--imgt-data",
    type=click.Path(exists=True, resolve_path=True, dir_okay=True),
    help="Where is the local IMGT data fastas. Should be in package repo",
)
@click.option(
    "--outpath",
    type=click.Path(resolve_path=True, file_okay=True, writable=True),
    help="Where to dump the SQLlite file? Default to reference/data/IMGT_REFERENCE.db",
)
@click.option(
    "--skip-imgt",
    is_flag=True,
    default=False,
    help="skip rebuilding of IMGT reference database and use tables present in sql engine. Will skip right to v and j gene assembly",
)
@click.option(
    "--skip-v",
    default=False,
    is_flag=True,
    help="skip v segment assembly and use whats in engine.",
)
@click.option(
    "--skip-j",
    is_flag=True,
    default=False,
    help="skip j segment assembly and use whats in engine.",
)
def make_imgt_database(verbose, imgt_fasta, imgt_data, outpath, skip_imgt, skip_v, skip_j):
    """[summary]

    Parameters
    ----------
    verbose : [type]
        [description]
    imgt_fasta : [type]
        [description]
    imgt_data : [type]
        [description]
    outpath : [type]
        [description]
    skip_imgt : [type]
        [description]
    skip_v : [type]
        [description]
    skip_j : [type]
        [description]

    Raises
    ------
    Exception
        [description]
    Exception
        [description]
    Exception
        [description]
    """
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)

    if not imgt_fasta:
        imgt_fasta = os.path.join(os.path.dirname(__file__), "data/IMGT_LOCAL/")
        click.echo(f"Attempting to locat IMGT fasta files {imgt_fasta}")

    if not os.path.exists(imgt_fasta):
        click.echo("No IMGT VQuest files exists, use `get-imgt` to download them or specify correct path")
        raise Exception("IMGT Missing")

    if imgt_data:
        if not os.path.exists(imgt_data):
            click.echo("Can't find %s, don't pass argument to use own file", imgt_data)
            raise Exception(f"IMGT Data path not found {imgt_data}")
        imgt_data = os.path.abspath(imgt_data)
    else:
        imgt_data = os.path.join(os.path.dirname(__file__), "data/IMGT_DB/")
        imgt_data = os.path.abspath(imgt_data)
        click.echo(f"Using IMGT data path {imgt_data}")
        if not os.path.exists(imgt_data):
            click.echo("Can't find %s, it's missing from package data", imgt_data)
            raise Exception("IMGT Data path not found")

    if outpath:
        click.echo(f"Writing to {outpath}")
        if not outpath.endswith(".db"):
            outpath = outpath + ".db"
        db_outpath = os.path.abspath(outpath)

    else:
        db_outpath = os.path.dirname(__file__) + "/data/IMGT_REFERENCE.db"
        click.echo(f"Default outpath {db_outpath}")

    engine = create_engine(f"sqlite:////{os.path.abspath(db_outpath)}")

    # Step 1  - Make the IMGT cohesive database from Vquest and IMGT-GB
    if not skip_imgt:
        make_imgt_db(imgt_fasta, imgt_data, engine)
        click.echo(f"Wrote IMGT DB segment data to {engine}")
    else:
        click.echo(f"skipping IMGT DB rebuild and using imgt_unique and imgt_all tables in {engine}")

    # Step 2 - Make V gene segments from data
    if not skip_v:
        generate_v_segment_data(engine)
        click.echo(f"Wrote VSegment segment data to {engine}")
    else:
        click.echo(f"skipping v segmentand using v_segment in {engine}")

    # Step 3 - Make J gene segments from data
    if not skip_j:
        generate_j_gene_data(engine)
        click.echo(f"Wrote JSegment segment data to {engine}")
    else:
        click.echo(f"skipping j segmentand using j_segment in {engine}")
    click.echo("Succeded")


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
    db_outpath = os.path.abspath(os.path.dirname(__file__) + "/data/IMGT_REFERENCE.db")
    click.echo(f"Default outpath {db_outpath}")

    if not os.path.exists(db_outpath):
        click.echo("Can't find imgt database")
        raise Exception(f"Can't find  IMGT database {db_outpath}")
    engine = create_engine(f"sqlite:////{os.path.abspath(db_outpath)}")

    if not outpath:
        outpath = os.path.abspath(os.path.dirname(__file__) + "/data/germlines")
        click.echo(f"No outpath specified, using {outpath}")

    # where to output internal data path
    internal_data_path = os.path.join(outpath, "internal_data")
    if not os.path.exists(internal_data_path):
        click.echo(f"Generating Internal data path {internal_data_path}")

    generate_internal_annotaion_file_from_db(engine, internal_data_path, only_functional=only_functional)
    click.echo(f"Generated Internal Data {outpath}/internal_data")

    blast_dir = os.path.join(outpath, "blastdb")
    if not os.path.exists(blast_dir):
        click.echo(f"Generating blast data path {blast_dir}")
    # # Send it to specialized function
    click.echo(f"Making blast data at {blast_dir}")
    # if only_functional:
    make_igblast_ref_database(engine, blast_dir, only_functional=only_functional)
    click.echo("Successfully made blast data")

    # # Auxiliary data
    aux_path = os.path.join(outpath, "aux_data")
    if not os.path.exists(aux_path):
        click.echo(f"Creating {aux_path}")
        os.makedirs(aux_path)
    # Send it to specialized function
    make_auxillary_file(engine, aux_path)
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
    "--imgt-fasta",
    type=click.Path(exists=True, dir_okay=True, readable=True, resolve_path=True),
    help="Where is the local IMGT file path is. Can generate with `get-imgt` script",
)
@click.option(
    "--imgt-db",
    type=click.Path(exists=True, dir_okay=True, readable=True, resolve_path=True),
    help="Where is the IMGT database file. This comes with the module by default, or you can run with parse-imgt-db",
)
@click.option(
    "--aux-path",
    type=str,
    help="Where is the auxillary data? this can be generated with igblast-setup, get-auxillary-db, or downloaded straight from igblast",
)
@click.option(
    "--outpath",
    type=str,
    help="where to dump the output genbank files?",
    default="./local_gene_reference",
)
def make_genebank_files(verbose, imgt_fasta, imgt_db, aux_path, outpath):
    """Make genebank files from IMGT Fasta files (vquest),IMGT DB files, and auxilary files

    Parameters
    ----------
    verbose : str
        the verbosity level, e.g vvv
    imgt_fasta : path
        path to imgt fasta files
    imgt_db : path
        path to imgt db
    aux_path : path
        path to IGBLAST generated auxilary files
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
    root_logger = logging.getLogger()
    logger = logging.getLogger(__name__)

    if not imgt_fasta:
        imgt_fasta = os.path.join(os.path.dirname(__file__), "data/IMGT_LOCAL/")
        click.echo(f"Attempting to file {imgt_fasta}")

    if not os.path.exists(imgt_fasta):
        # parse v quest fasta file
        click.echo("No IMGT VQuest files exists, use `get-imgt` to download them or specify correct path")
        raise Exception("IMGT Missing")

    if not imgt_db:
        # grab the imgt databases
        imgt_db = os.path.join(os.path.dirname(__file__), "data/IMGT_REFERENCE.csv.gz")
    if not os.path.exists(imgt_db):
        click.echo("No IMGT DB exists, use `parse-imgt-db` to generate")
        raise Exception("IMGT DB missingsMissing")

    if not aux_path:
        aux_path = os.path.join(os.path.dirname(__file__), "data/germlines/aux_data")
        click.echo("No aux path specied using package data")

    if not os.path.exists(aux_path):
        click.echo("No  aux path path `igblast-seutp` to generate or specify the correct path")
        raise Exception("Aux Data Missing")

    if not os.path.exists(outpath):
        logger.info("Creating %s", outpath)

    # No reason to use click echo over print except to show e can
    click.echo("Logging with level =>{}".format(logging.getLevelName(root_logger.getEffectiveLevel())))

    generate_genbank(imgt_fasta, imgt_db, aux_path, outpath)
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
    make_igblast_reference()