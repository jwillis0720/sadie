import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from urllib.parse import quote as url_quote

import pandas as pd
import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .blast import write_blast_db
from .models import GeneEntries, GeneEntry
from .yaml import YamlRef

# reference logger
logger = logging.getLogger("reference")

# G3 API Endpoint
_endpoint = "https://g3.jordanrwillis.com/api/v1/genes"


class G3Error(Exception):
    """Exception for G3 - helps with being specific"""

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


def _write_out_fasta(sequences: List[SeqRecord], outpath: Union[str, Path]) -> Path:
    """Conveinence function to write out a fasta file

    Parameters
    ----------
    sequences : List[SeqRecord]
        List of SeqRecord objects from Biopython
    outpath : Path
        Path to write fasta to

    Returns
    -------
    Path
        The fully qualified path
    """
    logger = logging.getLogger(__file__)
    if isinstance(outpath, str):
        outpath = Path(outpath)
    output_fasta = outpath.with_suffix(".fasta")
    logger.debug("output fasta {}".format(output_fasta))

    with open(output_fasta, "w") as f:
        for sequence in sequences:
            name = sequence.name
            seq = str(sequence.seq).replace(".", "")
            f.write(f">{name}\n{seq}\n")
    return output_fasta


def get_species_from_database(database_json: List[Dict[str, str]]) -> List[str]:
    """Get a list of available species from the G3 database

    Parameters
    ----------
    database_json : List[Dict]
        The G3 database response json. See get_database()
    Returns
    -------
    List[str]
        List of species
    """
    return list(set(map(lambda x: x["common"], database_json)))


def get_database(source: str, custom_endpoint: Union[str, None] = None) -> List[Dict[str, str]]:
    """Get the all entries of either the IMGT or Custom databases. This function calls on the G3 API

    Parameters
    ----------
    source : str
        Either 'imgt' or 'custom'

    Returns
    -------
    list: List of dictionaries where each dictionary is a gene entry

    Raises
    ------
    ValueError
        if source is not 'imgt' or 'custom'
    """
    if not custom_endpoint:
        # if a local endpoint for testint
        custom_endpoint = _endpoint
    request_str = f"{custom_endpoint}?source={source}&limit=-1"
    response = requests.get(request_str)
    if response.status_code != 200:
        # probably put something besides imgt or custom here
        error_log = f"Error loading database: {response.status_code}, {response.text}"
        logger.error(error_log)
        raise G3Error(error_log)
    response_json: List[Dict[str, str]] = response.json()
    return response_json


def get_loaded_database() -> Dict[str, List[Dict[str, str]]]:
    """Get G3 database , wrapped in this function so we don't call on response when module is loaded

    Returns
    -------
    dict
        Dictionary of G3 database. custom or imgt
    """
    return {"custom": get_database("custom"), "imgt": get_database("imgt")}


class Reference:
    """Reference class to handle reference databases for  sadie.airr and sadie.anarci"""

    def __init__(self, endpoint: str = _endpoint):
        """Initialize the reference object

        Parameters
        ----------
        endpoint : str, optional
           The endpoint API address to get the data. Defaults to the G3 API.
        """
        self.data: List[Dict[str, str]] = []
        self.endpoint = endpoint

    @property
    def endpoint(self) -> str:
        return self._endpoint

    @endpoint.setter
    def endpoint(self, endpoint: str) -> None:
        if requests.get(endpoint).status_code != 200:
            raise G3Error(f"{endpoint} is not a valid G3 API endpoint or is down")
        self._endpoint = endpoint

    def add_gene(self, gene: Dict[str, str]) -> None:
        """Add a single gene to the reference data

        Parameters
        ----------
        gene : dict
            ex. `gene` should contain the following keys: {'species', 'sub_species', 'gene', 'database'}

        Examples
        --------
        reference_object = Refrence()
        refrence_object.add_gene({"species": "human", "sub_species": "human", "gene": "IGHV1-2*01", "database": "imgt"})
        """

        # make gene model
        gene_valid = GeneEntry(**gene)

        # add dictionaries to list from G3
        self.data.append(self._get_gene(gene_valid))

    def add_genes(self, genes: Dict[str, Union[str, List[str]]]) -> None:
        """Add a List of genes to the reference data object

        Parameters
        ----------
        genes : List[Dict]
            _description_
        """
        genes_valid = GeneEntries(**genes)  # type: ignore
        self.data += self._get_genes(genes_valid)

    def get_dataframe(self) -> pd.DataFrame:
        """Return a pandas dataframe of the reference data"""
        return pd.json_normalize(self.data)

    def get_database_types(self) -> List[str]:
        """Return a list of the database types in the reference data"""
        return list(set((map(lambda x: x["source"], self.data))))

    @staticmethod
    def parse_yaml(yaml_path: Optional[Path] = None) -> "Reference":
        """Parse a yaml file into a reference file object

        Parameters
        ----------
        yaml_path : Path to yaml file

        Returns
        -------
        Reference - Reference Object
        """
        yaml_object = YamlRef(yaml_path)
        ref = Reference()

        # yaml iterator returns a list of genes for every species
        for genes_entry in yaml_object:
            logger.info(f"{genes_entry['database']}-{genes_entry['species']} quering from G3 gateway")
            ref.add_genes(genes_entry)
        return ref

    @staticmethod
    def read_file(path: Union[Path, str], type: str = "csv") -> "Reference":
        """Read file into a reference object

        Parameters
        ----------
        path : Union[Path,str]
            path to out file

        Examples
        --------
        # read csv
        reference = Reference.read_file("/path/to/file.csv") # can also be file.csv.gz

        # read json
        reference = Reference.read_file("/path/to/file.json") # can also be file.json.gz

        # read feather
        reference = Reference.read_file("/path/to/file.feather") # can also be file.feather

        Returns
        -------
        Reference - Reference Object

        Raises
        ------
        ValueError
            if a csv, json or feather is not passed
        """
        ref = Reference()
        if type == "csv":
            _data = pd.read_csv(path, index_col=0).to_dict(orient="records")
        elif type == "json":
            _data = pd.read_json(path, orient="records").to_dict(orient="records")
        elif type == "feather":
            _data = pd.read_feather(path).to_dict(orient="records")
        else:
            raise ValueError(f"{type} is not a valid file type")
        new_data: List[Dict[str, str]] = list(map(lambda x: {str(i): str(j) for i, j in x.items()}, _data))
        # for entry in _data:
        #     for i, k in entry.items():
        #         new_data.append({str(i): str(k)})
        # ref.data = new_data
        ref.data = new_data
        return ref

    def make_airr_database(self, output_path: Path) -> Path:
        """
        Make the igblastn database, internal database and auxilary database needed by igblast. On success
        return a path to the output database.

        Parameters
        ----------
        output_path : Path
            A path directory to output the database structure

        Returns
        -------
        Path
            On success return path of dumped database file.

        Examples
        --------
        ref_class = Reference()
        ref_class.add_gene({"species": "human", "gene": "IGHV1-69*01", "database": "imgt"})
        ref_class.add_gene({"species": "human", "gene": "IGHD3-3*01", "database": "imgt"})
        ref_class.to_airr_database("/path/to/output/")
        """
        if not self.data:
            # If empty make a reference from the yaml from object and call G3
            logger.warning("Reference data is empty - Generating from yaml")
            self.data = self.parse_yaml().data.copy()
        if isinstance(output_path, str):
            output_path = Path(output_path)
        # dataframe to internal annotation structure
        self._make_internal_annotaion_file(output_path)
        logger.info(f"Generated Internal Data {output_path}/internal_data")
        # dataframe to igblast annotation structure
        self._make_igblast_ref_database(output_path)
        logger.info(f"Generated Blast Data {output_path}/blast")
        # dataframe to igblast aux structure
        self._make_auxillary_file(output_path)
        logger.info(f"Generated Aux Data {output_path}/aux_db")
        return output_path

    def _g3_get(self, query: str) -> Tuple[int, List[Dict[str, str]]]:
        """Use the G3 Restful API

        Parameters
        ----------
        query : str
            query string - ig. https://g3.jordanrwillis.com/api/v1/genes?source=imgt&common=human&gene=IGHV1-69%2A01

        Returns
        -------
        Tuple[int, List[Dict[str,st]]]
            status code and response

        Raises
        ------
        G3Error
            if the response is 404 and we can't find the gene
        G3Error
            Any other response code that is not 200
        """
        response = requests.get(query)
        if response.status_code != 200:
            raise G3Error(f"{response.url} error G3 database response: {response.status_code}\n{response.text}")
        return response.status_code, response.json()

    def _get_gene(self, gene: GeneEntry) -> Dict[str, str]:
        """Get a single gene from the G3 Restful API using a GeneEntry Model

        Parameters
        ----------
        gene : GeneEntry
            The GeneEntry model

        Returns
        -------
         Single Json -> Dict response

        Raises
        ------
        ValueError
            If gene is not a GeneEntry model
        G3Error
            If more than one gene is found, i.e the list is longer than 1. Use _get_genes for more than 1.
        """
        if not isinstance(gene, GeneEntry):
            raise ValueError(f"{gene} is not GeneEntry")

        # change weird characters to url characters
        gene_url = url_quote(gene.gene)

        # we should never have more than one match thanks to the index
        query = f"{self.endpoint}?source={gene.database}&common={gene.sub_species}&gene={gene_url}"

        # use G3 get to return response and json
        status_code, response_json = self._g3_get(query)
        logger.debug(f"{gene.database}:{gene.species}:{gene.gene} database response: {status_code}")

        # put the species in sub species in because they are not a part of G3.
        logger.debug(f"have {len(response_json)} genes")
        response_data: Dict[str, str] = response_json[0]
        response_data["sub_species"] = str(gene.sub_species)
        # sub_species will never return None, only can be assigned NOne
        response_data["species"] = gene.species
        return response_data

    def _get_genes(self, genes: GeneEntries) -> List[Dict[str, str]]:
        """Get a list of genes from entries model. Similar to _get_gene but for multiple genes

        Parameters
        ----------
        genes : GeneEntries
            The GeneEntries model object.

        Returns
        -------
        List[dict]
            A list of Json-> Dict responses from G3


        Raises
        ------
        ValueError
            If genes is not  GeneEntries model
        """
        if not isinstance(genes, GeneEntries):
            raise ValueError(f"{genes} is not GeneEntries")

        # url query
        query = f"{self.endpoint}?source={genes.database}&common={genes.sub_species}&limit=-1"

        # get request as method for future async
        status_code, response_json = self._g3_get(query)
        logger.debug(f"{genes.database}:{genes.species} database response: {status_code}")

        # add the sub_species to the response json
        for document in response_json:
            document["sub_species"] = str(
                genes.sub_species
            )  # must cast this to str since validator converts an optional None to str
            document["species"] = genes.species
        # only get what the person reqeusted in the json.
        # this is faster than getting individual genes from the g3 api
        filtered_json = list(filter(lambda x: x["gene"] in genes.gene, response_json))
        return filtered_json

    def _make_igblast_ref_database(self, outpath: Union[Path, str]) -> None:
        """Generate the IgBlast reference database from the reference object

        Parameters
        ----------
        outpath : Union[Path, str]
            The output path to. example -> path/to/output.
            Then the database will dump to path/to/output/{custom,imgt}/{Ig,TCR}/blastdb
        """
        # The blast DB groups by V,D and J
        logger.debug("Generating from IMGT Internal Database File")

        # get the database as a dataframe
        database = self.get_dataframe()

        # first source, i.e. imgt, custom
        for group, group_df in database.groupby("source"):
            # second group is species, i.e. human, mouse
            for species, species_df in group_df.groupby("species"):
                receptor_blast_dir = Path(outpath) / Path(f"{group}/Ig/blastdb/")
                sub_species_keys = species_df["sub_species"].unique()
                if not receptor_blast_dir.exists():
                    receptor_blast_dir.mkdir()
                for segment, segment_df in species_df.groupby("gene_segment"):
                    chimera = False
                    if len(sub_species_keys) > 1:
                        chimera = True
                    else:
                        chimera = False

                    out_segment = os.path.join(receptor_blast_dir, f"{species}_{segment}")
                    if chimera:
                        # if chimera, join species and sub species with '|'
                        seqs = segment_df.apply(
                            lambda x: SeqRecord(Seq(str(x["sequence"])), name=x["sub_species"] + "|" + x["gene"]),
                            axis=1,
                        ).to_list()
                    else:
                        seqs = segment_df.apply(
                            lambda x: SeqRecord(Seq(str(x["sequence"])), name=x["gene"]), axis=1
                        ).to_list()

                    # write this to a fasta file
                    fasta_file = _write_out_fasta(seqs, out_segment)

                    # Convert fasta file to blast db
                    write_blast_db(fasta_file, Path(str(fasta_file).split(".fasta")[0]))
                    logger.info(f"Wrote blast for {fasta_file}")

    def _make_auxillary_file(self, outpath: Path) -> None:
        """Generate the auxillary file for the IgBlast reference database

        Parameters
        ----------
        outpath : Path
            The output path to. example -> path/to/output.
            Then the database will dump to path/to/output/{custom,imgt}/aux_db/{species}.aux

        Raises
        ------
        ValueError
            if the J region hasn't been added to the database, we refuse to make the aux file
        """

        # get dataframe
        database = self.get_dataframe()
        if database[database.label == "J-REGION"].empty:
            raise ValueError("No J-REGION found in reference object...make sure to add J def")

        # group by source
        for group, group_df in database.groupby("source"):
            receptor_aux_dir = Path(outpath).joinpath(f"{group}/aux_db")
            if not receptor_aux_dir.exists():
                logger.info(f"Creating {receptor_aux_dir}")
                receptor_aux_dir.mkdir(parents=True)
            # next group by species
            for species, species_df in group_df.groupby("species"):
                sub_species_keys = species_df["sub_species"].unique()
                chimera = False
                if len(sub_species_keys) > 1:
                    chimera = True

                aux_file_species = Path(str(receptor_aux_dir) + str(f"/{species}_gl.aux"))
                # get a DF with just common species name
                common_df = species_df[species_df["gene_segment"] == "J"].copy()
                # make sure we don't have any dangling J-REGION
                bad_remainders = common_df[(common_df["imgt.remainder"].isna())]
                if not bad_remainders.empty:
                    logger.warning(
                        f"Had to drop {bad_remainders.shape[0]} rows due to bad remainder for {group}-{species}"
                    )
                    common_df.drop(bad_remainders.index, inplace=True)

                # make columns of an aux databaee
                common_df = common_df[(common_df["imgt.cdr3_end"] != "")]
                common_df.loc[:, "reading_frame"] = common_df["imgt.reading_frame"].astype(int)
                common_df.loc[:, "left_over"] = common_df["imgt.remainder"].astype(int)
                common_df.loc[:, "end"] = common_df["imgt.cdr3_end"].astype(int) - 1
                common_df["marker"] = common_df["gene"].str.split("-").str.get(0).str[0:4].str[::-1].str[:2]
                if chimera:
                    # if chimera add '|' to the species name
                    common_df["gene"] = common_df["common"] + "|" + common_df["gene"]

                # write out the aux file with derived columns
                common_df[["gene", "reading_frame", "marker", "end", "left_over"]].to_csv(
                    aux_file_species, sep="\t", header=None, index=False
                )
                logger.info(f"Wrote aux to {aux_file_species}")

    def _make_blast_db_for_internal(self, df: pd.DataFrame, dboutput: Union[str, Path]) -> None:
        """Make a blast database from dataframe"""
        out_fasta = Path(dboutput).with_suffix(".fasta")
        logger.debug("Writing fasta to {}".format(out_fasta))
        with open(out_fasta, "w") as f:
            for id_, seq in zip(df["gene"], df["sequence"]):
                f.write(">{}\n{}\n".format(id_, seq))
        # get basename
        write_blast_db(out_fasta, dboutput)

    def _make_internal_annotaion_file(self, outpath: Path) -> None:
        """Generate the internal database file for IgBlast

        Parameters
        ----------
        outpath : Path
            The output path to. example -> path/to/output.
            Then the database will dump to path/to/output/{custom,imgt}/{Ig,TCR}/internal_data/{species}/{species}.ndm.imgt

        """
        logger.debug(f"Generating internal annotation file at {outpath}")
        # The internal data file structure goes {db_type}/Ig/internal_path/{species}/

        database = self.get_dataframe()
        for group, group_df in database.groupby("source"):
            # db_type eg. cutom, imgt

            # get a filtered database for V genes
            filtered_data = group_df.loc[group_df["gene_segment"] == "V"]

            # the species is the actual entity we are using for the annotation, e.g se09 or human
            for species, species_df in filtered_data.groupby("species"):
                # species level database
                species_internal_db_path = Path(os.path.join(outpath, group, "Ig", "internal_data", species))
                logger.debug(f"Found species {species}, using {group} database file")
                if not species_internal_db_path.exists():
                    logger.info(f"Creating {species_internal_db_path}")
                    species_internal_db_path.mkdir(parents=True)

                # maybe we requested a chimeric speicies that has more than one sub species, e.g a mouse model
                sub_species_keys = species_df["sub_species"].unique()

                # if we have hybrid species we shall name them with <species>|gene
                gene_df = species_df.copy()
                if len(sub_species_keys) > 1:
                    gene_df["gene"] = gene_df["common"] + "|" + gene_df["gene"]

                # subselect and order
                index_df = gene_df[
                    [
                        "gene",
                        "imgt.fwr1_start",
                        "imgt.fwr1_end",
                        "imgt.cdr1_start",
                        "imgt.cdr1_end",
                        "imgt.fwr2_start",
                        "imgt.fwr2_end",
                        "imgt.cdr2_start",
                        "imgt.cdr2_end",
                        "imgt.fwr3_start",
                        "imgt.fwr3_end",
                    ]
                ].copy()

                # makes everything an integer. sets gene to index so its not affected
                # add +1 to so we get 1-based indexing
                index_df = (index_df.set_index("gene") + 1).astype("Int64").reset_index()

                # drop anything where there is an na in the annotation idnex
                index_df = index_df.drop(index_df[index_df.isna().any(axis=1)].index)
                scheme = "imgt"
                internal_annotations_file_path = os.path.join(species_internal_db_path, f"{species}.ndm.{scheme}")

                # the segment key goes in a weird order
                if len(sub_species_keys) > 1:
                    segment = [i.split("|")[-1].split("-")[0][0:4][::-1][:2] for i in index_df["gene"]]
                else:
                    segment = [i.split("-")[0][0:4][::-1][:2] for i in index_df["gene"]]
                index_df["segment"] = segment
                index_df["weird_buffer"] = 0
                logger.info("Writing to annotation file {}".format(internal_annotations_file_path))
                index_df.to_csv(internal_annotations_file_path, sep="\t", header=False, index=False)
                logger.info("Wrote to annotation file {}".format(internal_annotations_file_path))

                # blast reads these suffixes depending on receptor
                suffix = "V"
                db_outpath = Path(str(species_internal_db_path) + f"/{species}_{suffix}")
                # Pass the dataframe and write out the blast database
                self._make_blast_db_for_internal(gene_df, db_outpath)
