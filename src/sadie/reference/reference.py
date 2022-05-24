"""This module houses the main Refernce object to manipulate the backend references"""
from __future__ import annotations


import logging
from pathlib import Path
from time import sleep
from typing import Any, Dict, List, Optional, Tuple, Union
from urllib.parse import quote as url_quote

import pandas as pd
import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sadie.reference.models import GeneEntries, GeneEntry
from sadie.reference.yaml import YamlRef
from sadie.reference.util import make_blast_db_for_internal, write_out_fasta, write_blast_db

# reference logger
logger = logging.getLogger("Reference")

# column typing from pandas stubs
Column = Union[int, str]


class G3Error(Exception):
    """Exception for G3 - helps with being specific"""

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class Reference:
    """Reference class to handle reference databases for  sadie.airr and sadie.numbering"""

    # G3 API Endpoint
    _endpoint = "https://g3.jordanrwillis.com/api/v1/genes"

    def __init__(self, endpoint: str = _endpoint):
        """Initialize the reference object

        Parameters
        ----------
        endpoint : str, optional
           The endpoint API address to get the data. Defaults to the G3 API.
        """
        self.data: List[Dict[Column, str] | Dict[str, str]] = []
        self.endpoint = endpoint

    @property
    def endpoint(self) -> str:
        return self._endpoint

    @endpoint.setter
    def endpoint(self, endpoint: str) -> None:
        _counter = 0
        while True:
            _get = requests.get(endpoint)
            if _get.status_code == 503:
                _counter += 1
                sleep(5)
                logger.info(f"Waiting for G3 API {endpoint} to be available --try: {_counter}")
            elif _get.status_code == 200:
                break
            else:
                raise G3Error(f"Error loading G3 API {endpoint}")
            if _counter > 5:
                raise G3Error(f"{endpoint} is not a valid G3 API endpoint or is down")

        logger.info(f"G3 API {endpoint} is available")
        self._endpoint = endpoint

    def add_gene(self, gene: Dict[str, str]) -> None:
        """Add a single gene to the reference data

        Parameters
        ----------
        gene : dict
            ex. `gene` should contain the following keys: {'species',  'gene', 'database'}

        Examples
        --------
        reference_object = Refrence()
        refrence_object.add_gene({"species": "human",  "gene": "IGHV1-2*01", "database": "imgt"})
        """

        # make gene model
        gene_valid = GeneEntry(**gene)

        # add dictionaries to list from G3
        self.data.append(self._get_gene(gene_valid))

    def add_genes(self, species: str, source: str, genes: List[str]) -> None:
        """Add a List of genes to the reference data object from a single species and database

        Parameters
        ----------
        species: str
        genes : List[str]
        database: str

        Examples
        --------
        ref_class = Reference()
        genes = []
        genes.append("IGHV1-69*01")
        genes.append("IGHD3-3*01"
        genes.append("IGHJ6*01")
        ref_class.add_genes('human','imgt',genes)
        """
        genes_valid = GeneEntries(species=species, source=source, genes=genes)
        self.data += self._get_genes(genes_valid)

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
        query = f"{self.endpoint}?source={gene.source}&common={gene.species}&gene={gene_url}"

        # use G3 get to return response and json
        status_code, response_json = self._g3_get(query)
        logger.debug(f"{gene.source}:{gene.species}:{gene.gene} database response: {status_code}")

        # put the species in sub species in because they are not a part of G3.
        logger.debug(f"have {len(response_json)} genes")
        response_data: Dict[str, str] = response_json[0]
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
        query = f"{self.endpoint}?source={genes.source}&common={genes.species}&limit=-1"

        # get request as method for future async
        status_code, response_json = self._g3_get(query)
        logger.debug(f"{genes.source}:{genes.species} database response: {status_code}")

        # this is faster than getting individual genes from the g3 api
        # @Todo, add a find_genes method to G3 rather than pulling all the data and filtering...
        filtered_json = list(filter(lambda x: x["gene"] in genes.genes, response_json))
        return filtered_json

    def get_dataframe(self) -> pd.DataFrame:
        """Return a pandas dataframe of the references data"""
        return pd.json_normalize(self.data)

    @staticmethod
    def from_dataframe(input_df: pd.DataFrame) -> "Reference":
        matched_indexes = pd.Index(
            [
                "_id",
                "source",
                "common",
                "gene",
                "label",
                "gene_segment",
                "receptor",
                "sequence",
                "latin",
                "imgt.sequence_gapped",
                "imgt.sequence_gapped_aa",
                "imgt.cdr3",
                "imgt.cdr3_aa",
                "imgt.fwr4",
                "imgt.fwr4_aa",
                "imgt.cdr3_start",
                "imgt.cdr3_end",
                "imgt.fwr4_start",
                "imgt.fwr4_end",
                "imgt.reading_frame",
                "imgt.ignored",
                "imgt.not_implemented",
                "imgt.expression",
                "imgt.expression_match",
                "imgt.remainder",
                "imgt.imgt_numbering",
                "imgt.sequence",
                "imgt.fwr1",
                "imgt.fwr1_aa",
                "imgt.fwr1_start",
                "imgt.fwr1_end",
                "imgt.cdr1",
                "imgt.cdr1_aa",
                "imgt.cdr1_start",
                "imgt.cdr1_end",
                "imgt.fwr2",
                "imgt.fwr2_aa",
                "imgt.fwr2_start",
                "imgt.fwr2_end",
                "imgt.cdr2",
                "imgt.cdr2_aa",
                "imgt.cdr2_start",
                "imgt.cdr2_end",
                "imgt.fwr3",
                "imgt.fwr3_aa",
                "imgt.fwr3_start",
                "imgt.fwr3_end",
                "imgt.imgt_functional",
                "imgt.contrived_functional",
                "chimera",
            ]
        )
        _diffs = matched_indexes.symmetric_difference(input_df.columns)
        if not _diffs.empty:
            raise ValueError(f"{_diffs} not in the dataframe")

        # fresh instance
        reference = Reference()

        # get dict as lis tof records
        input_list: List[Dict[Column, Any]] = input_df.to_dict(orient="records")

        # can't assign dirrectly so have to append to beat mypy
        for key in input_list:
            reference.data.append(key)
        return reference


class References:
    def __init__(self) -> None:
        self.references: Dict[str, Reference] = {}

    def add_reference(self, name: str, reference: Reference, overwrite: bool = False) -> None:
        if name in self.references.keys():
            if not overwrite:
                raise NameError(f"{name} exists in References. Use overwrite if this is your intention")
            else:
                logger.warning(f"overwriting {name}")
        logger.info(f"Adding {name} to references")
        self.references[name] = reference

    def get_dataframe(self) -> pd.DataFrame:
        """Return a pandas dataframe of the references data"""
        names_dataframe: List[pd.DataFrame] = []
        for name in self.references:
            # create a single reference object
            _ref: Reference = self.references[name]
            _df: pd.DataFrame = _ref.get_dataframe()

            # because a user could add genes multiple times, lets drop by unique id
            logger.info(f"dropping duplicated: {_df['_id'].duplicated().sum()}")
            _df = _df.drop_duplicates("_id")

            # insert name at beggining
            _df.insert(0, "name", name)  # type: ignore
            names_dataframe.append(_df)
        # concat all the dataframes
        concat_df = pd.concat(names_dataframe).reset_index(drop=True)

        # groupby names
        concat_df_groupby_name = concat_df.groupby("name")

        # within names, check how many species are there
        chimeric_gb = concat_df_groupby_name.apply(lambda x: len(x["common"].unique()) > 1)

        list_of_chimera: List[str] = chimeric_gb[chimeric_gb].index.to_list()
        logger.info(f"{list_of_chimera} are chimeric")

        # get the indexes which contain names that are to be chimerized
        indexes_to_chimera = concat_df[concat_df["name"].isin(list_of_chimera)].index

        # set all cherics to false
        concat_df["chimera"] = False
        concat_df.loc[indexes_to_chimera, "chimera"] = True

        # change the gene to common|gene
        concat_df.loc[indexes_to_chimera, "gene"] = concat_df.loc[indexes_to_chimera, ["common", "gene"]].apply(
            lambda x: "|".join(x), axis=1
        )
        return concat_df

    @staticmethod
    def from_yaml(yaml_path: Optional[Path] = None) -> "References":
        """Parse a yaml file into a references file object

        Parameters
        ----------
        yaml_path : Path to yaml file

        Returns
        -------
        Reference - Reference Object
        """
        yaml_ref_object = YamlRef(yaml_path)

        # the yaml object
        yaml_ref = yaml_ref_object.yaml

        # make emtpy references object
        references_object = References()

        # iterate through names
        for name in yaml_ref:
            reference_object = Reference()

            # iterate where they came from
            for source in yaml_ref.get(name):
                # iterate through species within source
                for species in yaml_ref.get(name).get(source):
                    logger.info(f"Adding {species} from {source} to {name}")
                    list_of_genes: List[str] = yaml_ref[name][source][species]
                    # add by list of genes per species given source
                    reference_object.add_genes(species, source, list_of_genes)
            references_object.add_reference(name, reference_object)
        return references_object

    def _make_igblast_ref_database(self, outpath: Union[Path, str]) -> None:
        """Generate the IgBlast reference database from the reference object

        Parameters
        ----------
        outpath : Union[Path, str]
            The output path to. example -> path/to/output.
            Then the database will dump to path/to/output/{Ig,TCR}/blastdb/{name}
        """
        # The blast DB groups by V,D and J
        logger.debug("Generating from IMGT Internal Database File")

        # get the database as a dataframe
        database = self.get_dataframe()
        if database[database.label == "D-REGION"].empty:
            raise ValueError("No D-REGION found in reference object...make sure to add D gene")

        # first name, i.e. "human" or "se09"
        groupby_dataframe = database.groupby("name")
        for name, group_df in groupby_dataframe:
            receptor_blast_dir = Path(outpath) / Path(f"Ig/blastdb/{name}/")
            if not receptor_blast_dir.exists():
                receptor_blast_dir.mkdir(parents=True)
            for segment, segment_df in group_df.groupby("gene_segment"):
                out_segment = receptor_blast_dir.joinpath(f"{name}_{segment}")
                seqs: List[SeqRecord] = segment_df.apply(
                    lambda x: SeqRecord(Seq(str(x["sequence"])), name=x["gene"]), axis=1
                ).to_list()

                # write this to a fasta file
                fasta_file = write_out_fasta(seqs, out_segment)

                # Convert fasta file to blast db
                write_blast_db(fasta_file, Path(str(fasta_file).split(".fasta")[0]))
                logger.info(f"Wrote blast for {fasta_file}")

    def _make_auxillary_file(self, outpath: Path) -> None:
        """Generate the auxillary file for the IgBlast reference database

        Parameters
        ----------
        outpath : Path
            The output path to. example -> path/to/output.
            Then the database will dump to path/to/output/aux_db/{scheme}/{name}.aux

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
        # for now we only have one scheme
        scheme = "imgt"
        receptor_aux_dir = Path(outpath).joinpath(f"aux_db/{scheme}")
        if not receptor_aux_dir.exists():
            logger.info(f"Creating {receptor_aux_dir}")
            receptor_aux_dir.mkdir(parents=True)
        for group, group_df in database.groupby("name"):
            aux_file_name = Path(str(receptor_aux_dir) + str(f"/{group}_gl.aux"))
            # get a DF with just common species name
            common_df = group_df[group_df["gene_segment"] == "J"].copy()
            # make sure we don't have any dangling J-REGION
            bad_remainders = common_df[(common_df["imgt.remainder"].isna())]
            if not bad_remainders.empty:
                logger.warning(f"Had to drop {bad_remainders.shape[0]} rows due to bad remainder for {group}")
                common_df.drop(bad_remainders.index, inplace=True)

            # make columns of an aux databaee common_df = common_df[(common_df["imgt.cdr3_end"] != "")]
            common_df.loc[:, "reading_frame"] = common_df["imgt.reading_frame"].astype(int)
            common_df.loc[:, "left_over"] = common_df["imgt.remainder"].astype(int)
            common_df.loc[:, "end"] = common_df["imgt.cdr3_end"].astype(int) - 1

            # JH, JK, JL
            common_df["marker"] = (
                common_df["gene"].str.split("|").str.get(-1).str.split("-").str.get(0).str[0:4].str[::-1].str[:2]
            )

            # write out the aux file with derived columns
            common_df[["gene", "reading_frame", "marker", "end", "left_over"]].to_csv(
                aux_file_name, sep="\t", header=None, index=False
            )
            logger.info(f"Wrote aux to {aux_file_name}")

    def _make_internal_annotaion_file(self, outpath: Path) -> None:
        """Generate the internal database file for IgBlast

        Parameters
        ----------
        outpath : Path
            The output path to. example -> path/to/output.
            Then the database will dump to path/to/output/{Ig,TCR}/internal_data/{name}/{name}.ndm.imgt
        """
        logger.debug(f"Generating internal annotation file at {outpath}")
        # The internal data file structure goes Ig/internal_path/{name}/

        database = self.get_dataframe()
        for name, group_df in database.groupby("name"):
            # get a filtered database for V genes
            filtered_data = group_df.loc[group_df["gene_segment"] == "V"].copy()

            # the species is the actual entity we are using for the annotation, e.g se09 or human
            name_internal_df_path = Path(outpath).joinpath(Path(f"Ig/internal_data/{name}/"))
            if not name_internal_df_path.exists():
                logger.info(f"Creating {name_internal_df_path}")
                name_internal_df_path.mkdir(parents=True)

            # subselect and order
            index_df = filtered_data[
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
            internal_annotations_file_path = name_internal_df_path.joinpath(f"{name}.ndm.{scheme}")

            segment = [i.split("|")[-1].split("-")[0][0:4][::-1][:2] for i in index_df["gene"]]
            index_df["segment"] = segment
            index_df["weird_buffer"] = 0
            logger.info(f"Writing to annotation file {internal_annotations_file_path}")
            index_df.to_csv(internal_annotations_file_path, sep="\t", header=False, index=False)
            # blast reads these suffixes depending on receptor
            suffix = "V"
            db_outpath = Path(str(name_internal_df_path) + f"/{name}_{suffix}")
            # Pass the dataframe and write out the blast database
            make_blast_db_for_internal(group_df, db_outpath)

    @staticmethod
    def from_json(path: Path | str) -> "References":
        """Read file into a reference object

        Parameters
        ----------
        path : Union[Path,str]
            path to out file

        Examples
        --------
        # read json
        reference = Reference.read_file("/path/to/file.json") # can also be file.json.gz

        Returns
        -------
        Reference - Reference Object
        """
        _data = pd.read_json(path, orient="records").astype(
            {"imgt.ignored": object, "imgt.not_implemented": object, "imgt.expression_match": object}
        )
        return References.from_dataframe(_data)

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
        if not self.references:
            # If empty make a reference from the yaml from object and call G3
            logger.warning("Reference data is empty - Generating from yaml")
            self.references = self.from_yaml().references.copy()
        if isinstance(output_path, str):
            output_path = Path(output_path)
        # dataframe to internal annotation structure
        self._make_internal_annotaion_file(output_path)
        logger.info(f"Generated Internal Data {output_path}/Ig/internal_data")
        # dataframe to igblast annotation structure
        self._make_igblast_ref_database(output_path)
        logger.info(f"Generated Blast Data {output_path}/Ig/blastdb")
        # dataframe to igblast aux structure
        self._make_auxillary_file(output_path)
        logger.info(f"Generated Aux Data {output_path}/aux_db")
        return output_path

    @staticmethod
    def from_dataframe(dataframe: pd.DataFrame) -> "References":
        """Read dataframe into a reference object

        Parameters
        ----------
        dataframe : pd.DataFrame
            dataframe of the Reference file

        Examples
        --------
        reference_df = pd.read_csv("/path/to/file.csv") # can also be file.csv.gz
        reference_object = Reference.from_dataframe(reference_df)

        Returns
        -------
        Reference - Reference Object

        Raises
        ------
        ValueError
            if pd.Dataframe is not suppplied
        """
        references = References()
        for name, name_df in dataframe.groupby("name"):
            name_df["gene"] = name_df["gene"].str.split("|").str[-1]
            ref = Reference().from_dataframe(
                name_df.drop(columns=["name"]).astype(
                    {"imgt.ignored": object, "imgt.not_implemented": object, "imgt.expression_match": object}
                )
            )
            references.add_reference(name, ref)
        return references

    def __repr__(self) -> str:
        return self.get_dataframe().__repr__()
