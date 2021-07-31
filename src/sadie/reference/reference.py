import logging
from pathlib import Path
from typing import List, Optional
from urllib.parse import quote as url_quote

import pandas as pd
import requests
from pydantic import BaseModel, validator
from sadie.reference.settings import IMGT_LOOKUP
from yaml import load

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

# reference logger
logger = logging.getLogger("reference")


class Error(Exception):
    pass


class G3Error(Exception):
    """Exception for G3"""

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}"


class YamlRef:
    def __init__(self, filepath=None):

        if not filepath:
            self.ref_path = Path(__file__).parent.joinpath("data/reference.yml")
        else:
            self.ref_path = filepath

        self.yaml = load(open(self.ref_path), Loader=Loader)

    def get_reference_types(self) -> set:
        """Return reference type. ex imgt or custom

        Returns
        -------
        set
            unique types
        """
        return set(self.yaml.keys())

    def get_species_keys(self, reference_type: str) -> list:
        """Return functional types from reference types and functional type

        Example
        -------
        yaml.get_species_keys('imgt','functional')
        >>> {'alpaca','cat','dog'...}

        Returns
        -------
        set
            unique types
        """
        return list(self.yaml.get(reference_type).keys())

    def get_sub_species(self, reference_type: str, species: str) -> list:
        """
        get sub species keys. For instance for humanized mouse, a sub species will be mouse and human

        Parameters
        ----------
        reference_type : str
            ex. imgt
        functional : str
            ex. all
        species : str
            cat

        Returns
        -------
        int
            number of sub species
        """
        return list(self.yaml.get(reference_type).get(species).keys())

    def get_genes(self, reference_type: str, species: str, sub_species: str) -> list:
        """Get the genes associated with these keys

        Parameters
        ----------
        reference_type : str
            ex. imgt
        species : str
            ex. human
        sub_species : str
            ex. human

        Returns
        -------
        list
            list of genes

        Examples
        --------
        object.get_genes('imgt'','human','human')
        >>> ['IGHV1-2*01','IGHV1-69*01'....]
        """
        return self.yaml.get(reference_type).get(species).get(sub_species)

    def get_gene_segment(self, reference_type: str, species: str, sub_species: str, gene_segment: str) -> list:
        """Get the genes associated with these keys

        Parameters
        ----------
        reference_type : str
            ex. imgt
        species : str
            ex. human
        sub_species : str
            ex. human
        gene_segment: str
            ex. V
        Returns
        -------
        list
            list of genes of the gene segment

        Examples
        --------
        object.get_gene_segment('imgt','human','human','V')
        >>> ['IGHV1-2*01','IGHV1-69*01'....]
        """
        return list(filter(lambda x: x[3] == gene_segment, self.get_genes(reference_type, species, sub_species)))

    def __repr__(self):
        return self.yaml.__repr__()

    def __iter__(self):
        for database in self.yaml:
            for species in self.yaml[database]:
                for sub_species in self.yaml[database][species]:
                    yield {
                        "database": database,
                        "species": species,
                        "sub_species": sub_species,
                        "gene": self.yaml[database][species][sub_species],
                    }


class Species(BaseModel):
    species: str

    @validator("species")
    def check_species(cls, v):
        # pylint: disable=no-self-argument
        available = list(IMGT_LOOKUP.keys()) + ["clk", "bat64", "hugl18", "se09", "se0916", "se16", "se684", "custom"]
        if v not in available:
            raise ValueError(f"{v} not in species list, have {available}, can use 'custom' for custom dataset")
        return v


class GeneEntry(BaseModel):
    """V,D or J Gene Entry with validation"""

    sub_species: Optional[str]
    species: str
    gene: str
    database: str

    @validator("sub_species")
    def check_sub_species(cls, v):
        # pylint: disable=no-self-argument
        return Species(**{"species": v}).species

    @validator("species")
    def check_species(cls, v, values):
        # pylint: disable=no-self-argument
        if values["sub_species"] is None:
            values["sub_species"] = Species(**{"species": v}).species
        return Species(**{"species": v}).species

    @validator("gene")
    def check_vgene(cls, v, values):
        # pylint: disable=no-self-argument
        if v[3] not in ["V", "D", "J"]:
            raise ValueError(f"gene must contain V,D or J at 3rd index, current have {v[3]} in {v} ")
        return v

    @validator("database")
    def check_database(cls, v):
        # pylint: disable=no-self-argument
        if v not in ["imgt", "custom"]:
            raise ValueError(f"{v} is not a valid database, chocies are 'imgt' or 'custom'")
        return v


class GeneEntries(BaseModel):
    """V,D or J Gene Entry with validation"""

    sub_species: Optional[str]
    species: str
    gene: List[str]
    database: str

    @validator("sub_species")
    def check_sub_species(cls, v):
        # pylint: disable=no-self-argument
        return Species(**{"species": v}).species

    @validator("species")
    def check_species(cls, v, values):
        # pylint: disable=no-self-argument
        if values["sub_species"] is None:
            values["sub_species"] = Species(**{"species": v}).species
        return Species(**{"species": v}).species

    @validator("gene", each_item=True)
    def check_vgene(cls, v, values):
        # pylint: disable=no-self-argument
        if v[3] not in ["V", "D", "J"]:
            raise ValueError(f"gene must contain V,D or J at 3rd index, current have {v[3]} in {v} ")
        return v

    @validator("database")
    def check_database(cls, v):
        # pylint: disable=no-self-argument
        if v not in ["imgt", "custom"]:
            raise ValueError(f"{v} is not a valid database, chocies are 'imgt' or 'custom'")
        return v


class Reference:
    """Factory class to generate reference data objects"""

    def __init__(self):
        self.data = []
        self.endpoint = "https://g3.jordanrwillis.com/api/v1/genes"

    def add_gene(self, gene: dict):
        """Add a single gene to the reference data

        Parameters
        ----------
        gene : dict
            ex. `gene` should contain the following keys: {'species', 'sub_species', 'gene', 'database'}

        Examples
        --------
        Reference.add_gene({"species": "human", "sub_species": "human", "gene": "IGHV1-2*01", "database": "imgt"})
        """
        gene_valid = GeneEntry(**gene)

        # add dictionaries to list from G3
        self.data.append(self._get_gene(gene_valid))

    def add_genes(self, genes: dict):
        genes_valid = GeneEntries(**genes)
        self.data += self._get_genes(genes_valid)

    def _g3_get(self, query):
        response = requests.get(query)
        if response.status_code != 200:
            if response.status_code == 404:
                raise G3Error(f"{response.url} not found in G3")
            raise G3Error(f"{response.url} error G3 database response: {response.status_code}\n{response.text}")
        return response.status_code, response.json()

    def _get_gene(self, gene: GeneEntry) -> dict:
        if not isinstance(gene, GeneEntry):
            raise ValueError(f"{gene} is not GeneEntry")
        gene_url = url_quote(gene.gene)
        query = f"{self.endpoint}?source={gene.database}&common={gene.sub_species}&gene={gene_url}"
        status_code, response_json = self._g3_get(query)
        logger.debug(f"{gene.database}:{gene.species}:{gene.gene} database response: {status_code}")
        if len(response_json) > 1:
            raise G3Error(f"{gene.database}:{gene.species}:{gene.gene} found more than one result")
        response_json[0]["sub_species"] = gene.sub_species
        response_json[0]["species"] = gene.species
        return response_json[0]

    def _get_genes(self, genes: GeneEntries) -> List[dict]:
        """Get multiple genes in one api call from G3 using GeneEntries Model"""
        if not isinstance(genes, GeneEntries):
            raise ValueError(f"{genes} is not GeneEntries")

        # url query
        query = f"{self.endpoint}?source={genes.database}&common={genes.sub_species}&limit=-1"

        # get request as method for future async
        status_code, response_json = self._g3_get(query)
        logger.debug(f"{genes.database}:{genes.species} database response: {status_code}")

        # add the sub_species to the response json
        for document in response_json:
            document["sub_species"] = genes.sub_species
            document["species"] = genes.species
        # only get what the person reqeusted in the json.
        # this is faster than getting individual genes from the g3 api
        filtered_json = list(filter(lambda x: x["gene"] in genes.gene, response_json))
        return filtered_json

    def get_dataframe(self) -> pd.DataFrame:
        """Return a pandas dataframe of the reference data"""
        return pd.json_normalize(self.data)

    def get_database_types(self) -> List[str]:
        """Return a list of the database types in the reference data"""
        return list(set((map(lambda x: x["source"], self.data))))

    @staticmethod
    def parse_yaml(yaml_path: Path = None) -> "Reference":
        """Parse a yaml file into a reference file object

        Parameters
        ----------
        yaml_path : Path to yaml file
        """
        yaml_object = YamlRef(yaml_path)
        ref = Reference()

        # yaml iterator returns a list of genes for every species
        for genes_entry in yaml_object:
            logger.info(f"{genes_entry['database']}-{genes_entry['species']} quering from G3 gateway")
            ref.add_genes(genes_entry)
        return ref

    @staticmethod
    def read_file(path: Path, type="csv") -> "Reference":
        """Read file into a reference object

        Parameters
        ----------
        path : Path to file
        """
        ref = Reference()
        if type == "csv":
            ref.data = pd.read_csv(path, index_col=0).to_dict(orient="records")
        elif type == "json":
            ref.data = pd.json(path).to_dict(orient="records")
        elif type == "feather":
            ref.data = pd.read_feather(path).to_dict(orient="records")
        else:
            raise ValueError(f"{type} is not a valid file type")
        return ref


def get_database(source: str) -> list:
    if source not in ["custom", "imgt"]:
        raise ValueError("Invalid database source, needs to be 'custom' or 'imt'")
    response = requests.get(f"https://g3.jordanrwillis.com/api/v1/genes?source={source}&limit=-1")
    logger.info(f"{source} database response: {response.status_code}")
    return response.json()


def get_loaded_database() -> dict:
    """Get G3 database , wrapped in this function so we don't call on response when module is loaded

    Returns
    -------
    dict
        Dictionary of G3 database. custom or imgt
    """
    return {"custom": get_database("custom"), "imgt": get_database("imgt")}
