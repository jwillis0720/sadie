import logging
import os
from pathlib import Path
from typing import Dict, List, Union
from urllib.parse import quote as url_quote

import pandas as pd
import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from yaml import load

from .blast import write_blast_db
from .models import GeneEntries, GeneEntry

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

# reference logger
logger = logging.getLogger("reference")

# G3 API Endpoint
_endpoint = "https://g3.jordanrwillis.com/api/v1/genes"


def _write_out_fasta(sequences: List[SeqRecord], outpath: Path) -> Path:
    logger = logging.getLogger(__file__)
    output_fasta = outpath + ".fasta"
    logger.debug("output fasta {}".format(output_fasta))

    # Im sure this will come back to haunt me, but if, we've seen the name, save that
    seen_short_names = []
    with open(output_fasta, "w") as f:
        for sequence in sequences:
            name = sequence.name
            seq = str(sequence.seq).replace(".", "")
            f.write(">{}\n{}\n".format(name, seq))
            seen_short_names.append(name)
    return output_fasta


def get_databases_types(database_json) -> List[str]:
    return list(set(map(lambda x: x["source"], database_json)))


def get_species_from_database(database_json) -> List[str]:
    return list(set(map(lambda x: x["common"], database_json)))


def _determine_left_over(X):
    nt = X[0]
    reading_frame = X[1]
    return len(nt[reading_frame:]) % 3


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
    def __init__(self, filepath: Union[None, Path, str] = None):
        if not filepath:
            self.ref_path = Path(__file__).parent.joinpath("data/reference.yml")
        else:
            self.ref_path = filepath

        self.yaml = load(open(self.ref_path), Loader=Loader)

    def get_database_types(self) -> set:
        """Return database types in current yaml

        Example
        -------
        yaml.get_reference_types()
        >>> {'imgt','custom'}

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
        """Iter method will step through the yaml file"""
        for database in self.yaml:
            for species in self.yaml[database]:
                for sub_species in self.yaml[database][species]:
                    yield {
                        "database": database,
                        "species": species,
                        "sub_species": sub_species,
                        "gene": self.yaml[database][species][sub_species],
                    }


class Reference:
    """Reference class to generate reference data objects. These are useful for annotation purposes"""

    def __init__(self):
        self.data = []
        self.endpoint = _endpoint

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

    def make_airr_database(self, output_path: Path) -> Path:
        """
        Make the igblast database, the internal database and the output database needed by igblast. On success
        return a path to the output database.

        Parameters
        ----------
        reference : Reference
            The reference object
        output_path : Path
            A path to dump the igblast reference structure

        Returns
        -------
        Path
            On success return path of dumped database file
        """ """
        Make the igblast database in an output path from the reference object

        Parameters
        ----------
        reference : Reference
            The reference object
        output_path : Path
            A path to dump the igblast reference structure

        Returns
        -------
        Path
            On success return path of dumped database file
        """
        if not self.data:
            raise ValueError("Reference data is empty")
        if isinstance(output_path, str):
            output_path = Path(output_path)
        self._make_internal_annotaion_file(output_path)
        logger.info(f"Generated Internal Data {output_path}/internal_data")
        self._make_igblast_ref_database(output_path)
        logger.info(f"Generated Blast Data {output_path}/blast")
        self._make_auxillary_file(output_path)
        logger.info(f"Generated Aux Data {output_path}/aux_db")
        return output_path

    def _make_igblast_ref_database(self, outpath: Union[Path, str]):
        # The blast DB groups by V,D and J
        logger.debug("Generating from IMGT Internal Database File")
        database = self.get_dataframe()
        for group, group_df in database.groupby("source"):
            for species, species_df in group_df.groupby("species"):
                receptor_blast_dir = os.path.join(outpath, f"{group}/Ig/blastdb/")
                sub_species_keys = species_df["sub_species"].unique()
                if not os.path.exists(receptor_blast_dir):
                    os.makedirs(receptor_blast_dir)
                for segment, segment_df in species_df.groupby("gene_segment"):
                    if len(sub_species_keys) > 1:
                        chimera = True
                    else:
                        chimera = False

                    out_segment = os.path.join(receptor_blast_dir, f"{species}_{segment}")
                    if chimera:
                        seqs = segment_df.apply(
                            lambda x: SeqRecord(Seq(str(x["sequence"])), name=x["sub_species"] + "|" + x["gene"]),
                            axis=1,
                        ).to_list()
                    else:
                        seqs = segment_df.apply(
                            lambda x: SeqRecord(Seq(str(x["sequence"])), name=x["gene"]), axis=1
                        ).to_list()
                    # returns a full fasta path
                    fasta_file = _write_out_fasta(seqs, out_segment)
                    write_blast_db(fasta_file, fasta_file.split(".fasta")[0])
                    logger.info(f"Wrote blast for {fasta_file}")

    def _make_auxillary_file(self, outpath: Path):
        """
        Given the imgt file structure object and aux path, make the aux data
        """
        database = self.get_dataframe()
        if database[database.label == "J-REGION"].empty:
            raise ValueError("No J-REGION found in reference object...make sure to add J def")
        for group, group_df in database.groupby("source"):
            receptor_aux_dir = os.path.join(outpath, f"{group}/aux_db")

            if not os.path.exists(receptor_aux_dir):
                logger.info(f"Creating {receptor_aux_dir}")
                os.makedirs(receptor_aux_dir)
            for species, species_df in group_df.groupby("species"):
                sub_species_keys = species_df["sub_species"].unique()
                if len(sub_species_keys) > 1:
                    chimera = True
                else:
                    chimera = False

                aux_file_species = os.path.join(receptor_aux_dir, f"{species}_gl.aux")
                common_df = species_df[species_df["gene_segment"] == "J"].copy()
                bad_remainders = common_df[(common_df["imgt.remainder"].isna())]
                if not bad_remainders.empty:
                    logger.warning(
                        f"Had to drop {bad_remainders.shape[0]} rows due to bad remainder for {group}-{species}"
                    )
                    common_df.drop(bad_remainders.index, inplace=True)
                common_df = common_df[(common_df["imgt.cdr3_end"] != "")]
                common_df.loc[:, "reading_frame"] = common_df["imgt.reading_frame"].astype(int)
                common_df.loc[:, "left_over"] = common_df["imgt.remainder"].astype(int)
                common_df.loc[:, "end"] = common_df["imgt.cdr3_end"].astype(int) - 1
                common_df["marker"] = common_df["gene"].str.split("-").str.get(0).str[0:4].str[::-1].str[:2]
                if chimera:
                    common_df["gene"] = common_df["common"] + "|" + common_df["gene"]
                common_df[["gene", "reading_frame", "marker", "end", "left_over"]].to_csv(
                    aux_file_species, sep="\t", header=None, index=False
                )
                logger.info(f"Wrote aux to {aux_file_species}")

    def _make_blast_db_for_internal(self, df, dboutput):
        """Make a blast database from dataframe"""
        out_fasta = dboutput + ".fasta"
        logger.debug("Writing fasta to {}".format(out_fasta))
        with open(out_fasta, "w") as f:
            for id_, seq in zip(df["gene"], df["sequence"]):
                f.write(">{}\n{}\n".format(id_, seq))
        out_db = out_fasta.split(".fasta")[0]
        write_blast_db(out_fasta, out_db)

    def _make_internal_annotaion_file(self, outpath: Path):
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
                species_internal_db_path = os.path.join(outpath, group, "Ig", "internal_data", species)
                logger.debug(f"Found species {species}, using {group} database file")
                if not os.path.exists(species_internal_db_path):
                    logger.info(f"Creating {species_internal_db_path}")
                    os.makedirs(species_internal_db_path)

                # maybe we requested a chimeric speicies that has more than one sub species, e.g a mouse model
                sub_species_keys = species_df["sub_species"].unique()

                # if we have hybrid species we shall name them with <species>|gene
                gene_df = species_df.copy()
                if len(sub_species_keys) > 1:
                    gene_df["gene"] = gene_df["common"] + "|" + gene_df["gene"]

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
                # suffix = "TV_V"
                DB_OUTPATH = os.path.join(species_internal_db_path, f"{species}_{suffix}")
                # Pass the dataframe and write out the blast database
                self._make_blast_db_for_internal(gene_df, DB_OUTPATH)


def get_database(source: str) -> List[Dict]:
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
