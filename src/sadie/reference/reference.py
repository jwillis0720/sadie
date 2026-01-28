"""This module houses the main Refernce object to manipulate the backend references"""

from __future__ import annotations

import logging
from pathlib import Path
from time import sleep
from typing import Any, Dict, List, Optional, Tuple, Union
from urllib.parse import quote as url_quote

import pandas as pd
import pyhmmer
import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from requests.exceptions import RequestException

from sadie.reference.models import GeneEntries, GeneEntry
from sadie.reference.util import (
    make_blast_db_for_internal,
    write_blast_db,
    write_out_fasta,
)
from sadie.reference.yaml import YamlRef

# reference logger
logger = logging.getLogger("Reference")

# column typing from pandas stubs
Column = Union[Union[int, str], str]


class G3Error(Exception):
    """Exception for G3 - helps with being specific"""

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)


class Reference:
    """Reference class to handle reference databases for  sadie.airr and sadie.numbering"""

    # G3 API Endpoint
    _endpoint = "https://g3.jordanrwillis.com/api/v1/genes"

    def __init__(self, endpoint: str = _endpoint, use_germlines: bool = False):
        """Initialize the reference object

        Parameters
        ----------
        endpoint : str, optional
           The endpoint API address to get the data. Defaults to the G3 API.
           Ignored if use_germlines is True.
        use_germlines : bool, optional
           If True, use local germlines module instead of G3 API. Defaults to False.
        """
        self.data: List[Dict[Column, str] | Dict[str, str]] = []
        self.use_germlines = use_germlines

        if not use_germlines:
            self.endpoint = endpoint
        else:
            # Skip G3 API validation when using germlines
            # Components imported locally in _get_gene/_get_genes
            self._endpoint = endpoint

    @property
    def endpoint(self) -> str:
        return self._endpoint

    @endpoint.setter
    def endpoint(self, endpoint: str) -> None:
        _counter = 0
        while True:
            try:
                _get = requests.get(endpoint)
            except RequestException as e:
                # Handle connection errors (DNS failures, connection refused, etc.)
                raise G3Error(f"Failed to connect to {endpoint}: {str(e)}")

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

        if self.use_germlines:
            from sadie.germlines import GermlineManager

            manager = GermlineManager(providers=[gene_valid.source])
            manager.validate_species(gene_valid.source, gene_valid.species)

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

        if self.use_germlines:
            from sadie.germlines import GermlineManager

            manager = GermlineManager(providers=[genes_valid.source])
            manager.validate_species(genes_valid.source, genes_valid.species)

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
        """Get a single gene from the G3 Restful API or germlines module using a GeneEntry Model

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

        # Use germlines module if enabled
        if self.use_germlines:
            from sadie.germlines import GermlineManager
            from sadie.germlines.g3_adapter import GermlineToG3Adapter

            # Create manager with explicit source (no priority fallback)
            manager = GermlineManager(providers=[gene.source])
            germline_gene = manager.get_gene_by_name(gene.gene, gene.species)

            if not germline_gene:
                raise G3Error(f"Gene {gene.gene} not found in {gene.source} database for species {gene.species}")

            # Transform to G3 format
            adapter = GermlineToG3Adapter()
            g3_dict = adapter.to_g3_format(germline_gene)
            logger.debug(f"Retrieved {gene.gene} from germlines module ({gene.source})")
            return g3_dict

        # Use G3 API (legacy path)
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
            A list of Json-> Dict responses from G3 or germlines module


        Raises
        ------
        ValueError
            If genes is not  GeneEntries model
        """
        if not isinstance(genes, GeneEntries):
            raise ValueError(f"{genes} is not GeneEntries")

        # Use germlines module if enabled
        if self.use_germlines:
            from sadie.germlines import GermlineManager
            from sadie.germlines.g3_adapter import GermlineToG3Adapter

            # Create manager with explicit source (no priority fallback)
            manager = GermlineManager(providers=[genes.source])
            adapter = GermlineToG3Adapter()

            results = []
            for gene_name in genes.genes:
                germline_gene = manager.get_gene_by_name(gene_name, genes.species)
                if germline_gene:
                    g3_dict = adapter.to_g3_format(germline_gene)
                    results.append(g3_dict)
                else:
                    logger.warning(f"Gene {gene_name} not found in {genes.source} database for {genes.species}")

            logger.debug(f"Retrieved {len(results)} genes from germlines module ({genes.source})")
            return results

        # Use G3 API (legacy path)
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
                "gene_curation_source",
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
        input_list: List[Dict[Column, Any]] = input_df.to_dict(orient="records")  # type: ignore

        # can't assign dirrectly so have to append to beat mypy
        for key in input_list:
            reference.data.append(key)
        return reference


class References:
    def __init__(self, default_output_path: Path | str | None = None) -> None:
        self.references: Dict[str, Reference] = {}
        if not default_output_path:
            self.default_output_path = Path(__file__).parent / "../airr/data/germlines"
        else:
            self.default_output_path = Path(default_output_path)
        logger.info(f"Default output path is {self.default_output_path}")
        self.reference_dataframe_path = self.default_output_path / ".references_dataframe.csv.gz"
        logger.info(f"Default dataframe is {self.default_output_path}")

    @property
    def reference_dataframe_path(self) -> Path:
        return self._reference_dataframe_path

    @reference_dataframe_path.setter
    def reference_dataframe_path(self, value: Path) -> None:
        self._reference_dataframe_path = value
        if self._reference_dataframe_path.exists():
            self.reference_dataframe = pd.read_csv(self._reference_dataframe_path, index_col=0)
        else:
            raise FileExistsError(
                f"Reference dataframe does not exist in default path {self._reference_dataframe_path}"
            )

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
        if not self.references:
            return self.reference_dataframe
        for name in self.references:
            # create a single reference object
            _ref: Reference = self.references[name]
            _df: pd.DataFrame = _ref.get_dataframe()

            # because a user could add genes multiple times, lets drop by unique id
            logger.info(f"dropping duplicated: {_df['_id'].duplicated().sum()}")
            _df = _df[~_df["_id"].duplicated()]
            # _df = _df.drop_duplicates("_id")

            # insert name at beggining
            _df.insert(0, "name", name)
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
        self.reference_dataframe = concat_df
        return concat_df

    @staticmethod
    def from_yaml(yaml_path: Optional[Path] = None, use_germlines: bool = False) -> "References":
        """Parse a yaml file into a references file object

        Parameters
        ----------
        yaml_path : Path
            Path to yaml file
        use_germlines : bool, optional
            If True, use local germlines module instead of G3 API. Defaults to False.

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
            reference_object = Reference(use_germlines=use_germlines)

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
        required_columns = [
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

        v_gene_df = database.loc[database["gene_segment"] == "V"].copy()
        missing_columns = [col for col in required_columns if col not in database.columns]
        if missing_columns:
            missing_genes = sorted(v_gene_df["gene"].dropna().unique().tolist())
            raise ValueError("Missing IMGT V-region position columns " f"{missing_columns} for genes: {missing_genes}")

        missing_positions = v_gene_df[v_gene_df[required_columns].isna().any(axis=1)]
        if not missing_positions.empty:
            missing_genes = sorted(missing_positions["gene"].dropna().unique().tolist())
            raise ValueError(
                f"Missing IMGT V-region positions for genes: {missing_genes}. "
                "Ensure IMGT-gapped sequences are available."
            )
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
        """Generate the internal database file for IgBlast.

        Creates combined VDJC FASTA and BLAST database for each reference name.
        The combined file is required for IgBLAST's complete_vdj calculation.

        Parameters
        ----------
        outpath : Path
            The output path to. example -> path/to/output.
            Then the database will dump to path/to/output/{Ig,TCR}/internal_data/{name}/{name}.ndm.imgt

        Notes
        -----
        IgBLAST requires a single file named {name}_V in internal_data/ that
        contains ALL gene segments (V, D, J, C) for proper complete_vdj calculation.
        The NDM file still only includes V gene annotations (framework/CDR regions).
        """
        logger.debug(f"Generating internal annotation file at {outpath}")
        # The internal data file structure goes Ig/internal_path/{name}/

        database = self.get_dataframe()
        for name, group_df in database.groupby("name"):
            # V genes for NDM file (framework/CDR regions)
            v_genes = group_df.loc[group_df["gene_segment"] == "V"].copy()

            # All segments for combined FASTA - IgBLAST needs V+D+J+C in internal_data
            all_segments = group_df.loc[group_df["gene_segment"].isin(["V", "D", "J", "C"])].copy()

            # the species is the actual entity we are using for the annotation, e.g se09 or human
            name_internal_df_path = Path(outpath).joinpath(Path(f"Ig/internal_data/{name}/"))
            if not name_internal_df_path.exists():
                logger.info(f"Creating {name_internal_df_path}")
                name_internal_df_path.mkdir(parents=True)

            # Generate NDM file from V genes only (framework/CDR region annotations)
            if not v_genes.empty:
                # subselect and order
                index_df = v_genes[
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

                # drop anything where there is an na in the annotation index
                index_df = index_df.drop(index_df[index_df.isna().any(axis=1)].index)
                scheme = "imgt"
                internal_annotations_file_path = name_internal_df_path.joinpath(f"{name}.ndm.{scheme}")

                segment = [i.split("|")[-1].split("-")[0][0:4][::-1][:2] for i in index_df["gene"]]
                index_df["segment"] = segment
                index_df["weird_buffer"] = 0
                logger.info(f"Writing to annotation file {internal_annotations_file_path}")
                index_df.to_csv(internal_annotations_file_path, sep="\t", header=False, index=False)

            # Build combined VDJC BLAST database
            # IgBLAST reads {name}_V from internal_data/ for its internal annotation
            suffix = "V"
            db_outpath = Path(str(name_internal_df_path) + f"/{name}_{suffix}")
            # Pass ALL segments (V+D+J+C) to the BLAST database
            make_blast_db_for_internal(all_segments, db_outpath)
            logger.info(f"Built combined VDJC database with {len(all_segments)} sequences for {name}")

    def _make_hmm_files(self, outpath: Path) -> None:
        """Generate HMM files for renumbering from gapped AA sequences.

        Creates Stockholm alignment files and HMM models from V and J gene
        gapped AA sequences for each name/chain combination.

        Parameters
        ----------
        outpath : Path
            The output path. HMM files will be created in:
            - outpath/stockholms/{name}_{chain}.sto
            - outpath/hmms/{name}_{chain}.hmm
        """
        # Create output directories
        stockholms_dir = Path(outpath) / "stockholms"
        hmms_dir = Path(outpath) / "hmms"
        stockholms_dir.mkdir(parents=True, exist_ok=True)
        hmms_dir.mkdir(parents=True, exist_ok=True)

        # Initialize pyhmmer components
        alphabet = pyhmmer.easel.Alphabet.amino()
        builder = pyhmmer.plan7.Builder(alphabet)
        background = pyhmmer.plan7.Background(alphabet)

        # Get reference dataframe
        database = self.get_dataframe()

        # Group by name (species/reference identifier)
        for name, name_df in database.groupby("name"):
            # Process each chain type (H, K, L)
            for chain in ["H", "K", "L"]:
                # Filter for V and J segments with matching chain type
                # Chain type is at position 2 in gene name (e.g., IGHV1-69*01 -> H)
                # Gene name format: IG{chain}{segment}...  (I=0, G=1, H=2, V=3)
                chain_mask = name_df["gene"].str.get(2) == chain
                segment_mask = name_df["gene_segment"].isin(["V", "J"])
                chain_df = name_df[chain_mask & segment_mask].copy()

                if chain_df.empty:
                    continue

                # Extract gapped AA sequences
                # Handle both pipe-separated (chimeric) and regular gene names
                pairs = []
                for _, row in chain_df.iterrows():
                    gene_name = str(row["gene"])
                    gapped_aa = row.get("imgt.sequence_gapped_aa")

                    # If no gapped AA, try to translate from gapped nucleotide
                    if pd.isna(gapped_aa) or not gapped_aa:
                        gapped_nt = row.get("imgt.sequence_gapped")
                        if not pd.isna(gapped_nt) and gapped_nt:
                            gapped_aa = self._translate_gapped_nt_to_aa(str(gapped_nt))

                    # Skip if still no gapped AA sequence
                    if not gapped_aa:
                        continue

                    # Clean gene name (remove pipe prefix for chimeric)
                    clean_name = gene_name.split("|")[-1] if "|" in gene_name else gene_name
                    pairs.append((clean_name, str(gapped_aa)))

                if not pairs:
                    logger.warning(f"No gapped AA sequences for {name} chain {chain}, skipping HMM")
                    continue

                # Write Stockholm alignment file
                sto_path = stockholms_dir / f"{name}_{chain}.sto"
                self._write_stockholm_file(pairs, name, chain, sto_path)

                # Build HMM from Stockholm
                try:
                    with pyhmmer.easel.MSAFile(
                        sto_path, digital=True, alphabet=alphabet, format="stockholm"
                    ) as msa_file:
                        msa = next(msa_file)
                        hmm, _, _ = builder.build_msa(msa, background)

                    # Save HMM to file
                    hmm_path = hmms_dir / f"{name}_{chain}.hmm"
                    with open(hmm_path, "wb") as f:
                        hmm.write(f)

                    logger.debug(f"Built HMM: {hmm_path}")

                except Exception as e:
                    logger.warning(f"Failed to build HMM for {name} chain {chain}: {e}")

    def _write_stockholm_file(
        self, pairs: List[Tuple[str, str]], name: str, chain: str, sto_path: Path
    ) -> None:
        """Write Stockholm alignment file.

        Parameters
        ----------
        pairs : List[Tuple[str, str]]
            List of (gene_name, gapped_aa_sequence) tuples
        name : str
            Reference name (species identifier)
        chain : str
            Chain type (H, K, or L)
        sto_path : Path
            Output path for Stockholm file
        """
        # Find max sequence length and max name length
        max_seq_len = max(len(seq) for _, seq in pairs)
        max_name_len = max(len(gene_name) for gene_name, _ in pairs)

        lines = ["# STOCKHOLM 1.0", f"#=GF ID {name}_{chain}", ""]

        for gene_name, seq in pairs:
            # Pad sequences to same length with gaps (.)
            padded_seq = seq.ljust(max_seq_len, ".")
            lines.append(f"{gene_name.ljust(max_name_len)}  {padded_seq}")

        # Add reference line (RF) matching the alignment length and terminator
        rf_label = "#=GC RF".ljust(max_name_len + 2)
        lines.append(f"{rf_label}  {'x' * max_seq_len}")
        lines.append("//")

        sto_path.write_text("\n".join(lines))

    def _translate_gapped_nt_to_aa(self, gapped_nt: str) -> Optional[str]:
        """
        Translate IMGT-gapped nucleotide sequence to gapped amino acid.

        IMGT gaps are represented as dots. We preserve gap positions while
        translating the nucleotide codons to amino acids.

        Parameters
        ----------
        gapped_nt : str
            IMGT-gapped nucleotide sequence (dots for gaps)

        Returns
        -------
        str or None
            Gapped amino acid sequence, or None if translation fails
        """
        # Remove gaps to get pure nucleotide sequence
        nt_ungapped = gapped_nt.replace(".", "").replace("-", "")

        # Must be multiple of 3 for translation
        # Truncate to nearest codon boundary
        codon_len = (len(nt_ungapped) // 3) * 3
        if codon_len < 3:
            return None

        nt_for_translation = nt_ungapped[:codon_len]

        try:
            aa_seq = str(Seq(nt_for_translation).translate())
        except Exception:
            return None

        # Now we need to insert gaps at the correct positions
        # IMGT numbering: every 3 nucleotide positions = 1 AA position
        # Gaps in NT sequence need to be mapped to AA positions

        # Build the gapped AA sequence
        aa_gapped_chars: List[str] = []
        aa_pos = 0

        i = 0
        while i < len(gapped_nt) and aa_pos < len(aa_seq):
            char = gapped_nt[i]
            if char in (".", "-"):
                # This is a gap - accumulate gaps until we have 3
                gap_count = 0
                while i < len(gapped_nt) and gapped_nt[i] in (".", "-"):
                    gap_count += 1
                    i += 1
                # For every 3 NT gaps, insert 1 AA gap
                aa_gaps = gap_count // 3
                aa_gapped_chars.extend(["."] * aa_gaps)
            else:
                # This is a nucleotide - consume 3 NTs, output 1 AA
                codon_chars = 0
                while i < len(gapped_nt) and codon_chars < 3:
                    if gapped_nt[i] not in (".", "-"):
                        codon_chars += 1
                    i += 1
                if aa_pos < len(aa_seq):
                    aa_gapped_chars.append(aa_seq[aa_pos])
                    aa_pos += 1

        return "".join(aa_gapped_chars) if aa_gapped_chars else None

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
        return References.from_dataframe(_data)  # type: ignore

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

        # Build HMM files for renumbering
        try:
            self._make_hmm_files(output_path)
            logger.info(f"Generated HMM files {output_path}/hmms/")
        except Exception as e:
            logger.warning(f"HMM generation skipped: {e}")

        self.default_output_path = Path(output_path)
        logger.debug(f"Regenerating frame to {self.reference_dataframe_path}")
        self.reference_dataframe = self.get_dataframe()
        _out = self.default_output_path / ".references_dataframe.csv.gz"
        logger.info(f"Writing out reference dataframe to {self.reference_dataframe_path}")
        self.reference_dataframe.to_csv(_out)
        self.reference_dataframe_path = _out
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
        references.reference_dataframe = dataframe
        return references

    def __repr__(self) -> str:
        return self.get_dataframe().__repr__()
