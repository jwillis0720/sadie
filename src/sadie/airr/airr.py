"""Main Objects for Interacting with Airr"""
import bz2
import gzip
import logging
import os
import shutil
import tempfile

# Std library
import warnings
from pathlib import Path
from typing import Generator, List, Tuple, Union

import filetype
import pandas as pd

# third party
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..reference import loaded_database
from .airrtable import AirrTable, ScfvAirrTable

# Module level
from .igblast import IgBLASTN, ensure_prefix_to

logger = logging.getLogger("AIRR")

# Get out of here with your partial codon warnigns
warnings.filterwarnings("ignore", "Partial codon")


class Error(Exception):
    """Base class for exceptions in this module."""


class BadDataSet(Error):
    """Exception raised for annotating a bad species

    Attributes:
    """

    def __init__(self, requested_type, accepted_types):
        super().__init__()
        self.requested_type = requested_type
        self.accepted_types = accepted_types

    def __str__(self):
        return "{} dataset, avail datasets{}".format(self.requested_type, sorted(self.accepted_types))


class BadRequstedFileType(Error):
    """Exception raised for unsupported file types

    Attributes:
    """

    def __init__(self, requested_type, accepted_types):
        super().__init__()
        self.requested_type = requested_type
        self.accepted_types = accepted_types

    def __str__(self):
        return "{} file passed, only accepts {}".format(self.requested_type, self.accepted_types)


class GermlineData:
    """
    The germline data paths are extremely cumbersome to workwith. This class will abstract away their paths to make it easier to fold into IgBLAST

    Examples
    --------
    >>> gd = GermlineData('human')
    >>> gd.base_dir
    /Users/jwillis/repos/sadie/airr/data/germlines
    >>> gd.v_gene_dir
    /Users/jwillis/repos/sadie/airr/data/germlines/blastdb/Ig/human/human_V'
    >>> gd.aux_path
    /Users/jwillis/repos/sadie/airr/data/germlines/aux_data/human_gl.aux
    """

    def __init__(
        self,
        species: str,
        database="imgt",
        functional="all",
        receptor="Ig",
    ):
        """

        Parameters
        ----------
        species : str
            The species of interest, e.g. human
        receptor : str, optional
            the receptor type, by default "Ig"
        """
        self.species = species
        self.base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/germlines"))
        self.blast_dir = os.path.join(self.base_dir, f"{database}/{functional}/{receptor}/blastdb/{species}_")
        self.v_gene_dir = self.blast_dir + "V"
        self.d_gene_dir = self.blast_dir + "D"
        self.j_gene_dir = self.blast_dir + "J"
        self.aux_path = os.path.join(self.base_dir, f"{database}/aux_db/{species}_gl.aux")
        self.igdata = os.path.join(self.base_dir, f"{database}/{functional}/{receptor}/")

    @property
    def base_dir(self) -> Path:
        """The base dir

        Returns
        -------
        Path
            The base directory path that contains all the germline data
        """
        return self._base_dir

    @base_dir.setter
    def base_dir(self, directory: str):
        _path = Path(directory)
        if not _path.exists():
            raise FileNotFoundError(f"Base directory, {directory} not found")
        self._base_dir = directory

    @property
    def blast_dir(self) -> Path:
        return self._blast_dir

    @blast_dir.setter
    def blast_dir(self, directory: str):
        # Must be a parent since this is not a valid path yet
        _path = Path(directory).parent
        if not _path.exists():
            raise FileNotFoundError(f"Blast directory, {directory} not found")
        self._blast_dir = directory

    @property
    def v_gene_dir(self) -> Path:
        """The V gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           human_V does not exists but it's the prefix to human_V.nod and other files used by blast
        """
        return self._v_gene_dir

    @v_gene_dir.setter
    def v_gene_dir(self, directory: str):
        _path = Path(directory)
        if not ensure_prefix_to(_path):
            raise FileNotFoundError(f"V gene directory glob, {directory} not found")
        self._v_gene_dir = _path

    @property
    def d_gene_dir(self) -> Path:
        """The D gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           ex: human_D does not exists but it's the prefix to human_D.nod and other files used by blast
        """
        return self._d_gene_dir

    @d_gene_dir.setter
    def d_gene_dir(self, directory: str):
        _path = Path(directory)
        if not ensure_prefix_to(_path):
            warnings.warn(f"D gene directory not found for {self.species}", UserWarning)
        self._d_gene_dir = _path

    @property
    def j_gene_dir(self) -> Path:
        """The J gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           ex: human_J does not exists but it's the prefix to human_j.nod and other files used by blast
        """
        return self._j_gene_dir

    @j_gene_dir.setter
    def j_gene_dir(self, directory: str):
        _path = Path(directory)
        if not ensure_prefix_to(_path):
            raise FileNotFoundError(f"J gene directory glob, {directory} not found")
        self._j_gene_dir = _path

    @property
    def aux_path(self) -> Path:
        """The auxillary data path used to reconstruct CDR3 regions.

        Returns
        -------
        Path
           the fully qualified path to the species auxilary data
           ex:/Users/jwillis/repos/sadie/airr/data/germlines/aux_data/human_gl.aux
        """
        return self._aux_path

    @aux_path.setter
    def aux_path(self, directory: str):
        _path = Path(directory)
        if not _path.exists():
            raise FileNotFoundError(f"J gene directory glob, {directory} not found")
        self._aux_path = _path

    @property
    def igdata(self) -> Path:
        return self._igdata

    @igdata.setter
    def igdata(self, directory: Path):
        _path = Path(directory)
        if not _path.exists():
            raise FileNotFoundError(f"IGDATA, {directory} not found")
        self._igdata = _path

    @staticmethod
    def get_available_datasets() -> list:
        """A static non-instantiated method to get a list of avaialble species with the builtin data

        Returns
        -------
        list
           available datasets (common_name, custom|imgt, functional|all)
        """

        functional = list(filter(lambda x: x["functional"] == "F", loaded_database))
        database_types_all_data = set(map(lambda x: (x["common"], x["source"], "all"), loaded_database))
        database_types_functional_data = set(map(lambda x: (x["common"], x["source"], "functional"), functional))
        return list(database_types_all_data) + list(database_types_functional_data)


class Airr:
    """Immune repertoire and single sequence data AIRR (adaptive immune receptor repertoire)

    Examples
    --------

    Only run on single sequence
    >>> pg9_seq = CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGT
              CCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGT
              GGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCC
              AGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATT
              TTTGTTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTA
              TTATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC
    >>> air_api = Airr("human")
    >>> airr_table = air_api.run_single("PG9", pg9_seq)
    >>> airr_table
    sequence_id                                           sequence locus  stop_codon  vj_in_frame  productive  ...  cdr3_end                                 np1 np1_length  np2 np2_length species
    0  GU272045.1  CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT...   IGH       False         True        True  ...       375  GGCTGGTGGGCCCGACTACCGTAATGGGTACAAC         34  NaN          0   human


    Or run on multiple sequcnes

    >>> pg9_multiple_seqs = list(SeqIO.parse('tests/fixtures/fasta_inputs/PG9_H_multiple.fasta','fasta'))
    >>> air_api = Airr("human")
    >>> air_api.run_multiple(pg9_multiple_seqs)
    sequence_id                                           sequence locus  stop_codon  vj_in_frame  productive  ...  cdr3_end                                 np1 np1_length  np2 np2_length species
    0  GU272045.1  CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT...   IGH       False         True        True  ...       375  GGCTGGTGGGCCCGACTACCGTAATGGGTACAAC         34  NaN          0   human
    1  GU272045.1  CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT...   IGH       False         True        True  ...       375  GGCTGGTGGGCCCGACTACCGTAATGGGTACAAC         34  NaN          0   human


    Or run directly on a file

    >>> air_api = Airr("human")
    >>> air_api.run_file('tests/fixtures/fasta_inputs/PG9_H_multiple.fasta')
    sequence_id                                           sequence locus  stop_codon  vj_in_frame  productive  ...  cdr3_end                                 np1 np1_length  np2 np2_length species
    0  GU272045.1  CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT...   IGH       False         True        True  ...       375  GGCTGGTGGGCCCGACTACCGTAATGGGTACAAC         34  NaN          0   human
    1  GU272045.1  CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT...   IGH       False         True        True  ...       375  GGCTGGTGGGCCCGACTACCGTAATGGGTACAAC         34  NaN          0   human
    """

    def __init__(
        self,
        species: str,
        igblast_exe="",
        adaptable=True,
        functional="functional",
        database="imgt",
        temp_directory=".airr_tmp",
    ):
        """Airr constructor

        Parameters
        ----------
        species : str
            the species to annotate against, ex. 'human
        igblast_exe : str, optional
            the executable. Can be a full path but has to be in $PATH, by default "igblastn"
        functional : str,
            run on only functional genes
        database : str,
           custom or imgt database

        Raises
        ------
        BadSpecies
            If you ask for a species that does not have a reference dataset, ex. robot
        """
        # Init igblast api module
        self._create_temp = False

        # quickly check if we have chosen bad species
        _available_datasets = GermlineData.get_available_datasets()
        _chosen_datasets = (species.lower(), database.lower(), functional.lower())
        if _chosen_datasets not in _available_datasets:
            raise BadDataSet(_chosen_datasets, _available_datasets)

        # if not, proceed
        self.igblast = IgBLASTN()
        self.species = species

        # Can catch bad executables here
        self.igblast.executable = igblast_exe

        # Set germline data and setup igblast
        self.germline_data = GermlineData(species, functional=functional, database=database)
        self.igblast.igdata = self.germline_data.igdata
        self.igblast.germline_db_v = self.germline_data.v_gene_dir
        self.igblast.germline_db_d = self.germline_data.d_gene_dir
        self.igblast.germline_db_j = self.germline_data.j_gene_dir
        self.igblast.aux_path = self.germline_data.aux_path
        self.igblast.organism = species
        self.igblast.temp_dir = temp_directory
        self.temp_directory = temp_directory
        if self.temp_directory:
            if not os.path.exists(temp_directory):
                self._create_temp = True
                os.makedirs(temp_directory)

        # do we try adaptable penalties
        self.adapt_penalty = adaptable
        self._liable_seqs = []

        # Init pre run check to make sure everything is good
        # We don't have to do this now as it happens at execution.
        self.igblast._pre_check()

    def run_single(self, seq_id: str, seq: str, scfv=False) -> Union[AirrTable, ScfvAirrTable]:
        """Run a single string sequence

        Parameters
        ----------
        seq_id : str
           the sequence_id of the string object, ex. "my_sequence"
        seq : str
            The string nucletodide sequence, ex. "CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT"

        Returns
        -------
        Union[AirrTable, ScfvAirrTable]
            Either a single airrtable for a single chain or an ScFV airrtable
        """
        if not isinstance(seq_id, str):
            raise TypeError(f"seq_id must be instance of str, passed {type(seq_id)}")

        if not scfv:
            query = ">{}\n{}".format(seq_id, seq)
            result = self.igblast.run_single(query)
            result.insert(2, "species", self.species)
            result = AirrTable(result)

            # There is liable sequences
            if (result["note"].str.lower() != "liable").all():
                return result
            else:
                self._liable_seqs = set(result[result["note"].str.lower() == "liable"].sequence_id)

                # If we allow adaption,
                if self.adapt_penalty:
                    logger.info("Relaxing penalities to resolve liabilities")
                    _tmp_v = self.igblast.v_penalty.value
                    _tmp_j = self.igblast.j_penalty.value
                    self.igblast.v_penalty = -2
                    self.igblast.j_penalty = -1
                    adaptable_result = self.igblast.run_single(query)
                    self.igblast.v_penalty = _tmp_v
                    self.igblast.j_penalty = _tmp_j
                    adaptable_result.insert(2, "species", self.species)
                    adaptable_result = AirrTable(adaptable_result)

                    # If we shifted from liable, return the adaptable results
                    if (adaptable_result["note"].str.lower() != "liable").all():
                        return adaptable_result
                return result
        else:
            with tempfile.NamedTemporaryFile(dir=self.temp_directory) as tmpfile:
                record = SeqRecord(Seq(seq), id=str(seq_id))
                SeqIO.write(record, tmpfile.name, "fasta")
                _results = self.run_file(tmpfile.name, scfv=True)
            return _results

    def run_dataframe(
        self,
        dataframe: pd.DataFrame,
        seq_id_field: Union[str, int],
        seq_field: Union[str, int],
        scfv=False,
        return_join=False,
    ) -> pd.DataFrame:
        """Pass dataframe and field and run airr.

        Parameters
        ----------
        dataframe : pd.DataFrame
            The input dataframe to run airr on

        seq_field: Union[str,int]
           The field in the dataframe to run airr on

        seq_id_field: Union[str,int]:
            The field that you want the "Sequence ID" in the airr table to correspond to.

        scfv : bool, optional
            if the fasta contains an H+L pair, by default False

        Returns
        -------
        pd.DataFrame
            [description]

        ToDo
        -------
        Default seq_id to be index. But have to account for it being a multi index
        """

        def _get_seq_generator():
            for seq_id, seq in zip(
                dataframe.reset_index()[seq_id_field],
                dataframe.reset_index()[seq_field],
            ):
                yield SeqRecord(id=str(seq_id), name=str(seq_id), description="", seq=Seq(seq))

        if return_join:
            dataframe[seq_id_field] = dataframe[seq_id_field].astype(str)
            _df = self.run_multiple(_get_seq_generator(), scfv=scfv).table
            # convert seq id field to stry stince sequence_id is cast to string
            return dataframe.merge(
                _df,
                left_on=seq_id_field,
                right_on="sequence_id",
            )
        else:
            return self.run_multiple(_get_seq_generator(), scfv=scfv)

    def run_multiple(self, seqrecords: List[SeqRecord], scfv=False) -> Union[AirrTable, ScfvAirrTable]:
        """Run multiple seq records

        Parameters
        ----------
        seqrecords : List[SeqRecord]
            A list of sequence records.

        Returns
        -------
        Union[AirrTable, ScfvAirrTable]
            Either a single airrtable for a single chain or an ScFV airrtable

        Raises
        ------
        TypeError
            if you don't pass a list of sequences
        """
        if not isinstance(seqrecords, (list, SeqIO.FastaIO.FastaIterator, Generator)):
            raise TypeError("seqrecords must be of type {} pased {}".format(list, type(seqrecords)))

        # write to tempfile
        with tempfile.NamedTemporaryFile(suffix=".fasta", dir=self.temp_directory) as temp_fasta:
            SeqIO.write(seqrecords, temp_fasta.name, "fasta")
            logger.info(f"Running tempfile {temp_fasta.name}")
            return self.run_file(temp_fasta.name, scfv=scfv)

    def run_file(self, file: Path, scfv=False) -> Union[AirrTable, Tuple[AirrTable, AirrTable]]:
        """Run airr annotator on a fasta file

        If it contains a scfv linked pair, it will annotate both heavy and light chain

        Parameters
        ----------
        file : Path
            a path to a fasta file. Can be uncompressed, bzip or gzip
        scfv : bool, optional
            if the fasta contains an H+L pair, by default False


        Returns
        -------
        Union[AirrTable, ScfvAirrTable]
            Either a single airrtable for a single chain or an ScFV AirrTable

        Raises
        ------
        BadRequstedFileType
            not a fasta or compressed file type
        """
        if isinstance(file, Path):
            # cast to str
            file = str(file)
        _filetype = filetype.guess(file)
        if _filetype:
            logger.info("Guess File Type is %s ", _filetype.extension)
            if _filetype.extension not in ["gz", "bz2"]:
                raise BadRequstedFileType(_filetype, ["bzip2", "gzip"])
            with tempfile.NamedTemporaryFile(delete=False, dir=self.temp_directory) as tmpfile:
                if _filetype.extension == "gz":
                    logger.info("File type is compressed gzip")
                    file_buffer = gzip.open(file)
                else:
                    logger.info("File type is compressed bzip2")
                    file_buffer = bz2.open(file)
                logger.debug(f"copying {file} to {tmpfile.name}")
                shutil.copyfileobj(file_buffer, tmpfile)
                logger.debug(f"copied {file} to {tmpfile.name}")
                file = tmpfile.name

        if scfv:
            logger.info("scfv file was passed")
            scfv_airr = self._run_scfv(file)
            if not scfv_airr.table.empty:
                scfv_airr.table.insert(2, "species", self.species)
            if _filetype:
                os.remove(file)
            if self._create_temp:
                shutil.rmtree(self.temp_directory)
            return scfv_airr

        else:
            logger.info(f"Running blast on {file}")
            result = self.igblast.run_file(file)
            logger.info(f"Ran blast on  {file}")
            result.insert(2, "species", self.species)
            result = AirrTable(result)
            if (result["note"].str.lower() == "liable").any():
                self._liable_seqs = set(result[result["note"].str.lower() == "liable"].sequence_id)
                # If we allow adaption,
                if self.adapt_penalty:
                    logger.info(f"Relaxing penalities to resolve liabilities for {len(self._liable_seqs)}")
                    _tmp_v = self.igblast.v_penalty.value
                    _tmp_j = self.igblast.j_penalty.value

                    # Set these to adaptable
                    self.igblast.v_penalty = -2
                    self.igblast.j_penalty = -1

                    # Set to false so we can call recursive
                    self.adapt_penalty = False
                    liable_dataframe = result.table[result.table["sequence_id"].isin(self._liable_seqs)]

                    # Will call without adaptive
                    adaptable_results = self.run_dataframe(liable_dataframe, "sequence_id", "sequence")

                    adaptable_not_liable = adaptable_results[adaptable_results["note"] != "liable"]
                    logger.info(f"Corrected {len(adaptable_not_liable)} sequences")
                    airr_table = result.table.set_index("sequence_id")
                    airr_table.update(adaptable_not_liable.set_index("sequence_id"))
                    airr_table = airr_table.reset_index()
                    self.igblast.v_penalty = _tmp_v
                    self.igblast.j_penalty = _tmp_j
                    self.adapt_penalty = True
                    result = AirrTable(airr_table)

            if _filetype:
                os.remove(file)
        return result

    def _run_scfv(self, file: Path) -> ScfvAirrTable:
        """An internal method to run a special scfv execution on paired scfv or other linked chains


        Returns
        -------
            ScfvAirrTable - A joined heavy light airr table
        """
        # Do one round of blast on a file
        result_a = self.igblast.run_file(file)

        # Now take out the results from the input sequence
        remaining_seq = (
            result_a[["sequence", "sequence_alignment"]]
            .fillna("")
            .apply(lambda x: x[0].replace(x[1].replace("-", ""), ""), axis=1)
        )
        # and get those ids
        remaining_id = result_a["sequence_id"]

        # Make some seqeuncing records
        seq_records = [SeqRecord(Seq(x), id=str(name)) for x, name in zip(remaining_seq, remaining_id)]
        with tempfile.NamedTemporaryFile() as tmpfile:
            SeqIO.write(seq_records, tmpfile.name, "fasta")
            # Now run airr again, but this time on the remaining sequencess
            result_b = self.igblast.run_file(tmpfile.name)

        airr_table_a = AirrTable(result_a).table
        airr_table_b = AirrTable(result_b).table

        # since we removed the seqeunce out of result B to run it, lets adjust the numerical columns
        adjuster = airr_table_a["sequence"].str.len() - airr_table_b["sequence"].str.len()
        for column in [
            "v_sequence_start",
            "v_sequence_end",
            "d_sequence_start",
            "d_sequence_end",
            "j_sequence_start",
            "j_sequence_end",
            "cdr1_start",
            "cdr1_end",
            "cdr2_start",
            "cdr2_end",
            "cdr3_start",
            "cdr3_end",
            "fwr1_start",
            "fwr2_start",
            "fwr3_start",
            "fwr4_start",
            "fwr1_end",
            "fwr2_end",
            "fwr3_end",
            "fwr4_end",
        ]:
            result_b.loc[:, column] = result_b[column] + adjuster

        # and also, the sequence should be the whole sequence, not just the subsequence
        result_b.loc[result_a.index, "sequence"] = result_a.sequence

        # Grab the Heavy Chains
        heavy_chain_table = pd.concat([result_a[result_a["locus"] == "IGH"], result_b[result_b["locus"] == "IGH"]])

        # Grab the Light Chains out of the set
        light_chain_table = pd.concat(
            [
                result_a[result_a["locus"].isin(["IGK", "IGL"])],
                result_b[result_b["locus"].isin(["IGK", "IGL"])],
            ]
        )

        # this ia a bit of an edge case but if eithere of the two chains are empty, we can fill it with
        # the sequence ids of the other
        if heavy_chain_table.empty:
            heavy_chain_table["sequence_id"] = light_chain_table["sequence_id"]
            heavy_chain_table["sequence"] = light_chain_table["sequence"]
        if light_chain_table.empty:
            light_chain_table["sequence_id"] = heavy_chain_table["sequence_id"]
            light_chain_table["sequence"] = heavy_chain_table["sequence"]
        # sometimes after removing the heavy or light chain, the matcher will find that same locus again, so we have to get uniques
        # so the best match will be the top one
        heavy_chain_table = heavy_chain_table.groupby(["sequence_id", "sequence"]).head(1)
        light_chain_table = light_chain_table.groupby(["sequence_id", "sequence"]).head(1)
        _heavy_airr = AirrTable(heavy_chain_table.reset_index(drop=True))
        _light_airr = AirrTable(light_chain_table.reset_index(drop=True))
        return ScfvAirrTable(_heavy_airr, _light_airr)

    @staticmethod
    def get_available_datasets() -> list:
        return GermlineData.get_available_datasets()

    @staticmethod
    def get_available_species() -> list:
        """get all species available

        Returns
        -------
        list
            Available species
        """
        return list(set(map(lambda x: x[0], GermlineData.get_available_datasets())))

    def __repr__(self):
        return self.igblast.__repr__()

    def __str__(self):
        return self.__repr__()

    def __del__(self):
        """Destructor"""
        if self._create_temp:
            # If we created a temp directory, let's blow it aways
            shutil.rmtree(self.temp_directory)
