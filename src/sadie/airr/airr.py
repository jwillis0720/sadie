"""SADIE Airr module"""
# Std library
import itertools
import logging
import os
import platform
import shutil
import tempfile
import warnings
from pathlib import Path
from types import GeneratorType
from typing import Generator, List, Tuple, Union


# third party
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.Interfaces import SequenceIterator
from Bio.SeqRecord import SeqRecord


# package/module level
from ..anarci import Anarci
from ..reference.yaml import YamlRef
from .airrtable import AirrTable, ScfvAirrTable
from .igblast import IgBLASTN, ensure_prefix_to
from .exceptions import BadIgBLASTExe, BadDataSet, BadRequstedFileType

logger = logging.getLogger("AIRR")

# Get out of here with your partial codon warnigns
warnings.filterwarnings("ignore", "Partial codon")


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
        y = YamlRef()
        db_types = []
        for database_type in y.yaml:
            for functional in y.yaml[database_type]:
                for common in y.yaml[database_type][functional]:
                    if (common, database_type, functional) not in db_types:
                        db_types.append((common, database_type, functional))
        return db_types


class Airr:
    """Immune repertoire data using AIRR standardss(adaptive immune receptor repertoire)

    Examples
    --------

    Run AIRR on single sequence
    >>> pg9_seq = CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGTCCCTGAGACTCTCCTGTGCAGCGT
                  CCGGATTCGACTTCAGTAGACAAGGCATGCACTGGGTCCGCCAGGCTCCAGGCCAGGGGCTGGAGTGGGT
                  GGCATTTATTAAATATGATGGAAGTGAGAAATATCATGCTGACTCCGTATGGGGCCGACTCAGCATCTCC
                  AGAGACAATTCCAAGGATACGCTTTATCTCCAAATGAATAGCCTGAGAGTCGAGGACACGGCTACATATT
                  TTGTTGAGAGAGGCTGGTGGGCCCGACTACCGTAATGGGTACAACTATTACGATTTCTATGATGGTTATT
                  ATAACTACCACTATATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCGAGC
    >>> air_api = Airr("human")
    >>> airr_table = air_api.run_single("PG9", pg9_seq)
    >>> airr_table
    sequence_id                                           sequence locus  stop_codon  vj_in_frame  productive  ...  cdr3_end                                 np1 np1_length  np2 np2_length species
    0  GU272045.1  CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT...   IGH       False         True        True  ...       375  GGCTGGTGGGCCCGACTACCGTAATGGGTACAAC         34  NaN          0   human


    Or run on multiple sequcnes

    >>> pg9_multiple_seqs = list(SeqIO.parse('tests/fixtures/fasta_inputs/PG9_H_multiple.fasta','fasta'))
    >>> air_api = Airr("human")
    >>> air_api.run_records(pg9_multiple_seqs)
    sequence_id                                           sequence locus  stop_codon  vj_in_frame  productive  ...  cdr3_end                                 np1 np1_length  np2 np2_length species
    0  GU272045.1  CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT...   IGH       False         True        True  ...       375  GGCTGGTGGGCCCGACTACCGTAATGGGTACAAC         34  NaN          0   human
    1  GU272045.1  CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT...   IGH       False         True        True  ...       375  GGCTGGTGGGCCCGACTACCGTAATGGGTACAAC         34  NaN          0   human


    Or run directly on a file

    >>> air_api = Airr("human")
    >>> air_api.run_fasta('tests/fixtures/fasta_inputs/PG9_H_multiple.fasta')
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
        v_gene_penalty=-1,
        d_gene_panalty=-1,
        j_gene_penalty=-1,
        temp_directory="",
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

        # If the temp directory is passed, it is important to keep track of it so we can delete it at the destructory
        self._create_temp = False

        # quickly check if we have chosen bad species
        _available_datasets = GermlineData.get_available_datasets()
        _chosen_datasets = (species.lower(), database.lower(), functional.lower())
        if _chosen_datasets not in _available_datasets:
            raise BadDataSet(_chosen_datasets, _available_datasets)

        # if not, proceed
        self.igblast = IgBLASTN()
        # Can catch bad executables here
        self.igblast.executable = igblast_exe

        # Properties of airr that will be shared with IgBlast class
        self.v_gene_penalty = v_gene_penalty
        self.d_gene_penalty = d_gene_panalty
        self.j_gene_penalty = j_gene_penalty

        # Properties that will be passed to germline Data Class.
        self.species = species
        self.functional = functional
        self.database = database

        # Set germline data and setup igblast
        self.germline_data = GermlineData(self.species, functional=functional, database=database)

        # This will set all the igblast params given the Germline Data class whcih validates them
        self.igblast.igdata = self.germline_data.igdata
        self.igblast.germline_db_v = self.germline_data.v_gene_dir
        self.igblast.germline_db_d = self.germline_data.d_gene_dir
        self.igblast.germline_db_j = self.germline_data.j_gene_dir
        self.igblast.aux_path = self.germline_data.aux_path
        self.igblast.organism = species
        self.igblast.v_penalty = self.v_gene_penalty
        self.igblast.d_penalty = self.d_gene_penalty
        self.igblast.j_penalty = self.j_gene_penalty

        self.igblast.temp_dir = temp_directory
        self.temp_directory = temp_directory

        # if we set the temp diretory, we need to create it
        if self.temp_directory:
            if not os.path.exists(temp_directory):
                os.makedirs(temp_directory)

        # if not, the tempfile class houses where the users system defaults to store temporary stuff
        else:
            self.temp_directory = tempfile.gettempdir()
            logger.info(f"Temp dir - {self.temp_directory}")

        # do we try adaptable penalties
        self.adapt_penalty = adaptable
        self._liable_seqs = []

        # Init pre run check to make sure everything is good
        # We don't have to do this now as it happens at execution.
        self.igblast._pre_check()

    @property
    def igblast(self) -> IgBLASTN:
        return self.igblast

    @property
    def executable(self) -> Path:
        return self._executable

    @executable.setter
    def executable(self, path: Path):
        _executable = "igblastn"
        # if the user wants us to find the executable
        if not path:
            system = platform.system().lower()
            igblastn_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                f"bin/{system}/{_executable}",
            )
            # check if its
            if os.path.exists(igblastn_path):
                # check if it's executable
                if shutil.which(igblastn_path):
                    self._executable = shutil.which(igblastn_path)
                else:
                    # If it's not check the access
                    _access = os.access(igblastn_path, os.X_OK)
                    raise BadIgBLASTExe(igblastn_path, f"is, this executable? Executable-{_access}")
            else:  # The package igblastn is not working
                logger.warning(
                    f"Can't find igblast executable in {igblastn_path}, with system {system} within package {__package__}. Trying to find system installed hmmer"
                )
                igblastn_path = shutil.which(_executable)
                if igblastn_path:
                    self._executable = igblastn_path
                else:
                    raise BadIgBLASTExe(
                        igblastn_path,
                        f"Can't find igblastn in package {__package__} or in path {os.environ['PATH']}",
                    )
        else:  # User specifed custome path
            logger.debug(f"User passed custom igblastn {path}")
            igblastn_path = shutil.which(path)
            if igblastn_path:
                self._executable = igblastn_path
            else:
                _access = os.access(igblastn_path, os.X_OK)
                raise BadIgBLASTExe(
                    igblastn_path,
                    f"Custom igblastn path is not executable {igblastn_path}, {_access} ",
                )
        self._executable = igblastn_path

    # Run methods below
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
                _results = self.run_fasta(tmpfile.name, scfv=True)
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
            _df = self.run_records(_get_seq_generator(), scfv=scfv).table
            # convert seq id field to stry stince sequence_id is cast to string
            return dataframe.merge(
                _df,
                left_on=seq_id_field,
                right_on="sequence_id",
            )
        else:
            return self.run_records(_get_seq_generator(), scfv=scfv)

    def run_records(
        self, seqrecords: Union[List[SeqRecord], SequenceIterator, Generator, itertools.chain], scfv=False
    ) -> Union[AirrTable, ScfvAirrTable]:
        """Run Airr annotation on seq records

        Parameters
        ----------
        seqrecords : Union[List[SeqRecord],SequenceIterator]
            A list of sequence records or a SequenceIterator from Bio.SeqIO.parse

        Returns
        -------
        Union[AirrTable, ScfvAirrTable]
            Either a single airrtable for a single chain or an ScFV airrtable

        Raises
        ------
        TypeError
            if you don't pass a list of sequences
        """

        # did they pass a sequence type iterator
        is_iterator = isinstance(seqrecords, SequenceIterator)
        is_list_of_seqs = False
        is_generator = isinstance(seqrecords, GeneratorType) or isinstance(seqrecords, itertools.chain)
        if isinstance(seqrecords, List):
            if all([isinstance(x, SeqRecord) for x in seqrecords]):
                is_list_of_seqs = True

        if not any([is_iterator, is_list_of_seqs, is_generator]):
            raise TypeError(
                f"seqrecords must be an instance of {SequenceIterator} or be a list of {SeqRecord} not {type(seqrecords)}"
            )

        # write to tempfile
        with tempfile.NamedTemporaryFile(suffix=".fasta", dir=self.temp_directory) as temp_fasta:
            SeqIO.write(seqrecords, temp_fasta.name, "fasta")
            logger.info("Running AIRR annotation on records")
            logger.debug(f"Running tempfile {temp_fasta.name}")
            results = self.run_fasta(temp_fasta.name, scfv=scfv)
        return results

    def run_fasta(self, file: Path, scfv=False) -> Union[AirrTable, Tuple[AirrTable, AirrTable]]:
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
            not a fasta file
        """
        if isinstance(file, Path):
            # cast to str
            file = str(file)

        if not Path(file).exists:
            raise FileNotFoundError(file)

        # make sure it's fasta - this will consume the generator but blast has its own fasta parser
        try:
            next(SeqIO.parse(file, "fasta"))
        except Exception:
            raise BadRequstedFileType("", "fasta")

        if scfv:
            logger.info("scfv file was passed")
            scfv_airr = self._run_scfv(file)
            if not scfv_airr.table.empty:
                scfv_airr.table.insert(2, "species", self.species)
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

        return result

    # private run methods
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

    @staticmethod
    def run_mutational_analysis(airrtable: AirrTable, scheme: str) -> AirrTable:
        """Run a mutational analysis given a numbering scheme. Returns an AirrTable with added mutational analysis columns

        This method is computationally expensive. So it's a stand alone static method. It will take in an airr table

        Parameters
        ----------
        airrtable : AirrTable
            An AirrTable class input
        scheme : str
            the numbering scheme: ex, 'martin','kabat','imgt','chothia'

        Returns
        -------
        AirrTable
            returns an airrtable with mutation and scheme fields containing the germline mutations

        Raises
        ------
        TypeError
            if input is not an airrtable
        """
        if not isinstance(airrtable, AirrTable):
            raise TypeError(f"{type(airrtable)} must be of type AirrTable")

        if not airrtable.table.index.is_monotonic_increasing:
            raise IndexError(f"{airrtable.table.index} must be monotonic increasing")

        # create anarci api
        logger.info("Running ANARCI on germline alignment")
        anarci_api = Anarci(scheme=scheme, allowed_chain=["H", "K", "L"])
        germline_results_anarci = anarci_api.run_dataframe(
            airrtable.table["germline_alignment_aa"]
            .str.replace("-", "")
            .to_frame()
            .join(airrtable.table["sequence_id"]),
            "sequence_id",
            "germline_alignment_aa",
        )
        logger.info("Running ANARCI on mature alignment")
        mature_results_anarci = anarci_api.run_dataframe(
            airrtable.table["sequence_alignment_aa"]
            .str.replace("-", "")
            .to_frame()
            .join(airrtable.table["sequence_id"]),
            "sequence_id",
            "sequence_alignment_aa",
        )
        logger.info("Getting ANARCI on alignment tables")
        sets_of_lists = [
            set(germline_results_anarci["Id"]),
            set(mature_results_anarci["Id"]),
            set(airrtable["sequence_id"]),
        ]
        sets_of_lists = sorted(sets_of_lists, key=lambda x: len(x))
        common_results = set.intersection(*sets_of_lists)

        logger.info(f"Can run mutational analysis on {len(common_results)} out of {len(airrtable)} results")
        germline_results_anarci = germline_results_anarci.loc[germline_results_anarci["Id"].isin(common_results), :]
        germline_results_anarci_at = germline_results_anarci.get_alignment_table()
        mature_results_anarci = mature_results_anarci.loc[mature_results_anarci["Id"].isin(common_results), :]
        mature_results_anarci_at = mature_results_anarci.get_alignment_table()
        lookup_dataframe = (
            mature_results_anarci_at.drop(["chain_type", "scheme"], axis=1)
            .set_index("Id")
            .transpose()
            .join(
                germline_results_anarci_at.drop(["chain_type", "scheme"], axis=1).set_index("Id").transpose(),
                lsuffix="_mature",
                rsuffix="_germ",
            )
        )
        lookup_dataframe = lookup_dataframe[sorted(lookup_dataframe.columns)].fillna("-")
        mutation_arrays = []
        logger.info(f"Finding mutations on {len(mature_results_anarci)} sequences")
        for x in mature_results_anarci["Id"]:
            germ_tag = x + "_germ"
            mat_tag = x + "_mature"

            # get section of dataframe for only the two we are interested in
            lookup_specific = lookup_dataframe[[germ_tag, mat_tag]]

            # mutation array are all the mutations in a list
            mutation_array = lookup_specific[
                lookup_specific.apply(lambda x: x[0] != x[1] and x[0] != "X", axis=1)
            ].apply(lambda x: x[0] + x.name + x[1], axis=1)
            if mutation_array.empty:
                mutation_array = []
            else:
                mutation_array = mutation_array.to_list()
            mutation_arrays.append(mutation_array)

        mature_results_anarci["mutations"] = mutation_arrays
        return AirrTable(
            airrtable.table.merge(
                mature_results_anarci.rename({"Id": "sequence_id"}, axis=1)[
                    [
                        "sequence_id",
                        "scheme",
                        "mutations",
                    ]
                ],
                on="sequence_id",
            )
        )

    def __repr__(self):
        return self.igblast.__repr__()

    def __str__(self):
        return self.__repr__()
