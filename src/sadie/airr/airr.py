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
from typing import Generator, List, Union


# third party
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.Interfaces import SequenceIterator
from Bio.SeqRecord import SeqRecord


# package/module level
from sadie.airr.airrtable import AirrTable, LinkedAirrTable
from sadie.airr.igblast import IgBLASTN, GermlineData
from sadie.airr.exceptions import BadIgBLASTExe, BadDataSet, BadRequstedFileType

logger = logging.getLogger("AIRR")
warnings.filterwarnings("ignore", "Partial codon")


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
    >>> air_api.run_multiple(pg9_multiple_seqs)
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
        igblast_exe: Union[Path, str] = "",
        adaptable: bool = False,
        database: str = "imgt",
        v_gene_penalty: int = -1,
        d_gene_penalty: int = -1,
        j_gene_penalty: int = -2,
        allow_vdj_overlap: bool = False,
        correct_indel: bool = True,
        temp_directory: Union[str, Path, None] = None,
    ):
        """Airr constructor

        Parameters
        ----------
        species : str
            the species to annotate against, ex. 'human
        igblast_exe : str, optional
            override sadie package executable has to be in $PATH
        adaptable : bool, optional
            turn on adaptable penalties, by default True
        database : str, optional
            run on custom or imgt database, by default "imgt"
        v_gene_penalty : int, optional
            the penalty for mismatched v gene nt, by default -1
        d_gene_penalty : int, optional
            the penalty for mismatched d gene nt, by default -1
        j_gene_penalty : int, optional
            the penalty for mismatched j gene nt, by default -2
        allow_vdj_overlap : bool, optional
            allow vdj overlap genes, by default False
        correct_indel : bool, optional
            correct the indel gaps in the gemrline_aa vs mature_aa alignments in the airrtable
        temp_directory : Union[str,Path,None], optional
            the temporary working directory, by default uses your enviroments tempdir

        Raises
        ------
        BadSpecies
            If you ask for a species that does not have a reference dataset, ex. robot or database=custom,species=human
        """

        # If the temp directory is passed, it is important to keep track of it so we can delete it at the destructory
        self._create_temp = False

        # the setter handles all the logic behind choosign the correct executables
        self.executable = igblast_exe
        self.igblast = IgBLASTN(self.executable)

        # Properties of airr that will be shared with IgBlast class
        self._v_gene_penalty = v_gene_penalty
        self._d_gene_penalty = d_gene_penalty
        self._j_gene_penalty = j_gene_penalty
        self._allow_vdj_overlap = allow_vdj_overlap
        self.adapt_penalty = adaptable

        # if we set allow_vdj, we have to adjust penalties
        if self._allow_vdj_overlap:
            if self.adapt_penalty:
                logger.warning("allow vdj overlap will be overwritten by adapt_penalty")
            if self._d_gene_penalty != -4:
                self._d_gene_penalty = -4
                logger.warning("Allow V(D)J overlap, d_gene_penalty set to -4")
            if self._j_gene_penalty != -3:
                self._j_gene_penalty = -3
                logger.warning("Allow V(D)J overlap, j_gene_penalty set to -3")

        # if we set adapt_penalty
        if self.adapt_penalty:
            # start with -1 if adaptable pentalty
            if self._allow_vdj_overlap:
                logger.warning("Adaptable pentalty resetting J gene penalty to -1")
            self._v_gene_penalty = -1
            self._j_gene_penalty = -1

        # Properties that will be passed to germline Data Class.
        # Pass theese as private since germline class will handle setter logic
        self._species = species
        self._database = database

        # Check if this requested dataset is available
        _available_datasets = GermlineData.get_available_datasets()
        _chosen_datasets = (species.lower(), database.lower())
        if _chosen_datasets not in _available_datasets:
            raise BadDataSet(_chosen_datasets, _available_datasets)

        # set the germline data
        self.germline_data = GermlineData(self.species, database=database)

        # This will set all the igblast params given the Germline Data class whcih validates them
        self.igblast.igdata = self.germline_data.igdata
        self.igblast.germline_db_v = self.germline_data.v_gene_dir
        self.igblast.germline_db_d = self.germline_data.d_gene_dir
        self.igblast.germline_db_j = self.germline_data.j_gene_dir
        self.igblast.aux_path = self.germline_data.aux_path
        self.igblast.organism = self.species

        # setting penalties
        self.igblast.v_penalty = self._v_gene_penalty
        self.igblast.d_penalty = self._d_gene_penalty
        self.igblast.j_penalty = self._j_gene_penalty
        self.igblast.allow_vdj_overlap = self._allow_vdj_overlap

        # set local instance and igblast temp dir instance
        self.temp_directory = temp_directory
        self.igblast.temp_dir = self.temp_directory

        # set airr specific attributes
        self.correct_indel = correct_indel
        self.liable_seqs = []

        # Init pre run check to make sure everything is good
        # We don't have to do this now as it happens at execution.
        self.igblast._pre_check()

    @property
    def adapt_penalty(self) -> bool:
        return self._adapt_penalty

    @adapt_penalty.setter
    def adapt_penalty(self, adapt_penalty: bool):
        if not isinstance(adapt_penalty, bool):
            raise TypeError(f"Adapt_penalty must be bool, not {type(adapt_penalty)}")
        self._adapt_penalty = adapt_penalty

    @property
    def species(self) -> str:
        return self._species

    @property
    def temp_directory(self) -> Union[None, str, Path]:
        return self._temp_directory

    @temp_directory.setter
    def temp_directory(self, temp_directory: Union[None, str, Path]):
        # if we set the temp diretory, we need to create it
        if temp_directory:
            if not os.path.exists(temp_directory):
                os.makedirs(temp_directory)
            self._temp_directory = temp_directory
        # if not, the tempfile class houses where the users system defaults to store temporary stuff
        else:
            self._temp_directory = tempfile.gettempdir()
            logger.info(f"Temp dir - {self._temp_directory}")

    @property
    def igblast(self) -> IgBLASTN:
        """Get IgBLAST instance

        Returns
        -------
        IgBLASTN
            igblastn instance
        """
        return self._igblast

    @igblast.setter
    def igblast(self, igblast: IgBLASTN):
        if not isinstance(igblast, IgBLASTN):
            raise TypeError(f"{igblast} must be an instance of {IgBLASTN}")
        self._igblast = igblast

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
                    igblastn_path = shutil.which(igblastn_path)
                else:
                    # If it's not check the access
                    _access = os.access(igblastn_path, os.X_OK)
                    raise BadIgBLASTExe(igblastn_path, f"is, this executable? Executable-{_access}")
            else:  # The package igblastn is not working
                logger.warning(
                    f"Can't find igblast executable in {igblastn_path}, with system {system} within package {__package__}. Trying to find system installed hmmer"
                )
                igblastn_path = shutil.which(_executable)
                if not igblastn_path:
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
                    igblastn_path, f"Custom igblastn path is not executable {igblastn_path}, {_access} "
                )
        self._executable = igblastn_path

    def run_single(self, seq_id: str, seq: str, scfv=False) -> Union[AirrTable, LinkedAirrTable]:
        """Run a single string sequence

        Parameters
        ----------
        seq_id : str
           the sequence_id of the string object, ex. "my_sequence"
        seq : str
            The string nucletodide sequence, ex. "CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT"

        Returns
        -------
        Union[AirrTable, ]
            Either a single airrtable for a single chain or an ScFV airrtable
        """
        if not isinstance(seq_id, str):
            raise TypeError(f"seq_id must be instance of str, passed {type(seq_id)}")

        records = [SeqRecord(Seq(seq), id=seq_id, name=seq_id)]
        return self.run_records(records, scfv=scfv)

    def run_dataframe(
        self,
        dataframe: pd.DataFrame,
        seq_id_field: Union[str, int],
        seq_field: Union[str, int],
        scfv=False,
        return_join=False,
    ) -> AirrTable:
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

        if not isinstance(dataframe, pd.DataFrame):
            raise TypeError(f"{dataframe} must be instance of pd.DataFrame")

        if dataframe[seq_id_field].duplicated().any():
            raise TypeError(
                f"{dataframe[dataframe[seq_id_field].duplicated()][seq_id_field]} is duplicated. Nees to be unique"
            )

        def _get_seq_generator():
            for seq_id, seq in zip(
                dataframe.reset_index()[seq_id_field],
                dataframe.reset_index()[seq_field],
            ):
                if not isinstance(seq, str):
                    raise ValueError(f"{seq_id} needs to have string. passed {type(seq)}, drop nan's first")

                yield SeqRecord(id=str(seq_id), name=str(seq_id), description="", seq=Seq(seq))

        if return_join:
            dataframe[seq_id_field] = dataframe[seq_id_field].astype(str)
            _df = self.run_records(_get_seq_generator(), scfv=scfv)
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
    ) -> Union[AirrTable, LinkedAirrTable]:
        """Run Airr annotation on seq records

        Parameters
        ----------
        seqrecords : Union[List[SeqRecord],SequenceIterator]
            A list of sequence records or a SequenceIterator from Bio.SeqIO.parse

        Returns
        -------
        Union[AirrTable, LinkedAirrTable]
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

    def run_fasta(self, file: Path, scfv=False) -> Union[AirrTable, LinkedAirrTable]:
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
        Union[AirrTable, LinkedAirrTable]
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
            if not scfv_airr.empty:
                scfv_airr.insert(2, "species", self.species)
            if self.correct_indel:
                scfv_airr.correct_indel()
            return scfv_airr

        else:
            logger.info(f"Running blast on {file}")
            result = self.igblast.run_file(file)
            logger.info(f"Ran blast on  {file}")
            result.insert(2, "species", self.species)
            result = AirrTable(result)
            if result["liable"].any():
                # If we allow adaption,
                if self.adapt_penalty:
                    # recurse
                    for v_penalty in range(-2, -4, -1):
                        for j_penalty in range(-1, -3, -1):
                            _start_df = result[result["liable"]].copy().reset_index().astype({"index": str})
                            if _start_df.empty:
                                logger.info("Corrected all liabilities")
                                break
                            adaptable_api = Airr(
                                self.species,
                                igblast_exe=self.executable,
                                database=self._database,
                                v_gene_penalty=v_penalty,
                                d_gene_penalty=self._d_gene_penalty,
                                j_gene_penalty=j_penalty,
                                allow_vdj_overlap=False,
                                correct_indel=self.correct_indel,
                                temp_directory=self.temp_directory,
                                adaptable=False,
                            )
                            adapt_results = adaptable_api.run_dataframe(_start_df, "index", "sequence")
                            adapt_results = (
                                adapt_results.rename({"sequence_id": "index"}, axis=1)
                                .astype({"index": int})
                                .set_index("index")
                                .rename({"index": ""})
                            )
                            logger.info(
                                f"{len(adapt_results[~adapt_results['liable']])} corrected with adapting penalties v={v_penalty},j={j_penalty}"
                            )
                            if adapt_results[~adapt_results["liable"]].empty:
                                logger.info("Skipping update...")
                                continue
                            result.update(adapt_results[~adapt_results["liable"]])

                    # finally try to correct liabilites with VDJ overlap
                    _start_df = result[result["liable"]].copy().reset_index().astype({"index": str})
                    if _start_df.empty:
                        logger.info("Corrected all liabilities")
                    else:
                        adaptable_api = Airr(
                            self.species,
                            igblast_exe=self.executable,
                            database=self._database,
                            d_gene_penalty=self._d_gene_penalty,
                            allow_vdj_overlap=True,
                            correct_indel=self.correct_indel,
                            temp_directory=self.temp_directory,
                            adaptable=False,
                        )
                        adapt_results = adaptable_api.run_dataframe(_start_df, "index", "sequence")
                        adapt_results = (
                            adapt_results.rename({"sequence_id": "index"}, axis=1)
                            .astype({"index": int})
                            .set_index("index")
                            .rename({"index": ""})
                        )
                        logger.info(
                            f"{len(adapt_results[~adapt_results['liable']])} corrected with adapting penalties v={v_penalty},j={j_penalty}"
                        )
                        if adapt_results[~adapt_results["liable"]].empty:
                            logger.info("Skipping update...")
                        else:
                            result.update(adapt_results[~adapt_results["liable"]])

        if self.correct_indel:
            result.correct_indel()
        return result

    # private run methods
    def _run_scfv(self, file: Path) -> LinkedAirrTable:
        """An internal method kito run a special scfv execution on paired scfv or other linked chains


        Returns
        -------
             - A joined heavy light airr table
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

        airr_table_a = AirrTable(result_a)
        airr_table_b = AirrTable(result_b)

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
            [result_a[result_a["locus"].isin(["IGK", "IGL"])], result_b[result_b["locus"].isin(["IGK", "IGL"])]]
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
        linked_table = _heavy_airr.merge(_light_airr, suffixes=["_heavy", "_light"], on="sequence_id")
        linked_table = LinkedAirrTable(linked_table)
        return linked_table

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
