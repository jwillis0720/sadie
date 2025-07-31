"""SADIE Airr module"""
from __future__ import annotations

# Std library
import itertools
import logging
import os
import platform
import shutil
import tempfile
import warnings
from multiprocessing import cpu_count
from pathlib import Path
from types import GeneratorType
from typing import Generator, Iterator, List, Optional, Set, Union

# third party
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.Interfaces import SequenceIterator
from Bio.SeqRecord import SeqRecord

# package/module level
from sadie.airr.airrtable import AirrTable, LinkedAirrTable
from sadie.airr.exceptions import BadDataSet, BadIgBLASTExe, BadRequstedFileType
from sadie.airr.igblast import GermlineData, IgBLASTN
from sadie.reference.reference import References

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
        reference_name: str,
        igblast_exe: Path | str = "",
        adaptable: bool = False,
        v_gene_penalty: int = -1,
        d_gene_penalty: int = -1,  # -2 is the default for NCBI web IgBLAST
        j_gene_penalty: int = -2,
        allow_vdj_overlap: bool = False,
        correct_indel: bool = True,
        temp_directory: Optional[str | Path] = None,
        num_cpus: int = -1,
        receptor: str = "Ig",
        scheme: str = "imgt",
        references: Optional[References] = None,
        debug: bool = False,
        # Additional IgBLAST parameters with SADIE defaults
        num_alignments_v: int = 3,
        num_alignments_d: int = 3,
        num_alignments_j: int = 3,
        extend_align5end: bool = True,
        extend_align3end: bool = True,
        min_d_match: int = 5,
        word_size: int = 5,
        gap_open: int = 5,
        gap_extend: int = 2,
        coerce: bool = False,
    ):
        """Airr constructor

        Note: SADIE uses the same penalty parameters as igblastn (V=-1, D=-1, J=-2),
        maintains the same number of alignments (3 for V/D/J), and keeps the minimum D matches at the default 5 nucleotides.
        The main differences are that SADIE still enables both 5' and 3' end extension by default (igblastn: disabled),
        adds explicit gap penalties (open=5, extend=2), and includes the coerce option for handling missing allele annotations.
        SADIE also retains its adaptive penalty adjustment feature when sequences fail initial annotation.

        Parameters
        ----------
        reference_name : str | Reference
            the reference name to run annotate against. These can be composed of IMGT or custom database genes, ex 'human'
            or you can pass the Reference object from a custom database.
        igblast_exe : str
            override sadie package executable has to be in $PATH
        adaptable : bool
            turn on adaptable penalties, by default True
        v_gene_penalty : int
            the penalty for mismatched v gene nt, by default -1
        d_gene_penalty : int
            the penalty for mismatched d gene nt, by default -1
        j_gene_penalty : int
            the penalty for mismatched j gene nt, by default -2
        allow_vdj_overlap : bool
            allow vdj overlap genes, by default False
        correct_indel : bool, optional
            correct the indel gaps in the gemrline_aa vs mature_aa alignments in the airrtable
        temp_directory : str|Path|None
            the temporary working directory, by default uses your enviroments tempdir
        references: Optional[References] = None
            A refernces class with custom references in it. If None, will default to SADIE shipped references
        debug : bool
            if True, print the IgBLAST command before execution, by default False
        num_alignments_v : int
            Number of V gene alignments to show, by default 3
        num_alignments_d : int
            Number of D gene alignments to show, by default 3
        num_alignments_j : int
            Number of J gene alignments to show, by default 3
        extend_align5end : bool
            Extend alignment at 5' end, by default True
        extend_align3end : bool
            Extend alignment at 3' end, by default True
        min_d_match : int
            Minimum D gene nucleotide matches, by default 5
        word_size : int
            Word size for alignment, by default 5
        gap_open : int
            Gap opening penalty, by default 5
        gap_extend : int
            Gap extension penalty, by default 2
        coerce : bool
            Accept the highest scored allele that exists in the auxiliary files when exact match not found, by default False
        """

        # If the temp directory is passed, it is important to keep track of it so we can delete it at the destructory
        self._create_temp = False
        self.references = references
        self.debug = debug

        # Store IgBLAST parameters
        self.num_alignments_v = num_alignments_v
        self.num_alignments_d = num_alignments_d
        self.num_alignments_j = num_alignments_j
        self.extend_align5end = extend_align5end
        self.extend_align3end = extend_align3end
        self.min_d_match = min_d_match
        self.word_size = word_size
        self.gap_open = gap_open
        self.gap_extend = gap_extend

        # check for package executable, then path igblastn in path
        if not igblast_exe:
            try:
                self.executable = self._get_package_igblast_exe()
            except BadIgBLASTExe:
                logger.error(f"Could not find igblast executable {igblast_exe}, searching path")
                try:
                    # cast this to string for type checking
                    self.executable = Path(str(shutil.which("igblastn")))
                except BadIgBLASTExe as e:
                    logger.error(f"package igblast and path igblast not working: {os.environ['PATH']}")
                    raise e
            logger.debug(f"Using igblastn from {self.executable}")
        else:
            self.executable = Path(igblast_exe)

        # give executable to IgBLASTN
        self.igblast = IgBLASTN(self.executable)
        self.igblast.debug = debug

        # Apply IgBLAST parameters
        self.igblast.num_v = num_alignments_v
        self.igblast.num_d = num_alignments_d
        self.igblast.num_j = num_alignments_j
        self.igblast.extend_5 = extend_align5end
        self.igblast.extend_3 = extend_align3end
        self.igblast.min_d_match = min_d_match
        self.igblast.word_size = word_size
        self.igblast.gap_open = gap_open
        self.igblast.gap_extend = gap_extend

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
                warnings.warn("allow vdj overlap will be overwritten by adapt_penalty", UserWarning)
            if self._d_gene_penalty != -4:
                self._d_gene_penalty = -4
                logger.warning("Allow V(D)J overlap, d_gene_penalty set to -4")
            if self._j_gene_penalty != -3:
                self._j_gene_penalty = -3
                logger.warning("Allow V(D)J overlap, j_gene_penalty set to -3")

        # if we set adapt_penalty
        if self.adapt_penalty:
            # start with -1 if adaptable pentalty
            if self._allow_vdj_overlap and self._j_gene_penalty != -1:
                logger.warning("Adaptable pentalty resetting J gene penalty to -1")
                warnings.warn("Adaptable penalty resetting J gene penalty to -1", UserWarning)
            self._v_gene_penalty = -1
            self._j_gene_penalty = -1

        # set temp_directory so it will never be None
        if not temp_directory:
            if "TMPDIR" in os.environ:
                temp_directory = Path(os.environ["TMPDIR"])
            else:
                temp_directory = Path(tempfile.gettempdir()) / "airr/"

        # set igblast temp after instance temp
        self.temp_directory = Path(temp_directory)
        self.igblast.temp_dir = self.temp_directory

        # Properties that will be passed to germline Data Class.
        # Pass theese as private since germline class will handle setter logic
        self._name = reference_name
        self.scheme = scheme  # Store scheme as instance attribute
        if isinstance(references, References):
            _custom_avail = list(references.get_dataframe()["name"].unique())
            if self.name not in _custom_avail:
                raise BadDataSet(self.name, _custom_avail)

            # always make
            out_data_path = references.make_airr_database(self.temp_directory / "germlines/")

            # set the germline database
            self.germline_data = GermlineData(reference_name, receptor, out_data_path)
        else:
            self._name = reference_name

            # Check if this requested dataset is available
            _available_datasets = GermlineData.get_available_datasets()
            if reference_name not in _available_datasets:
                raise BadDataSet(reference_name, list(_available_datasets))

            # set the germline data, None will use default germline
            self.germline_data = GermlineData(self.name, receptor, None, scheme)
        # This will set all the igblast params given the Germline Data class whcih validates them
        self.igblast.igdata = self.germline_data.igdata
        self.igblast.germline_db_v = self.germline_data.v_gene_dir
        self.igblast.germline_db_d = self.germline_data.d_gene_dir
        self.igblast.germline_db_j = self.germline_data.j_gene_dir
        self.igblast.germline_db_c = self.germline_data.c_gene_dir
        self.igblast.aux_path = self.germline_data.aux_path
        self.igblast.organism = self.name

        if num_cpus == -1:
            self.igblast.num_threads = cpu_count()
        else:
            self.igblast.num_threads = num_cpus

        # setting penalties
        self.igblast.v_penalty = self._v_gene_penalty
        self.igblast.d_penalty = self._d_gene_penalty
        self.igblast.j_penalty = self._j_gene_penalty
        self.igblast.allow_vdj_overlap = self._allow_vdj_overlap

        # set airr specific attributes
        self.correct_indel = correct_indel
        self.coerce = coerce
        self.liable_seqs: List[str | int] = []

        # Init pre run check to make sure everything is good
        # We don't have to do this now as it happens at execution.
        self.igblast.pre_check()

    @property
    def adapt_penalty(self) -> bool:
        return self._adapt_penalty

    @adapt_penalty.setter
    def adapt_penalty(self, adapt_penalty: bool) -> None:
        if not isinstance(adapt_penalty, bool):
            raise TypeError(f"Adapt_penalty must be bool, not {type(adapt_penalty)}")
        self._adapt_penalty = adapt_penalty

    @property
    def name(self) -> str:
        return self._name

    @property
    def temp_directory(self) -> Path:
        return self._temp_directory

    @temp_directory.setter
    def temp_directory(self, temp_directory: str | Path) -> None:
        # if we set the temp diretory, we need to create it
        if not os.path.exists(temp_directory):
            os.makedirs(temp_directory)
            self._create_temp = True
            logger.info(f"Temp dir - {temp_directory}")

        self._temp_directory = Path(temp_directory)
        # check for write access
        if not os.access(str(self._temp_directory), os.W_OK):
            raise IOError(f"Temp directory {self._temp_directory} is not writable")

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
    def igblast(self, igblast: IgBLASTN) -> None:
        if not isinstance(igblast, IgBLASTN):
            raise TypeError(f"{igblast} must be an instance of {IgBLASTN}")
        self._igblast = igblast

    @property
    def executable(self) -> Path:
        return self._executable

    @executable.setter
    def executable(self, path: Union[Path, str]) -> None:
        if isinstance(path, str):
            path = Path(path)
        if not path.exists():
            raise BadIgBLASTExe(path, "Exe does not exist")

        # finally check that it's executable
        _access = os.access(path, os.X_OK)
        if not _access:
            raise BadIgBLASTExe(path, f"is, this executable? Executable-{_access}")

        # set path as executable
        self._executable = path

    def run_single(self, seq_id: str, seq: str, scfv: bool = False) -> Union[AirrTable, LinkedAirrTable]:
        """Run a single string sequence

        Parameters
        ----------
        seq_id : str
           the sequence_id of the string object, ex. "my_sequence"
        seq : str
            The string nucletodide sequence, ex. "CAGCGATTAGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGTCGT"

        Returns
        -------
        Union[AirrTable, LinkedAirrTable]
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
        scfv: bool = False,
        return_join: bool = False,
    ) -> Union[AirrTable, pd.DataFrame]:
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
            raise ValueError(
                f"{dataframe[dataframe[seq_id_field].duplicated()][seq_id_field]} is duplicated. Nees to be unique"
            )

        def _get_seq_generator() -> Iterator[SeqRecord]:
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
        self,
        seqrecords: Union[List[SeqRecord], SequenceIterator, Iterator[SeqRecord], Generator[SeqRecord, None, None]],
        scfv: bool = False,
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

    def run_fasta(self, file: Union[Path, str], scfv: bool = False) -> Union[AirrTable, LinkedAirrTable]:
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

        if not Path(file).exists():
            raise FileNotFoundError(file)

        # make sure it's fasta - this will consume the generator but blast has its own fasta parser
        try:
            next(SeqIO.parse(file, "fasta"))
        except Exception:
            raise BadRequstedFileType("", ["fasta"])

        if scfv:
            logger.info("scfv file was passed")
            scfv_airr = self._run_scfv(file)
            if not scfv_airr.empty:
                scfv_airr.insert(2, "reference_name", pd.Series([self.name] * len(scfv_airr)))  #
            if self.correct_indel:
                scfv_airr.correct_indel()
            return scfv_airr

        else:
            logger.info(f"Running blast on {file}")
            result = self.igblast.run_file(Path(file))
            logger.info(f"Ran blast on  {file}")

            # this is worthless since query
            result.insert(2, "reference_name", pd.Series([self.name] * len(result)))
            result = AirrTable(result)
            result["v_penalty"] = self._v_gene_penalty
            result["d_penalty"] = self._d_gene_penalty
            result["j_penalty"] = self._j_gene_penalty

            # Apply coercion if enabled
            if self.coerce:
                result = self._apply_allele_coercion(result)
            if result["liable"].any():
                # If we allow adaption,
                if self.adapt_penalty:
                    # recurse
                    for v_penalty in range(-2, -4, -1):
                        for j_penalty in range(-1, -3, -1):
                            _start_df = pd.DataFrame(result[result["liable"]].copy().reset_index()).astype(
                                {"index": str}
                            )
                            if _start_df.empty:
                                logger.info("Corrected all liabilities")
                                break
                            adaptable_api = Airr(
                                self.name,
                                igblast_exe=self.executable,
                                v_gene_penalty=v_penalty,
                                d_gene_penalty=self._d_gene_penalty,
                                j_gene_penalty=j_penalty,
                                allow_vdj_overlap=False,
                                correct_indel=self.correct_indel,
                                temp_directory=self.temp_directory,
                                adaptable=False,
                                references=self.references,
                                debug=self.debug,
                                # Pass through IgBLAST parameters
                                num_alignments_v=self.num_alignments_v,
                                num_alignments_d=self.num_alignments_d,
                                num_alignments_j=self.num_alignments_j,
                                extend_align5end=self.extend_align5end,
                                extend_align3end=self.extend_align3end,
                                min_d_match=self.min_d_match,
                                word_size=self.word_size,
                                gap_open=self.gap_open,
                                gap_extend=self.gap_extend,
                                coerce=self.coerce,
                            )
                            adapt_results = pd.DataFrame(adaptable_api.run_dataframe(_start_df, "index", "sequence"))
                            adapt_results = pd.DataFrame(
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
                    _start_df = pd.DataFrame(result[result["liable"]].copy().reset_index()).astype({"index": str})
                    if _start_df.empty:
                        logger.info("Corrected all liabilities")
                    else:
                        adaptable_api = Airr(
                            self.name,
                            igblast_exe=self.executable,
                            d_gene_penalty=self._d_gene_penalty,
                            allow_vdj_overlap=True,
                            correct_indel=self.correct_indel,
                            temp_directory=self.temp_directory,
                            adaptable=False,
                            references=self.references,
                            debug=self.debug,
                            # Pass through IgBLAST parameters
                            num_alignments_v=self.num_alignments_v,
                            num_alignments_d=self.num_alignments_d,
                            num_alignments_j=self.num_alignments_j,
                            extend_align5end=self.extend_align5end,
                            extend_align3end=self.extend_align3end,
                            min_d_match=self.min_d_match,
                            word_size=self.word_size,
                            gap_open=self.gap_open,
                            gap_extend=self.gap_extend,
                            coerce=self.coerce,
                        )
                        adapt_results = adaptable_api.run_dataframe(_start_df, "index", "sequence")
                        adapt_results = (
                            pd.DataFrame(adapt_results.rename({"sequence_id": "index"}, axis=1))
                            .astype({"index": int})
                            .set_index("index")
                            .rename({"index": ""})
                        )
                        logger.info(
                            f"{len(adapt_results[~adapt_results['liable']])} corrected with adapting penalties v={v_penalty},j={j_penalty}"
                        )
                        if not adapt_results[~adapt_results["liable"]].empty:
                            result.update(adapt_results[~adapt_results["liable"]])

        if self.correct_indel:
            result.correct_indel()

        return result

    # private run methods
    def _run_scfv(self, file: Union[Path, str]) -> LinkedAirrTable:
        """An internal method kito run a special scfv execution on paired scfv or other linked chains


        Returns
        -------
             - A joined heavy light airr table
        """
        # Do one round of blast on a file
        result_a = self.igblast.run_file(Path(file))

        # Now take out the results from the input sequence
        remaining_seq = (
            result_a[["sequence", "sequence_alignment"]]
            .fillna("")
            .apply(lambda x: x[0].replace(x[1].replace("-", ""), ""), axis=1)
        )
        # and get those ids
        remaining_id = result_a["sequence_id"]

        # Make some seqeuncing records
        seq_records: List[SeqRecord] = [SeqRecord(Seq(x), id=str(name)) for x, name in zip(remaining_seq, remaining_id)]
        with tempfile.NamedTemporaryFile() as tmpfile:
            SeqIO.write(seq_records, tmpfile.name, "fasta")
            # Now run airr again, but this time on the remaining sequencess
            result_b: pd.DataFrame = self.igblast.run_file(tmpfile.name)

        airr_table_a = AirrTable(result_a)
        airr_table_b = AirrTable(result_b)

        # since we removed the seqeunce out of result B to run it, lets adjust the numerical columns
        adjuster = airr_table_a["sequence"].str.len() - airr_table_b["sequence"].str.len()  # type: ignore
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
        result_a_heavy = result_a[result_a["locus"] == "IGH"]
        result_b_heavy = result_b[result_b["locus"] == "IGH"]
        heavy_chain_table = pd.concat([result_a_heavy, result_b_heavy])

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
        linked_table = _heavy_airr.merge(_light_airr, suffixes=("_heavy", "_light"), on="sequence_id")
        linked_table = LinkedAirrTable(linked_table)
        return linked_table

    def _get_package_igblast_exe(self) -> Path:
        """Get package igblastn"""

        # this patch allows us to mock a non-existant igblast that is not in package
        if "IGBLASTN_MONEKY_PATCH" in os.environ:
            executable = os.environ["IGBLASTN_MONEKY_PATCH"]
        else:
            executable = "igblastn"
        system = platform.system().lower()
        return Path(os.path.dirname(os.path.abspath(__file__))) / "bin" / system / executable

    @staticmethod
    def get_available_datasets() -> Set[str]:
        return GermlineData.get_available_datasets()

    @staticmethod
    def get_available_species() -> List[str]:
        """get all species available

        Returns
        -------
        list
            Available species
        """
        return list(GermlineData.get_available_datasets())

    def __repr__(self) -> str:
        return self.igblast.__repr__()

    def _apply_allele_coercion(self, result: AirrTable) -> AirrTable:
        """Apply allele coercion to use available alleles when exact match not found

        Parameters
        ----------
        result : AirrTable
            The AIRR table with IgBLAST results

        Returns
        -------
        AirrTable
            Updated AIRR table with coerced alleles and comments
        """
        # Read the NDM file to get available alleles
        ndm_path = self.germline_data.igdata / "internal_data" / self.name / f"{self.name}.ndm.{self.scheme}"
        if not ndm_path.exists():
            logger.warning(f"NDM file not found at {ndm_path}, skipping coercion")
            return result

        # Load available alleles from NDM file
        available_alleles = set()
        with open(ndm_path, "r") as f:
            for line in f:
                if line.strip():
                    allele = line.split("\t")[0]
                    available_alleles.add(allele)

        # Initialize comments column if it doesn't exist
        if "comments" not in result.columns:
            result["comments"] = ""

        # Process each row
        for idx, row in result.iterrows():
            # Process V, D, and J calls
            for gene_type in ["v", "d", "j"]:
                call_field = f"{gene_type}_call"
                if pd.isna(row[call_field]) or not row[call_field]:
                    continue

                # Get all called alleles with scores
                calls = str(row[call_field]).split(",")
                if not calls:
                    continue

                # Check if top scoring allele is in available set
                top_allele = calls[0]
                if top_allele in available_alleles:
                    continue  # No coercion needed

                # Find best available allele
                gene_family = top_allele.split("*")[0] if "*" in top_allele else top_allele
                best_available = None

                # Look through all called alleles to find one that's available
                for called_allele in calls:
                    if called_allele in available_alleles:
                        best_available = called_allele
                        break

                # If no exact match found, look for same gene family
                if not best_available:
                    for avail_allele in available_alleles:
                        if avail_allele.startswith(gene_family + "*"):
                            best_available = avail_allele
                            break

                # Apply coercion if we found an available allele
                if best_available and best_available != top_allele:
                    # Get scores if available
                    score_field = f"{gene_type}_score"
                    top_score = row.get(score_field, "unknown")

                    # Update the call to use the available allele
                    new_calls = [best_available] + [c for c in calls if c != best_available]
                    result.at[idx, call_field] = ",".join(new_calls)

                    # Add comment about the coercion
                    comment = f"FR/CDR boundaries derived from {best_available} due to missing annotation for top-scoring {top_allele}"
                    if top_score != "unknown":
                        comment += f" (score: {top_score})"

                    # Add to comments field
                    existing_comment = result.at[idx, "comments"]
                    if existing_comment and not pd.isna(existing_comment) and existing_comment:
                        result.at[idx, "comments"] = f"{existing_comment}; {comment}"
                    else:
                        result.at[idx, "comments"] = comment

                    logger.info(f"Coerced {top_allele} to {best_available} for sequence {row.get('sequence_id', idx)}")

        return result

    def __str__(self) -> str:
        return self.__repr__()
