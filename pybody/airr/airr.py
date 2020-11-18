"""Main Objects for Interacting with Airr"""
import bz2
import glob
import gzip
import logging
import os
import shutil
import tempfile

# Std library
import warnings
from pathlib import Path
from typing import List, Tuple, Union

import filetype
import pandas as pd

# third party
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .airrtable import AirrTable, ScfvAirrTable

# Module level
from .igblast import IgBLASTN

logger = logging.getLogger("AIRR")

# Get out of here with your partial codon warnigns
warnings.filterwarnings("ignore", "Partial codon")


class Error(Exception):
    """Base class for exceptions in this module."""


class BadSpecies(Error):
    """Exception raised for annotating a bad species

    Attributes:
    """

    def __init__(self, requested_type, accepted_types):
        super().__init__()
        self.requested_type = requested_type
        self.accepted_types = accepted_types

    def __str__(self):
        return "{} species, avail speciesl{}".format(self.requested_type, self.accepted_types)


class BadRequstedFileType(Error):
    """Exception raised for not finiding the igblast module

    Attributes:
    """

    def __init__(self, requested_type, accepted_types):
        super().__init__()
        self.requested_type = requested_type
        self.accepted_types = accepted_types

    def __str__(self):
        return "{} file passed, only accepts {}".format(self.requested_type, self.accepted_types)


class EagerError(Error):
    def __init__(self):
        super().__init__()

    def __str__(self):
        return "Pleas run Airr first"


class GermlineData:
    """
    The germline data paths are extremely cumbersome to workwith. This class will abstract away their paths to make it easier to fold into IgBLAST

    Examples
    --------
    >>> gd = GermlineData('human')
    >>> gd.base_dir
    /Users/jwillis/repos/pybody/airr/data/germlines
    >>> gd.v_gene_dir
    /Users/jwillis/repos/pybody/airr/data/germlines/blastdb/Ig/human/human_V'
    >>> gd.aux_path
    /Users/jwillis/repos/pybody/airr/data/germlines/aux_data/human_gl.aux
    """

    def __init__(self, species: str, receptor="Ig"):
        """

        Parameters
        ----------
        species : str
            The species of interest, e.g. human
        receptor : str, optional
            the receptor type, by default "Ig"
        """
        self._base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/germlines"))

        blast_dir = os.path.join(
            self._base_dir, "blastdb/{receptor}/{species}/{species}_".format(receptor=receptor, species=species)
        )
        self._v_gene_dir = blast_dir + "V"
        self._d_gene_dir = blast_dir + "D"
        self._j_gene_dir = blast_dir + "J"
        self._aux_path = os.path.join(self.base_dir, f"aux_data/{species}_gl.aux")
        self._available_species = glob.glob(
            os.path.join(
                self.base_dir,
                "aux_data",
            )
            + "/*"
        )

    @property
    def base_dir(self) -> Path:
        """The base dir

        Returns
        -------
        Path
            The base directory path that contains all the germline data
        """
        return self._base_dir

    @property
    def v_gene_dir(self) -> str:
        """The V gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           ex. /Users/jwillis/repos/pybody/airr/data/germlines/blastdb/Ig/human/human_V
           human_V does not exists but it's the prefix to human_V.nod and other files used by blast
        """
        return self._v_gene_dir

    @v_gene_dir.setter
    def v_gene_dir(self, directory):
        self._v_gene_dir = directory

    @property
    def d_gene_dir(self):
        """The D gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           ex. /Users/jwillis/repos/pybody/airr/data/germlines/blastdb/Ig/human/human_D
           human_D does not exists but it's the prefix to human_V.nod and other files used by blast
        """
        return self._d_gene_dir

    @d_gene_dir.setter
    def d_gene_dir(self, directory):
        self._d_gene_dir = directory

    @property
    def j_gene_dir(self):
        return self._j_gene_dir

    @j_gene_dir.setter
    def j_gene_dir(self, directory):
        """The J gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           ex. /Users/jwillis/repos/pybody/airr/data/germlines/blastdb/Ig/human/human_J
           human_D does not exists but it's the prefix to human_V.nod and other files used by blast
        """
        self._j_gene_dir = directory

    @property
    def aux_path(self) -> Path:
        """The auxillary data path used to reconstruct CDR3 regions.

        Returns
        -------
        Path
           the fully qualified path to the species auxilary data
           ex:/Users/jwillis/repos/pybody/airr/data/germlines/aux_data/human_gl.aux
        """
        return self._aux_path

    @aux_path.setter
    def aux_path(self, directory):
        if not os.path.exists(directory):
            raise FileNotFoundError(directory)
        self._aux_path = directory

    @property
    def available_species(self) -> list:
        """list of available species

        Returns
        -------
        list
            avail species
        """
        return self._available_species

    @staticmethod
    def get_available_species() -> list:
        """A static non-instantiated method to get a list of avaialble species with the builtin data

        Returns
        -------
        list
           available species
        """
        base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/germlines/aux_data/"))
        globs = [os.path.basename(i).split("_")[0] for i in glob.glob(base_path + "/*")]
        return globs


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

    def __init__(self, species: str, igblast_exe="igblastn", temp_directory="."):
        """Airr constructor

        Parameters
        ----------
        species : str
            the species to annotate against, ex. 'human
        igblast_exe : str, optional
            the executable. Can be a full path but has to be in $PATH, by default "igblastn"

        Raises
        ------
        BadSpecies
            If you ask for a species that does not have a reference dataset, ex. robot
        """
        if species.lower() not in GermlineData.get_available_species():
            raise BadSpecies(species, GermlineData.get_available_species())
        # Init igblast api module
        self.igblast = IgBLASTN()
        self.species = species

        # Can catch bad executables here
        self.igblast.executable = igblast_exe

        # Set germline data and setup igblast
        self.germline_data = GermlineData(species)
        self.igblast.igdata = self.germline_data.base_dir
        self.igblast.germline_db_v = self.germline_data.v_gene_dir
        self.igblast.germline_db_d = self.germline_data.d_gene_dir
        self.igblast.germline_db_j = self.germline_data.j_gene_dir
        self.igblast.aux_path = self.germline_data.aux_path
        self.igblast.organism = species
        self.temp_directory = temp_directory

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
        if not scfv:
            query = ">{}\n{}".format(seq_id, seq)
            result = self.igblast.run_single(query)
            result.loc[:, "species"] = self.species
        else:
            with tempfile.NamedTemporaryFile(dir=self.temp_directory) as tmpfile:
                record = SeqRecord(Seq(seq), id=str(seq_id))
                SeqIO.write(record, tmpfile.name, "fasta")
                return self.run_file(tmpfile.name, scfv=True)
        return AirrTable(result)

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
        if not isinstance(seqrecords, (list, SeqIO.FastaIO.FastaIterator)):
            raise TypeError("seqrecords must be of type {} pased {}".format(list, type(seqrecords)))

        # write to tempfile
        with tempfile.NamedTemporaryFile(suffix="fasta", dir=".") as temp_fasta:
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
        _filetype = filetype.guess(file)
        if _filetype:
            logger.info("Guess File Type is %s ", _filetype.extension)
            with tempfile.NamedTemporaryFile(delete=False, dir=self.temp_directory) as tmpfile:
                if _filetype.extension == "gz":
                    logger.info("File type is compressed gzip")
                    file_buffer = gzip.open(file)
                elif _filetype.extension == "bz2":
                    logger.info("File type is compressed bzip2")
                    file_buffer = bz2.open(file)
                else:
                    raise BadRequstedFileType(_filetype, ["bzip2", "gzip"])
                logger.debug(f"copying {file} to {tmpfile.name}")
                shutil.copyfileobj(file_buffer, tmpfile)
                logger.debug(f"copied {file} to {tmpfile.name}")
                file = tmpfile.name

        if scfv:
            logger.info("scfv file was passed")
            scfv_airr = self._run_scfv(file)
            if not scfv_airr.table.empty:
                scfv_airr.loc[:, "species"] = self.species
            return scfv_airr

        else:
            logger.info(f"Running blast on {file}")
            result = self.igblast.run_file(file)
            logger.info(f"Ran blast on  {file}")
            result.loc[:, "species"] = self.species
            result = AirrTable(result)
        return result

    def _run_multi_file(self, file):
        def run_single_instance(x):
            result = self.igblast.run_file(x)

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
        return ScfvAirrTable(_heavy_airr, _light_airr)

    @staticmethod
    def get_available_species() -> list:
        """get all species available

        Returns
        -------
        list
            Available species
        """
        return GermlineData.get_available_species()

    def __repr__(self):
        return self.igblast.__repr__()

    def __str__(self):
        return self.__repr__()


if __name__ == "__main__":
    api_air = Airr("human")
    at = api_air.run_single(
        "2545",
        """CCGCCCGGCCGCTGCTATGCGGCATCAAGTAGATTGATAAGCCCTTTCTTTGATCCAGTTTGGAAAGTGGGGTCCCATCA
           CGCTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAACCTGAAGATTTTGCAACTTACTA
           CTGTCAACAGAGTTACACATTCCCCCGTTCATTTGGCGGAGGTACCAAGGTGGAGATCAAACGTACGGTTGCCGCTCCTT
           CTGTATTCATATTTCCGCCCTCCGATGAACAGCTTAAATCGGGCACTGCTTCGGTAGTCTGCCTTCTGAATAATTTCTAT
           CCCCGCGAGGCCAAGGTGCAATGGAAAGTCGACAATGCACTGCAAAGTGGAAACTCGCAAGAAAGCGTCACCGAACAGGA
           CAGTAAGATTCCACCTATAGCCTGTCATCGACACTTACCCTGAGTAAGGCTGATTACGAAAAGCACAAGGTTTACGCTTG
           CGAAGTAACTCACCAGGGCCTCTCAAGCCCTGTTACAAAGTCATTTAACAGAGGGGAATGCTAATTCTAGATAATTATCA
           AGGAGACAGTCATAATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCTGCCCAACCAGCCATGGCC
           GAGGTGCAGCTGCTCGAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTG
           TGCAGCCTCTGGATTTACATTTCGTCGTTATGCTATGAGTTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCA
           GCGCCATAAGCGGCTCAGGTGGGTCAACAAAATACGCTGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCC
           AAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTGTATTACTGTGCACGTCGTGCCGGTCG
           TTGGCTGCAGTCTCCTTATTATTATTATGGGATGGATGTATGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAGCAAGTA
           CCAAGGGGCCTTCAGTTTTTCCGCTCGCTCCAAGTTCTAAATCTACTTCCGGTGGAACTGCTGCTCTGGGGTGCCTGGTT
           AAAGACTATTTCCCAGAACCCGTGACTGTAAGTTGGAACAGCGGAGCATTAACCTCAGGAGTGCACACATTCCCGGCCGT
           ATTGCAAAGTTCTGGCCTGTACTCACTCTCTTCTGTTGTAACGGTTCCATCCAGCTCTTTGGGCACCCAAACCTATATAT
           GCAACGTGAATCACAAACCGTCAAACACGAAAGTCGACAAAAAGGTGGAGCCGAAAACTAGTCACCATCACCACCATCAT
           GGCGCATATCCGTATGATGTGCCGGACTATGCTTCTTAGGGCCAGGCCGGCCAGGAGGCTCG""".replace(
            "\n", ""
        ).replace(
            " ", ""
        ),
        scfv=True,
    )
    print(at.to_genbank("test.gb"))
