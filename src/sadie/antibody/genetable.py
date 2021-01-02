import logging
import os

import pandas as pd
from numpy import ndarray

from .exception import AmbiguousGene, BadGene, NoSpecies
from .segment import (
    CDR1AA,
    CDR1NT,
    CDR2AA,
    CDR2NT,
    CDR3AA,
    CDR3NT,
    FrameWork1AA,
    FrameWork1NT,
    FrameWork2AA,
    FrameWork2NT,
    FrameWork3AA,
    FrameWork3NT,
    FrameWork4AA,
    FrameWork4NT,
)

logger = logging.getLogger(__name__)


class GeneTable:
    """
    Gene table base class, parses a gene table into pythonic objects

    Attributes
    ----------
    table_path: path
        path to table in posix file system format

    available_species: list
        list of available species

    gene_table: GeneTable
        gene_table dataframe

    """

    def __init__(self, path):
        """Gene Table constructor

        Parameters
        ----------
        path : path
            The posix file system path to dataframe file
        """
        self.table_path = path

        # Load gene table into memory on first instance
        self.gene_table = self.table_path

    @property
    def available_species(self) -> ndarray:
        """Availble species

        Returns
        -------
        ndarray
           available species list
        """
        return self.gene_table["species"].unique()

    @property
    def available_loci(self) -> ndarray:
        """get available loci of table

        Returns
        -------
        ndarray
            available loci
        """
        return self.gene_table["full"].str[0:3].unique()

    @property
    def table_path(self) -> str:
        """table path

        Returns
        -------
        str
            string table path
        """
        return self._table_path

    @table_path.setter
    def table_path(self, path):
        """Setter for table path

        Parameters
        ----------
        path : path
            path to table

        Raises
        ------
        FileNotFoundError
            raised if file not found
        """
        if not os.path.exists(path):
            raise FileNotFoundError("Can't find Gene Table %s" % path)
        self._table_path = path

    @property
    def gene_table(self) -> pd.DataFrame:
        """Gene table attribute

        Returns
        -------
        pd.Dataframe
            pandas dataframe of gene table
        """
        return self._gene_table

    @gene_table.setter
    def gene_table(self, path):
        """set new gene table

        Parameters
        ----------
        path : path
            path to gene table
        """
        self._gene_table = (
            pd.read_csv(path, index_col=0).fillna("").sort_values(["species", "full"]).reset_index(drop=True)
        )

    def get_gene(self, gene, species) -> pd.Series:
        """returns the row corresponding to the gene and species

        Parameters
        ----------
        gene : str
            the gene to lookup, e.g. "IGHV3-15*01"
        species : str
            the species to lookup, .e.g "human"

        Returns
        -------
        pd.Series
            the row object of the gene species

        Raises
        ------
        NoSpecies
           if no species is found
        BadGene
            if no gene is found
        AmbiguousGene
            if there are multiple alleles for a gene
        """
        if species not in self.available_species:
            raise NoSpecies(species, self.available_species)

        # trim down to just the species of interest
        species_gene_table = self.gene_table.set_index("species").loc[species]
        # Check for valid clone name nome
        if "*" in gene:
            # if the user input an allele, lookup by full name
            species_gene_table = species_gene_table.set_index("full")
            if gene not in species_gene_table.index:
                if gene.split("*")[0] in species_gene_table.reset_index().set_index("gene").index:
                    # try to resolve the gene not being found by trimming allle
                    logger.warning(f"{gene} not found but {gene.split('*')[0]} is found")
                    gene = gene.split("*")[0]
                    species_gene_table = species_gene_table.reset_index().set_index("gene")
                else:
                    # if can't find gene
                    raise BadGene(species, gene, species_gene_table.index)
        else:
            # if the user didn't specify an allele, set lookup by gene
            species_gene_table = species_gene_table.set_index("gene")
            if gene not in species_gene_table.index:
                # if we cant find
                raise BadGene(species, gene, species_gene_table.index)

        # Just get row of interest
        row_lookup = species_gene_table.loc[gene]
        if isinstance(row_lookup, pd.DataFrame):
            # If the user was too general and there is more than one entry, we have to ask for an allele too
            _alternatives = list(row_lookup["full"])
            raise AmbiguousGene(species, gene, _alternatives)

        # return the pd.Series of a single object
        return row_lookup

    def get_available_genes(self, species: str, first_allele=True, locus="all", by="full") -> ndarray:
        """Get available genes given a species

        Parameters
        ----------
        species : str
            lookup by species, ex. 'human'
        first_allele : bool, optional
            only display the first allele, by default True
        locus : str, optional, ['all','IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']
            locus, by default 'all'
        by : str,  ['full','gene']
            lookup by full name or gene name name, by default "full"

        Returns
        -------
        ndarray
           arrayy of available genes

        Raises
        ------
        NoSpecies
            if species is not found
        """
        if species not in self.available_species:
            raise NoSpecies(species, self.available_species)

        if (locus not in list(self.available_loci)) and locus != "all":
            raise ValueError(f"locus not in {self.available_loci} or all")

        # Lookup by species
        _species_lookup = self._gene_table.set_index("species").loc[species]

        # If they specify a locus, just get those
        if locus != "all":
            _locus_lookup = _species_lookup[_species_lookup["full"].str[0:3] == locus]
        else:
            _locus_lookup = _species_lookup

        # if they just want the first allele, groupby and return head
        if first_allele:
            return _locus_lookup.groupby("gene").head(1)[by].to_numpy()
        else:
            return _locus_lookup[by].to_numpy()


class VGeneTable(GeneTable):
    """The V gene table object, should be loaded into memory on module level call

    Examples
    --------
    You can get the full lookup with the getgene method

    >>>v_gene = VGeneTable()
    >>>v_gene.get_gene("IGHV3-15*01","human")
    gene                                                IGHV3-15
    fwr1_nt    GAGGTGCAGCTGGTGGAGTCTGCCGGAGCCTTGGTACAGCCTGGGG...
    cdr1_nt                             GGATTCACTTGCAGTAACGCCTGG
    fwr2_nt    ATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTG...
    cdr2_nt                       ATTAAAAGCAAAGCTAATGGTGGGACAACA
    fwr3_nt    GACTACGCTGCACCTGTGAAAGGCAGATTCACCATCTCAAGAGTTG...
    cdr3_nt                                             ACCACAGA
    fwr1_aa                            EVQLVESAGALVQPGGSLRLSCAAS
    cdr1_aa                                             GFTCSNAW
    fwr2_aa                                    MSWVRQAPGKGLEWVGR
    cdr2_aa                                           IKSKANGGTT
    fwr3_aa               DYAAPVKGRFTISRVDSKNTLYLQMNSLKTEDTAVYYC
    cdr3_aa                                                   TT
    Name: IGHV3-15*01, dtype: object

    Or you can use the series to lookup by attribute

    >>>v_gene = VGeneTable()
    >>>v_gene.get_gene("IGHV3-15*01","human")['cdr2_aa']
    "IKSKTDGGTT"
    """

    def __init__(self):
        # use the file name of this file to find path
        _path = os.path.abspath(__file__)
        _basename = os.path.dirname(_path)
        super().__init__(os.path.join(_basename, "data/VSEGMENT.csv.gz"))


class JGeneTable(GeneTable):
    """The J gene table object, should be loaded into memory on module level call

    Examples
    --------
    You can get the full lookup with the get_gene method

    >>>j_gene = JGeneTable()
    >>>j_gene_table.get_gene("IGHJ6*01","human")
    gene                                    IGHJ6
    cdr3_nt         ATTACTACTACTACTACGGTATGGACGTC
    cdr3_aa                             YYYYYGMDV
    fwr4_nt    TGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG
    fwr4_aa                           WGQGTTVTVSS
    Name: IGHJ6*01, dtype: object

    Or you can use the series to lookup by attribute

    >>>j_gene = JGeneTable()
    >>>j_gene.get_gene("IGHJ6*01","human")['cdr3_nt']
    ATTACTACTACTACTACGGTATGGACGTC
    """

    def __init__(self):
        # use the file name of this file to find path
        _path = os.path.abspath(__file__)
        _basename = os.path.dirname(_path)
        super().__init__(os.path.join(_basename, "data/JSEGMENT.csv.gz"))


# Load these at module time so we only load it once
VTABLE = VGeneTable()
JTABLE = JGeneTable()


class VGene:
    """V gene class for handling the V gene segements and retrieving their sequence

    Examples
    --------
    >>>v_gene_object = antibody.VGene("IGHV3-23*01", "human")
    >>>v_gene_object.cdr1_aa
    GFTFSSYA
    >>>v_gene_object["cdr1_aa"]
    GFTFSSYA
    """

    def __init__(self, name, species, region_assign="imgt"):
        """V gene constructor

        Parameters
        ----------
        name : str
            VGene name, ex. IGHV1-69
        species : str
            The derived species, ex. human
        """
        self._name = name
        self._species = species
        self._row_lookup = VTABLE.get_gene(name, species)

        # nt
        self._fwr1_nt = FrameWork1NT(self.row_lookup["fwr1_nt"].upper())
        self._fwr2_nt = FrameWork2NT(self.row_lookup["fwr2_nt"].upper())
        self._fwr3_nt = FrameWork3NT(self.row_lookup["fwr3_nt"].upper())
        self._cdr1_nt = CDR1NT(self.row_lookup["cdr1_nt"].upper())
        self._cdr2_nt = CDR2NT(self.row_lookup["cdr2_nt"].upper())
        self._cdr3_nt = CDR3NT(self.row_lookup["cdr3_nt"].upper())

        # AA
        self._fwr1_aa = FrameWork1AA(self.row_lookup["fwr1_aa"].upper())
        self._fwr2_aa = FrameWork2AA(self.row_lookup["fwr2_aa"].upper())
        self._fwr3_aa = FrameWork3AA(self.row_lookup["fwr3_aa"].upper())
        self._cdr1_aa = CDR1AA(self.row_lookup["cdr1_aa"].upper())
        self._cdr2_aa = CDR2AA(self.row_lookup["cdr2_aa"].upper())
        self._cdr3_aa = CDR3AA(self.row_lookup["cdr3_aa"].upper())
        self._v_full = (
            self._fwr1_aa.aa
            + self._cdr1_aa.aa
            + self._fwr2_aa.aa
            + self._cdr2_aa.aa
            + self._fwr3_aa.aa
            + self._cdr3_aa.aa
        )

        if region_assign != "imgt":
            raise NotImplementedError(f"Other region definitions {region_assign} not implmented yet")
            # anarci_result = anarci.Anarci(scheme="imgt", region_assign=region_assign).run_single(
            #     self._name, self._v_full
            # )
            # self._fwr1_aa = FrameWork1AA(anarci_result.framework1_aa.upper().replace("-", ""))
            # self._fwr2_aa = FrameWork2AA(anarci_result.framework2_aa.upper().replace("-", ""))
            # self._fwr3_aa = FrameWork3AA(anarci_result.framework3_aa.upper().replace("-", ""))
            # self._cdr1_aa = CDR1AA(anarci_result.cdr1_aa.upper().replace("-", ""))
            # self._cdr2_aa = CDR2AA(anarci_result.cdr2_aa.upper().replace("-", ""))
            # self._cdr3_aa = CDR3AA(anarci_result.cdr3_aa.upper().replace("-", ""))
        # Set the ranges according to the length of each amino acid segment:
        _start = 0
        for segment in [
            self._fwr1_aa,
            self._cdr1_aa,
            self._fwr2_aa,
            self._cdr2_aa,
            self._fwr3_aa,
            self._cdr3_aa,
        ]:
            segment.start_index = _start
            _start = _start + len(segment)

        # Set the ranges according to the length of each segment:
        _start = 0
        for segment in [
            self._fwr1_nt,
            self._cdr1_nt,
            self._fwr2_nt,
            self._cdr2_nt,
            self._fwr3_nt,
            self._cdr3_nt,
        ]:
            segment.start_index = _start
            _start = _start + len(segment)

        self._v_gene_nt = "".join(
            map(
                lambda x: str(x),
                [
                    self._fwr1_nt,
                    self._cdr1_nt,
                    self._fwr2_nt,
                    self._cdr2_nt,
                    self._fwr3_nt,
                    self._cdr3_nt,
                ],
            )
        )

        self._v_gene_aa = "".join(
            map(
                lambda x: str(x),
                [
                    self._fwr1_aa,
                    self._cdr1_aa,
                    self._fwr2_aa,
                    self._cdr2_aa,
                    self._fwr3_aa,
                    self._cdr3_aa,
                ],
            )
        )

    @property
    def name(self):
        return self._name

    @property
    def species(self):
        return self._species

    @property
    def row_lookup(self):
        return self._row_lookup

    @property
    def v_gene_nt(self):
        return self._v_gene_nt

    @property
    def v_gene_aa(self):
        return self._v_gene_aa

    @property
    def fwr1_nt(self):
        return self._fwr1_nt

    @property
    def fwr2_nt(self):
        return self._fwr2_nt

    @property
    def fwr3_nt(self):
        return self._fwr3_nt

    @property
    def cdr1_nt(self):
        return self._cdr1_nt

    @property
    def cdr2_nt(self):
        return self._cdr2_nt

    @property
    def cdr3_nt(self):
        return self._cdr3_nt

    @property
    def fwr1_aa(self):
        return self._fwr1_aa

    @property
    def fwr2_aa(self):
        return self._fwr2_aa

    @property
    def fwr3_aa(self):
        return self._fwr3_aa

    @property
    def cdr1_aa(self):
        return self._cdr1_aa

    @property
    def cdr2_aa(self):
        return self._cdr2_aa

    @property
    def cdr3_aa(self):
        return self._cdr3_aa

    @property
    def segment_dictionary(self):
        return {
            "fwr1_nt": self._fwr1_aa,
            "fwr2_nt": self._fwr2_nt,
            "fwr3_nt": self._fwr3_nt,
            "cdr1_nt": self._cdr1_nt,
            "cdr2_nt": self._cdr2_nt,
            "cdr3_nt": self._cdr3_nt,
            "fwr1_aa": self._fwr1_aa,
            "fwr2_aa": self._fwr2_aa,
            "fwr3_aa": self._fwr3_aa,
            "cdr1_aa": self._cdr1_aa,
            "cdr2_aa": self._cdr2_aa,
            "cdr3_aa": self._cdr3_aa,
        }

    def __repr__(self):
        _rep = """<{}>{}->{}
        fwr1_nt-{}
        fwr1_aa-{}
        cdr1_nt-{}
        cdr1_aa-{}
        fwr2_nt-{}
        fwr2_aa-{}
        cdr2_nt-{}
        cdr2_aa-{}
        fwr3_nt-{}
        fwr3_aa-{}
        cdr3_nt-{}
        cdr3_aa-{}
        """.format(
            self.__class__.__name__,
            self.species,
            self.name,
            self.fwr1_nt,
            self.fwr1_aa,
            self.cdr1_nt,
            self.cdr1_aa,
            self.fwr2_nt,
            self.fwr2_aa,
            self.cdr2_nt,
            self.cdr2_aa,
            self.fwr3_nt,
            self.fwr3_aa,
            self.cdr3_nt,
            self.cdr3_aa,
        )
        return _rep

    def __str__(self):
        return self.__repr__()

    def __getitem__(self, lookup):
        return self.segment_dictionary[lookup]


class JGene:
    """
    Examples
    --------
    >>>j_gene_object = antibody.JGene("IGHJ6*01","human")
    >>>j_gene_object.cdr3_aa
    YYYYYGMD
    >>>j_gene_object["cdr1_aa"]
    YYYYYGMD
    """

    def __init__(self, name, species):
        """J gene constructor

        Parameters
        ----------
        name : str
            JGene name, ex. IGHJ6*01
        species : str
            The derived species, ex. human
        """
        self._name = name
        self._species = species
        self._row_lookup = JTABLE.get_gene(name, species)

        # # Now we can the corresponding AA segments for Germline Gene Segment:
        self._cdr3_nt = CDR3NT(self._row_lookup["cdr3_nt"].upper())
        self._cdr3_aa = CDR3AA(self._row_lookup["cdr3_aa"].upper())
        self._fwr4_nt = FrameWork4NT(self._row_lookup["fwr4_nt"].upper())
        self._fwr4_aa = FrameWork4AA(self._row_lookup["fwr4_aa"].upper())

    @property
    def name(self):
        return self._name

    @property
    def species(self):
        return self._species

    @property
    def row_lookup(self):
        return self._row_lookup

    @property
    def cdr3_nt(self):
        return self._cdr3_nt

    @property
    def cdr3_aa(self):
        return self._cdr3_aa

    @property
    def fwr4_nt(self):
        return self._fwr4_nt

    @property
    def fwr4_aa(self):
        return self._fwr4_aa

    @property
    def segment_dictionary(self):
        return {
            "cdr3_nt": self._cdr3_nt,
            "cdr3_aa": self._cdr3_aa,
            "fwr4_nt": self._fwr4_nt,
            "fwr4_aa": self._fwr4_aa,
        }

    def __repr__(self):
        _rep = """<{}>{}->{}
        cdr3_nt-{}
        cdr3_aa-{}
        fwr4_nt-{}
        fwr4_aa-{}
        """.format(
            self.__class__.__name__,
            self.species,
            self.name,
            self.cdr3_nt,
            self.cdr3_aa,
            self.fwr4_nt,
            self.fwr4_aa,
        )
        return _rep

    def __str__(self):
        return self.__repr__()

    def __getitem__(self, lookup):
        return self.segment_dictionary[lookup]
