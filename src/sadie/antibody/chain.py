"""higher level antibody objects"""
import json

# Third Party
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


# Lib level
from .segment import (
    FrameWork1NT,
    FrameWork2NT,
    FrameWork3NT,
    FrameWork4NT,
    FrameWork1AA,
    FrameWork2AA,
    FrameWork3AA,
    FrameWork4AA,
    CDR1NT,
    CDR2NT,
    CDR3NT,
    CDR1AA,
    CDR2AA,
    CDR3AA,
)
from sadie.utility.util import format_alignment
from .genetable import VGene, JGene


class AntibodyChainAA:
    """heavy or light chain object of amino acid sequences

    Examples
    --------

    >>> heavy_chain_aa = antibody.AntibodyChainAA(
    ...   fwr1_aa="QVQLKESGPGLVQPSQTLSLTCTVS",
    ...    cdr1_aa="GLSLTSNS",
    ...    fwr2_aa="VSWIRQPPGKGLEWMGV",
    ...    cdr2_aa="IWSNGGT",
    ...    fwr3_aa="DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
    ...    cdr3_aa="ASIYYYDADYLHWYFDF",
    ...    fwr4_aa="WGPGTMVTVSS",
    ...    v_gene="IGHV2-47",
    ...    j_gene="IGHJ1",
    ...    species="rat")

    Show segmented VDJ recombination of amino acids

    >>> heavy_chain_aa.get_segmented_vdj_aa()
    QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS

    Get segmented alignments between gemrline and mature

    >>> heavy_chain_aa.get_segmented_alignment_aa()
    antibodychain QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC AR--------YYWYFDF WGPGTMVTVSS
    IGHV2-47      ......................... ........ ................. ....... ......E.....N................P........ .SIYYYDADYLH..... ........... IGHJ1"

    Get germline of specifc segments

    >>> heavy_chain_aa.fwr1_aa_germline
    QVQLKESGPGLVQPSQTLSLTCTVS

    >>> print(heavy_chain_aa.fwr1_aa.get_formatted_alignment())
    germline    QVQLKESGPGLVQPSQTLSLTCTVS
    target      .........................
    QVQLKESGPGLVQPSQTLSLTCTVS
    """

    def __init__(
        self,
        name="antibodychain",
        fwr1_aa="",
        cdr1_aa="",
        fwr2_aa="",
        cdr2_aa="",
        fwr3_aa="",
        cdr3_aa="",
        fwr4_aa="",
        v_gene="",
        j_gene="",
        species="",
        leader="",
        tail="",
        region_def="imgt",
    ):
        """Constructor for AntibodyAA object

        Parameters
        ----------

        name : str, optional
           name of object, by default "antibodychain"
        fwr1_aa : str
           framework 1 amino acid sequence
        cdr1_aa : str
            cdr1 amino acid sequence
        fwr2_aa : str
            framework 2 amino acid sequence
        cdr2_aa : str
            cdr 2 amino acid sequence
        fwr3_aa : str
            framework 3 amino acid sequence
        cdr3_aa : str
            cdr3 amino acid sequence
        fwr4_aa : str
            framework 4 amino acid sequence
        v_gene : str
           v gene name, ex. "IGHV3-15"
        j_gene : str
           J gene name, ex. "IGHJ6*01"
        species : str
           species name, ex "human"
        region_def: how to define the region
        """

        # name and amino acid assignments
        self._name = name
        self._fwr1_aa = FrameWork1AA(fwr1_aa)
        self._cdr1_aa = CDR1AA(cdr1_aa)
        self._fwr2_aa = FrameWork2AA(fwr2_aa)
        self._cdr2_aa = CDR2AA(cdr2_aa)
        self._fwr3_aa = FrameWork3AA(fwr3_aa)
        self._cdr3_aa = CDR3AA(cdr3_aa)
        self._fwr4_aa = FrameWork4AA(fwr4_aa)
        self._v_gene = VGene(v_gene, species, region_assign=region_def)
        self._j_gene = JGene(j_gene, species)
        self._leader = leader
        self._tail = tail
        self._species = species
        # chain AA list in right order
        self._chain_aa = [
            self._fwr1_aa,
            self._cdr1_aa,
            self._fwr2_aa,
            self._cdr2_aa,
            self._fwr3_aa,
            self._cdr3_aa,
            self._fwr4_aa,
        ]

        # also make a dictionary for __item__
        self._chain_dict_aa = {
            "fwr1_aa": self._fwr1_aa,
            "cdr1_aa": self._cdr1_aa,
            "fwr2_aa": self._fwr2_aa,
            "cdr2_aa": self._cdr2_aa,
            "fwr3_aa": self._fwr3_aa,
            "cdr3_aa": self._cdr3_aa,
            "fwr4_aa": self._fwr4_aa,
        }

        # Set germline
        self._fwr1_aa.germline = self._v_gene.fwr1_aa
        self._fwr2_aa.germline = self._v_gene.fwr2_aa
        self._fwr3_aa.germline = self._v_gene.fwr3_aa
        self._fwr4_aa.germline = self._j_gene.fwr4_aa

        self._cdr1_aa.germline = self._v_gene.cdr1_aa
        self._cdr2_aa.germline = self._v_gene.cdr2_aa
        self._cdr3_aa.germline = (self._v_gene.cdr3_aa, self._j_gene.cdr3_aa)
        # vdj aa
        self._vdj_aa = "".join([str(i) for i in self._chain_aa])

        # Set the ranges according to the length of each segment:
        _start = len(self.leader)
        for segment in self._chain_aa:
            segment.start_index = _start
            _start = _start + len(segment)

    @property
    def name(self) -> str:
        """Get name of Antibody Object

        Returns
        -------
        str
           name
        """
        return self._name

    @property
    def leader(self) -> str:
        """Get any leading sequences to this antibody

        Returns
        -------
        str
           the leading sequence of this antibody
        """
        return self._leader

    @property
    def tail(self) -> str:
        """Get any tailing sequence to this antibody

        Returns
        -------
        str
           the leading sequence of this antibody
        """
        return self._tail

    @property
    def species(self) -> str:
        """Get species name

        Returns
        -------
        str
            common species name
        """
        return self._species

    @property
    def fwr1_aa(self) -> str:
        """Get framework 1 amino acid sequence

        Returns
        -------
        str
           amino acid sequence
        """
        return self._fwr1_aa

    @property
    def fwr1_aa_germline(self) -> str:
        """Get germline framework 1 amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._fwr1_aa.germline

    @property
    def cdr1_aa(self) -> str:
        """Get cdr1 amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._cdr1_aa

    @property
    def cdr1_aa_germline(self) -> str:
        """Get cdr1 germline amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._cdr1_aa.germline

    @property
    def fwr2_aa(self) -> str:
        """Get framework 2 amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._fwr2_aa

    @property
    def fwr2_aa_germline(self) -> str:
        """Get framework 2 germline amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._fwr2_aa.germline

    @property
    def cdr2_aa(self) -> str:
        """Get cdr2 amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._cdr2_aa

    @property
    def cdr2_aa_germline(self):
        """Get cdr2 2 germline amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._cdr2_aa.germline

    @property
    def fwr3_aa(self):
        """Get framework 3 amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._fwr3_aa

    @property
    def fwr3_aa_germline(self):
        """Get framework 3  germline amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._fwr3_aa.germline

    @property
    def cdr3_aa(self) -> str:
        """Get cdr3  amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._cdr3_aa

    @property
    def cdr3_aa_germline_v(self) -> str:
        """Get cdr3 germline amino acid sequences that overlap with the v segment

        Returns
        -------
        str
            amino acid sequence
        """
        return self._v_gene.cdr3_aa

    @property
    def cdr3_aa_germline_j(self) -> str:
        """Get cdr3 germline amino acid sequences that overlap with the j segment

        Returns
        -------
        str
            amino acid sequence
        """
        return self._j_gene.cdr3_aa

    @property
    def fwr4_aa(self) -> str:
        """Get framework 4 amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._fwr4_aa

    @property
    def fwr4_aa_germline(self) -> str:
        """Get framework 4 germline amino acid sequence

        Returns
        -------
        str
            amino acid sequence
        """
        return self._fwr4_aa.germline

    @property
    def vdj_aa(self) -> str:
        """return the amino acid sequence of the complete vdj


        Returns
        -------
        str
            vdj_aa
        """
        return self._vdj_aa

    def get_segmented_vdj_aa(self) -> str:
        """get segmented VDJ recombination of amino acids


        Returns
        -------
        str
           segmented vdj amino acid sequence

        Examples
        --------

        Show segmented VDJ recombination of amino acids

        >>> heavy_chain_aa.get_segmented_vdj_aa()
        QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS
        """
        return " ".join([i.aa for i in self._chain_aa])

    def get_segmented_alignment_aa(self, line_length=100, ljust=2, maxid=30) -> str:
        """Perhaps the handiest feature, get the segmented alignment between the germline amino acid sequence and the and the input amino
        acid sequecne

        Parameters
        ----------
        line_length : int, optional
            the line length of the format before breaking to next line, by default 100
        ljust : int, optional
            buffer between id and seq, by default 2
        maxid : int, optional
            truncate id after this length, by default 30

        Returns
        -------
        str
           string alignment

        Examples
        --------

        Get segmented alignments between gemrline and mature

        >>> print(heavy_chain_aa.get_segmented_alignment_aa())
        antibodychain  QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC
        IGHV2-47|IGHJ6 ......................... ........ ................. ....... ......E.....N................P........

        antibodychain   AR--------YYWYFDF WGPGTMVTVSS
        IGHV2-47|IGHJ6  .SIYYYDADYLH..... ...........
        """

        # title for the top alignment and bottom alignment
        bottom_title = self.name
        top_title = self._v_gene.name + "|" + self._j_gene.name

        # Which title is logner
        ljust = max(len(top_title), len(bottom_title)) + ljust

        # get the title of the J segment
        _top_alignment = []
        _bottom_alignment = []

        # Go through the segments and align then
        for seg in self._chain_aa:

            template = seg.alignment[0]
            target = seg.alignment[1]

            template_string = template
            target_string = ""
            # After you align the segment, add it to the list
            for t, j in zip(template, target):
                if t == j:
                    target_string += "."
                else:
                    target_string += j
            _top_alignment.append(template_string)
            _bottom_alignment.append(target_string)

        alignment_object = MultipleSeqAlignment(
            [
                SeqRecord(Seq(" ".join(_top_alignment)), id=top_title),
                SeqRecord(Seq(" ".join(_bottom_alignment)), id=bottom_title),
            ]
        )

        return format_alignment(alignment_object, line_length, ljust, maxid)

    def get_json(self, indent=4) -> json:
        """Return json serialized string object

        This is useful for serializing the object to pass around

        Parameters
        ----------
        indent : int, optional
            how much to indent the json string, by default 4

        Returns
        -------
        json
            json encoded object

        Example
        -------

        >>> print(antibodyaa.get_json())
        {
            "name": "antibodychain",
            "fwr1_aa": "QVQLKESGPGLVQPSQTLSLTCTVS",
            "fwr2_aa": "VSWIRQPPGKGLEWMGV",
            "fwr3_aa": "DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
            "fwr4_aa": "WGPGTMVTVSS",
            "cdr1_aa": "GLSLTSNS",
            "cdr2_aa": "IWSNGGT",
            "cdr3_aa": "ASIYYYDADYLHWYFDF",
            "v_gene": "IGHV2-47",
            "species": "rat",
            "j_gene": "IGHJ1"
        }
        """
        dictionary = {
            "name": self.name,
            "fwr1_aa": self._fwr1_aa.__str__(),
            "fwr2_aa": self._fwr2_aa.__str__(),
            "fwr3_aa": self._fwr3_aa.__str__(),
            "fwr4_aa": self._fwr4_aa.__str__(),
            "cdr1_aa": self._cdr1_aa.__str__(),
            "cdr2_aa": self._cdr2_aa.__str__(),
            "cdr3_aa": self._cdr3_aa.__str__(),
            "v_gene": self._v_gene.name,
            "species": self._species,
            "j_gene": self._j_gene.name,
            "leader": self.leader,
            "tail": self.tail,
        }
        return json.dumps(dictionary, indent=indent)

    def to_json(self, file: str):
        """Dump object to file in JSON

        Parameters
        ----------
        file : path
            str file path
        """
        data = json.loads(self.get_json())
        with open(file, "w") as outfile:
            json.dump(data, outfile)

    @staticmethod
    def read_json(file: str) -> "AntibodyChainAA":
        """
        Read json file into AntibodyChainAA object

        Parameters
        ----------
        file : str path
            file path string

        Returns
        -------
        AntibodyChainAA
           AntibodyChainAA object
        """
        return AntibodyChainAA(**json.load(open(file)))

    @staticmethod
    def from_json(json_object: json) -> "AntibodyChainAA":
        """
        Read json string and return AntibodyChainAA object

        Returns
        -------
        AntibodyChainAA
           AntiboydChainAA Object
        """
        return AntibodyChainAA(**json.loads(json_object))

    def __str__(self):
        return self._vdj_aa

    def __repr__(self):
        # _repr = [f"<{self.__class__.__name__}>"]
        _repr = [f"{self._name}"]
        for x in self._chain_aa:
            _repr.append("{} {}-{}:{}".format(x.__class__.__name__, x.start, x.end, str(x)))
        return "\n".join(_repr)

    def __getitem__(self, lookup):
        return self._chain_dict_aa[lookup]

    def __eq__(self, other):
        return all(
            [
                self._chain_dict_aa == other._chain_dict_aa,
                self._v_gene.name == other._v_gene.name,
                self._j_gene.name == other._j_gene.name,
            ]
        )


class AntibodyChainNT(AntibodyChainAA):
    """heavy or light chain object of nucleotide sequences


    Examples
    --------

    >>> chain_nt = antibody.AntibodyChainNT(
    ... fwr1_nt="CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA",
    ... cdr1_nt="GGCCTATCGCTCACCAGCAACTCC",
    ... fwr2_nt="GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG",
    ... cdr2_nt="ATTTGGTCCAACGGTGGCACC",
    ... fwr3_nt="GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC",
    ... cdr3_nt="GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC",
    ... fwr4_nt="TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC",
    ... v_gene="IGHV2-47",
    ... j_gene="IGHJ3",
    ... species="rat")

    Show segmented VDJ recombination of nucleotides

    >>> chain_nt.get_segmented_vdj_nt()
    CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA GGCCTATCGCTCACCAGCAACTCC GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG ATTTGGTCCAACGGTGGCACC GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC

    Show segmented VDJ recombination of aa

    >>> chain_nt.get_segmented_vdj_aa()
    QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS


    Get segmented alignments between gemrline and mature nucleotides

    >>> chain_nt.get_segmented_alignment_nt()
    antibodychain CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT GGGTTATCATTAACCAGCAATAGT GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA ATATGGAGTAATGGAGGCACA GATTATAATTCAGCTATCAAATCCCGACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTTCTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT GCCAGAAA--------------------------ACAATTGGTTTGCTTAC TGGGGCCAAGGCACTCTGGTCACTGTCTCTTCAG
    IGHV2-47      ..G.....G..G......AGC..C......T....G...........A..T..T.....G..A.....C..G..A ..CC....GC.C........CTCC ..C........A..T.....G.....C..A.................T..G ..T...TCC..C..T.....C ..C..C..C..C......G.GAG...CT..TCT....A.C.C........T..ATC...G........G.........TCG..T...C....G..T..G..T........T..C ...TCC.TTTATTACTATGACGCTGACTACCTCC..TGG.AC..C.A..T. .......CC......A....G..C..GAGC..-C IGHJ3

    Get segmented alignments between gemrline and mature aa

    >>> chain_nt.get_segmented_alignment_aa()
    antibodychain QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC AR--------YYWYFDF WGPGTMVTVSS
    IGHV2-47      ......................... ........ ................. ....... ......E.....N................P........ .SIYYYDADYLH..... ........... IGHJ1"

    Get specifc segments for amino acids or nucletodies for mature or germline segments

    >>> chain_nt.fwr1_aa
    FrameWork1AA 1-26: QVQLKESGPGLVQPSQTLSLTCTVSe
    QVQLKESGPGLVQPSQTLSLTCTVS

    >>> chain_nt.fwr1_nt
    FrameWork1NT 1-76: CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA

    >>> chain_nt.fwr1_nt_germline
    CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT

    >>> chain_nt.fwr1_aa_germline
    QVQLKESGPGLVQPSQTLSLTCTVS
    """

    def __init__(
        self,
        name="antibodychain",
        fwr1_nt="",
        cdr1_nt="",
        fwr2_nt="",
        cdr2_nt="",
        fwr3_nt="",
        cdr3_nt="",
        fwr4_nt="",
        v_gene="",
        j_gene="",
        species="",
        leader="",
        tail="",
    ):
        """Constructor for AntibodyNT object

        Parameters
        ----------
        name : str, optional
           name of object, by default "antibodychain"
        fwr1_nt : str
           framework 1 nucleotide sequence
        cdr1_nt : str
            cdr1 nucleotide sequence
        fwr2_nt : str
            framework 2 nucleotide sequence
        cdr2_nt : str
            cdr 2 nucleotide sequence
        fwr3_nt : str
            framework 3 nucleotide sequence
        cdr3_nt : str
            cdr3 nucleotide sequence
        fwr4_nt : str
            framework 4 nucleotide sequence
        v_gene : str
           v gene name, ex. "IGHV3-15"
        j_gene : str
           J gene name, ex. "IGHJ6*01"
        species : str
           species name, ex "human"
        leader : str
            leading nucleotide sequence
        tail : str
            tailing nucleotide sequence
        """
        self._name = name
        self._fwr1_nt = FrameWork1NT(fwr1_nt)
        self._cdr1_nt = CDR1NT(cdr1_nt)
        self._fwr2_nt = FrameWork2NT(fwr2_nt)
        self._cdr2_nt = CDR2NT(cdr2_nt)
        self._fwr3_nt = FrameWork3NT(fwr3_nt)
        self._cdr3_nt = CDR3NT(cdr3_nt)
        self._fwr4_nt = FrameWork4NT(fwr4_nt)
        self._leader = leader
        self._tail = tail
        super().__init__(
            self.name,
            self._fwr1_nt.aa,
            self._cdr1_nt.aa,
            self._fwr2_nt.aa,
            self._cdr2_nt.aa,
            self._fwr3_nt.aa,
            self._cdr3_nt.aa,
            self._fwr4_nt.aa,
            v_gene,
            j_gene,
            species,
            Seq(self._leader).translate(),
            Seq(self._tail).translate(),
        )
        self._chain_nt = [
            self._fwr1_nt,
            self._cdr1_nt,
            self._fwr2_nt,
            self._cdr2_nt,
            self._fwr3_nt,
            self._cdr3_nt,
            self._fwr4_nt,
        ]
        self._chain_dict_nt = {
            "fwr1_nt": self._fwr1_nt,
            "cdr1_nt": self._cdr1_nt,
            "fwr2_nt": self._fwr2_nt,
            "cdr2_nt": self._cdr2_nt,
            "fwr3_nt": self._fwr3_nt,
            "cdr3_nt": self._cdr3_nt,
            "fwr4_nt": self._fwr4_nt,
        }
        # Set germline
        self._fwr1_nt.germline = self._v_gene.fwr1_nt
        self._fwr2_nt.germline = self._v_gene.fwr2_nt
        self._fwr3_nt.germline = self._v_gene.fwr3_nt
        self._fwr4_nt.germline = self._j_gene.fwr4_nt

        self._cdr1_nt.germline = self._v_gene.cdr1_nt
        self._cdr2_nt.germline = self._v_gene.cdr2_nt
        self._cdr3_nt.germline = (self._v_gene.cdr3_nt, self._j_gene.cdr3_nt)

        # V portion must come after germline assignemtn
        self._v_nt = [
            self._fwr1_nt,
            self._cdr1_nt,
            self._fwr2_nt,
            self._cdr2_nt,
            self._fwr3_nt,
            self._cdr3_nt.v_portion,
        ]
        # vdj nt
        self._vdj_nt = "".join([str(i) for i in self._chain_nt])
        self._v_segment = "".join([str(i) for i in self._v_nt])

        # Set the ranges according to the length of each segment:
        _start = len(self._leader)
        for segment in self._chain_nt:
            if not segment:
                continue
            segment.start_index = _start
            _start = _start + len(segment)

    @property
    def vdj_nt(self) -> str:
        return self._vdj_nt

    @property
    def v_region(self) -> str:
        return self._v_segment

    @v_region.setter
    def v_region(self, v_region: str):
        self._v_segment = v_region

    @property
    def leader(self) -> str:
        return self._leader

    @property
    def tail(self) -> str:
        return self._leader

    @property
    def fwr1_nt(self) -> str:
        """Get framework 1 nucleotide sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._fwr1_nt

    @property
    def fwr1_nt_germline(self) -> str:
        """Get framework 1 nucleotide germline sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._fwr1_nt.germline

    @property
    def cdr1_nt(self):
        """Get cdr1 nucleotide sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._cdr1_nt

    @property
    def cdr1_nt_germline(self):
        """Get cdr1 nucleotide germline sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._cdr1_nt.germline

    @property
    def fwr2_nt(self):
        """Get framework 2 nucleotide sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._fwr2_nt

    @property
    def fwr2_nt_germline(self):
        """Get framework 2 nucleotide germline sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._fwr2_nt.germline

    @property
    def cdr2_nt(self):
        """Get cdr2 nucleotide sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._cdr2_nt

    @property
    def cdr2_nt_germline(self) -> str:
        """Get cdr2 nucleotide germline sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._cdr2_nt.germline

    @property
    def fwr3_nt(self) -> str:
        """Get framework 3 nucleotide sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._fwr3_nt

    @property
    def fwr3_nt_germline(self) -> str:
        """Get framework 3 nucleotide germline sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._fwr3_nt.germline

    @property
    def cdr3_nt(self) -> str:
        """Get cdr3 nucleotide sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._cdr3_nt

    @property
    def cdr3_nt_germline_v(self) -> str:
        """Get cdr3  nucleotide germline sequence from the V gene segment

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._v_gene.cdr3_nt

    @property
    def cdr3_nt_germline_j(self):
        """Get cdr3  nucleotide germline sequence from the J gene segment

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._j_gene.cdr3_nt

    @property
    def fwr4_nt(self):
        """Get framework 4 nucleotide sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._fwr4_nt

    @property
    def fwr4_nt_germline(self) -> str:
        """Get framework 4 nucleotide germline sequence

        Returns
        -------
        str
           nucleotide sequence
        """
        return self._j_gene.fwr4_nt

    def get_segmented_vdj_nt(self) -> str:
        """get segmented VDJ recombination of nucletodie


        Returns
        -------
        str
           segmented vdj nucletodie sequence

        Examples
        --------
         Show segmented VDJ recombination of nucletodie

        >>> chain_nt.get_segmented_alignment_nt()
        antibodychain CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT GGGTTATCATTAACCAGCAATAGT GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA ATATGGAGTAATGGAGGCACA GATTATAATTCAGCTATCAAATCCCGACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTTCTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT GCCAGAAA--------------------------ACAATTGGTTTGCTTAC TGGGGCCAAGGCACTCTGGTCACTGTCTCTTCAG
        IGHV2-47      ..G.....G..G......AGC..C......T....G...........A..T..T.....G..A.....C..G..A ..CC....GC.C........CTCC ..C........A..T.....G.....C..A.................T..G ..T...TCC..C..T.....C ..C..C..C..C......G.GAG...CT..TCT....A.C.C........T..ATC...G........G.........TCG..T...C....G..T..G..T........T..C ...TCC.TTTATTACTATGACGCTGACTACCTCC..TGG.AC..C.A..T. .......CC......A....G..C..GAGC..-C IGHJ3
        """
        return " ".join([i.nt for i in self._chain_nt])

    def get_segmented_alignment_nt(self, line_length=100, ljust=2, maxid=30) -> str:
        """Perhaps the handiest feature, get the segmented alignment between the germline nucletoide sequence and the and the input nucletodie sequence

        Parameters
        ----------
        line_length : int, optional
            the line length of the format before breaking to next line, by default 100
        ljust : int, optional
            buffer between id and seq, by default 2
        maxid : int, optional
            truncate id after this length, by default 30

        Returns
        -------
        str
           string alignment

        Examples
        --------
        Get segmented alignments between germline and mature

        >>> chain_nt.get_segmented_alignment_nt()
        antibodychain  CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT GGGTTATCATTAACCAGCAATAGT GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA ATATGGAGTAATGGAGGCACA GATTATAATTCAGCTATCAAATCCCGACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTTCTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT GCCAGAAA--------------------------ACAATTGGTTTGCTTAC TGGGGCCAAGGCACTCTGGTCACTGTCTCTTCAG
        IGHV2-47|IGHJ3 ..G.....G..G......AGC..C......T....G...........A..T..T.....G..A.....C..G..A ..CC....GC.C........CTCC ..C........A..T.....G.....C..A.................T..G ..T...TCC..C..T.....C ..C..C..C..C......G.GAG...CT..TCT....A.C.C........T..ATC...G........G.........TCG..T...C....G..T..G..T........T..C ...TCC.TTTATTACTATGACGCTGACTACCTCC..TGG.AC..C.A..T. .......CC......A....G..C..GAGC..-C"""

        # title for the top alignment and bottom alignment
        bottom_title = self.name
        top_title = self._v_gene.name + "|" + self._j_gene.name

        # Which title is logner
        ljust = max(len(top_title), len(bottom_title)) + ljust

        # get the title of the J segment
        _top_alignment = []
        _bottom_alignment = []

        # Go through the segments and align then
        for seg in self._chain_nt:

            template = seg.alignment[0]
            target = seg.alignment[1]

            template_string = template
            target_string = ""
            # After you align the segment, add it to the list
            for t, j in zip(template, target):
                if t == j:
                    target_string += "."
                else:
                    target_string += j
            _top_alignment.append(template_string)
            _bottom_alignment.append(target_string)

        alignment_object = MultipleSeqAlignment(
            [
                SeqRecord(Seq(" ".join(_top_alignment)), id=str(top_title)),
                SeqRecord(Seq(" ".join(_bottom_alignment)), id=str(bottom_title)),
            ]
        )

        return format_alignment(alignment_object, line_length, ljust, maxid)

    def get_json(self, indent=4) -> json:
        """Return json serialized string object

        This is useful for serializing the object to pass around

        Parameters
        ----------
        indent : int, optional
            how much to indent the json string, by default 4

        Returns
        -------
        json
            json encoded object

        Example
        -------

        >>> print(antibodynt.get_json())
        {
            "name": "antibodychain",
            "fwr1_nt": "CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA",
            "fwr2_nt": "GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG",
            "fwr3_nt": "GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC",
            "fwr4_nt": "TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC",
            "cdr1_nt": "GGCCTATCGCTCACCAGCAACTCC",
            "cdr2_nt": "ATTTGGTCCAACGGTGGCACC",
            "cdr3_nt": "GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC",
            "v_gene": "IGHV2-47",
            "species": "rat",
            "j_gene": "IGHJ3"
        }
        """
        dictionary = {
            "name": self.name,
            "fwr1_nt": self._fwr1_nt.__str__(),
            "fwr2_nt": self._fwr2_nt.__str__(),
            "fwr3_nt": self._fwr3_nt.__str__(),
            "fwr4_nt": self._fwr4_nt.__str__(),
            "cdr1_nt": self._cdr1_nt.__str__(),
            "cdr2_nt": self._cdr2_nt.__str__(),
            "cdr3_nt": self._cdr3_nt.__str__(),
            "v_gene": self._v_gene.name,
            "species": self._species,
            "j_gene": self._j_gene.name,
            "leader": str(self.leader),
            "tail": str(self.tail),
        }
        return json.dumps(dictionary, indent=indent)

    def to_json(self, file: str):
        """Dump object to file in JSON

        Parameters
        ----------
        file : path
            str file path
        """
        data = json.loads(self.get_json())
        with open(file, "w") as outfile:
            json.dump(data, outfile)

    @staticmethod
    def read_json(file: str) -> "AntibodyChainNT":
        """
        Read json file into AntibodyChainNT object

        Parameters
        ----------
        file : str path
            file path string

        Returns
        -------
        AntibodyChainNT
           AntibodyChainNT object
        """
        return AntibodyChainNT(**json.load(open(file)))

    @staticmethod
    def from_json(json_object: json) -> "AntibodyChainNT":
        """
        Read json string and return AntibodyChainNT object

        Returns
        -------
        AntibodyChainNT
           AntiboydChainNT Object
        """
        return AntibodyChainNT(**json.loads(json_object))

    def __str__(self):
        return self._vdj_nt

    def __repr__(self):
        # _repr = [f"<{self.__class__.__name__}>\n{self._name}"]
        _repr = [f"{self._name}"]
        for x in self._chain_nt:
            _repr.append("{} {}-{}:{}".format(x.__class__.__name__, x.start, x.end, str(x)))
        return "\n".join(_repr)

    def __getitem__(self, lookup):
        return self._chain_dict_nt[lookup]

    def __eq__(self, other):
        return all(
            [
                self._chain_dict_nt == other._chain_dict_nt,
                self._chain_dict_aa == other._chain_dict_aa,
                self._v_gene.name == other._v_gene.name,
                self._j_gene.name == other._j_gene.name,
            ]
        )


class HeavyChainAA(AntibodyChainAA):
    """heavy or chain object of amino acid sequences

    Examples
    --------

    >>> heavy_chain_aa = antibody.HeavyChainAA(
    ...   fwr1_aa="QVQLKESGPGLVQPSQTLSLTCTVS",
    ...    cdr1_aa="GLSLTSNS",
    ...    fwr2_aa="VSWIRQPPGKGLEWMGV",
    ...    cdr2_aa="IWSNGGT",
    ...    fwr3_aa="DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
    ...    cdr3_aa="ASIYYYDADYLHWYFDF",
    ...    fwr4_aa="WGPGTMVTVSS",
    ...    v_gene="IGHV2-47",
    ...    j_gene="IGHJ1",
    ...    species="rat")

    Show segmented VDJ recombination of amino acids

    >>> heavy_chain_aa.get_segmented_vdj_aa()
    QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS

    Get segmented alignments between gemrline and mature

    >>> heavy_chain_aa.get_segmented_alignment_aa()
    antibodychain QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC AR--------YYWYFDF WGPGTMVTVSS
    IGHV2-47      ......................... ........ ................. ....... ......E.....N................P........ .SIYYYDADYLH..... ........... IGHJ1"

    Get mature or germline of specifc segments

    >>> heavy_chain_aa.fwr3_aa
    DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC

    >>> heavy_chain_aa.fwr3_aa_germline
    DYNSAIISRLSISRDTNKSQVFLKMNSLQPEDTAMYFC
    """

    def __init__(
        self,
        name="heavy_chain",
        fwr1_aa="",
        cdr1_aa="",
        fwr2_aa="",
        cdr2_aa="",
        fwr3_aa="",
        cdr3_aa="",
        fwr4_aa="",
        v_gene="",
        j_gene="",
        species="",
        locus="IGH",
        leader="",
        tail="",
        region_def="imgt",
    ):
        """Constructor for HeavyChainAA object

        Parameters
        ----------
        name : str, optional
           name of object, by default "heavy_chain"
        fwr1_aa : str
           framework 1 amino acid sequence
        cdr1_aa : str
            cdr1 amino acid sequence
        fwr2_aa : str
            framework 2 amino acid sequence
        cdr2_aa : str
            cdr 2 amino acid sequence
        fwr3_aa : str
            framework 3 amino acid sequence
        cdr3_aa : str
            cdr3 amino acid sequence
        fwr4_aa : str
            framework 4 amino acid sequence
        v_gene : str
           v gene name, ex. "IGHV3-15"
        j_gene : str
           J gene name, ex. "IGHJ6*01"
        species : str
           species name, ex "human"
        locus : str
            locus name, defaults "IGH"
        leader : str
            leading amino acid sequence, defaults ""
        tail : str
            tailing amino acid sequence, defaults ""
        """
        super().__init__(
            name,
            fwr1_aa,
            cdr1_aa,
            fwr2_aa,
            cdr2_aa,
            fwr3_aa,
            cdr3_aa,
            fwr4_aa,
            v_gene,
            j_gene,
            species,
            leader,
            tail,
            region_def,
        )
        self._locus = locus

    @property
    def locus(self) -> str:
        """get locus name

        Returns
        -------
        str
            locus
        """
        return self._locus

    def get_json(self, indent=4):
        """Return json serialized string object

        This is useful for serializing the object to pass around

        Parameters
        ----------
        indent : int, optional
            how much to indent the json string, by default 4

        Returns
        -------
        json
            json encoded object

        Example
        -------

        >>> print(heavy_chain_aa.get_json())
        {
            "name": "heavy_chain",
            "fwr1_aa": "QVQLKESGPGLVQPSQTLSLTCTVS",
            "fwr2_aa": "VSWIRQPPGKGLEWMGV",
            "fwr3_aa": "DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
            "fwr4_aa": "WGPGTMVTVSS",
            "cdr1_aa": "GLSLTSNS",
            "cdr2_aa": "IWSNGGT",
            "cdr3_aa": "ASIYYYDADYLHWYFDF",
            "v_gene": "IGHV2-47",
            "species": "rat",
            "j_gene": "IGHJ1"
            "locus:"IGH"
            "leader":"",
            "tail":""
        }
        """
        dictionary = json.loads(super().get_json())
        dictionary["locus"] = "IGH"
        return json.dumps(dictionary, indent=indent)

    @staticmethod
    def read_json(file) -> "HeavyChainAA":
        """
        Read json file into HeavyChainAA object

        Parameters
        ----------
        file : str path
            file path string

        Returns
        -------
        HeavyChainAA
           AntibodyChainAA object
        """
        return HeavyChainAA(**json.load(open(file)))

    @staticmethod
    def from_json(json_object) -> "HeavyChainAA":
        """
        Read json string and return HeavyChainAA object

        Returns
        -------
        HeavyChainAA
           HeavyChainAA Object
        """
        return HeavyChainAA(**json.loads(json_object))


class KappaChainAA(AntibodyChainAA):
    """kappa chain object of amino acid sequences

    Examples
    --------
    >>> kappa_chain_aa = antibody.KappaChainAA(
    ... fwr1_aa="DIQMTQSPASLSASLGETVSIECLAS",
    ... cdr1_aa="EGISNS",
    ... fwr2_aa="LAWYQLKPGKSPQFLI",
    ... cdr2_aa="YATS",
    ... fwr3_aa="SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVYYC",
    ... cdr3_aa="QQGYKFPLT",
    ... fwr4_aa="FGSGTKLKIK",
    ... v_gene="IGKV12S11",
    ... j_gene="IGKJ5*01",
    ... species="rat",
    )
    Show segmented VDJ recombination of amino acids

    >>> kappa_chain_aa.get_segmented_vdj_aa()
    DIQMTQSPASLSASLGETVSIECLAS EGISNS LAWYQLKPGKSPQFLI YATS SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVYYC QQGYKFPLT FGSGTKLKIK

    Get segmented alignments between gemrline and mature

    >>> kappa_chain_aa.get_segmented_alignment_aa()
    kappa_chain DIQMTQSPHSLSASLGETVSIECLAS EGISNY LAWYQQKPGKSPQLLIY YA-S SLQDGVPSRFSGSGSGTQYSLKISNMQPEDEGVYYC QQGYKYPLT FGSGTKLEIK
    IGKV12S11   ........A................. .....S .....L.......F..- ..T. ........................G........... .....F... .......K.. IGKJ5*01

    Get mature or germline of specifc segments

    >>> kappa_chain_aa.fwr1_aa
    FrameWork1AA 1-27: DIQMTQSPASLSASLGETVSIECLAS

    >>> kappa_chain_aa.fwr1_aa_germline
    DIQMTQSPHSLSASLGETVSIECLAS

    >>> kappa_chain_aa.fwr2_aa
    FrameWork2AA 33-49: LAWYQLKPGKSPQFLI

    >>> kappa_chain_aa.fwr2_aa_germline
    LAWYQQKPGKSPQLLIY
    """

    def __init__(
        self,
        name="kappa_chain",
        fwr1_aa="",
        cdr1_aa="",
        fwr2_aa="",
        cdr2_aa="",
        fwr3_aa="",
        cdr3_aa="",
        fwr4_aa="",
        v_gene="",
        j_gene="",
        species="",
        locus="IGK",
        leader="",
        tail="",
        region_def="imgt",
    ):
        """Constructor for KappaChainAA object

        Parameters
        ----------
        name : str, optional
           name of object, by default "KappaChain"
        fwr1_aa : str
           framework 1 amino acid sequence
        cdr1_aa : str
            cdr1 amino acid sequence
        fwr2_aa : str
            framework 2 amino acid sequence
        cdr2_aa : str
            cdr 2 amino acid sequence
        fwr3_aa : str
            framework 3 amino acid sequence
        cdr3_aa : str
            cdr3 amino acid sequence
        fwr4_aa : str
            framework 4 amino acid sequence
        v_gene : str
           v gene name, ex. "IGKV3-11"
        j_gene : str
           J gene name, ex. "IGKJ1*01"
        species : str
           species name, ex "human"
        locus : str
            locus name, defaults "IGK"
        leader : str
            leading amino acid sequence, defaults ""
        tail : str
            tailing amino acid sequence, defaults ""
        """
        super().__init__(
            name,
            fwr1_aa,
            cdr1_aa,
            fwr2_aa,
            cdr2_aa,
            fwr3_aa,
            cdr3_aa,
            fwr4_aa,
            v_gene,
            j_gene,
            species,
            leader,
            tail,
            region_def,
        )
        self.locus = locus

    def get_json(self, indent=4) -> json:
        """Return json serialized string object

        This is useful for serializing the object to pass around

        Parameters
        ----------
        indent : int, optional
            how much to indent the json string, by default 4

        Returns
        -------
        json
            json encoded object

        Example
        -------

        >>> print(kappa_chain_aa.get_json())
        {
            "name": "kappa_chain",
            "fwr1_aa": "DIQMTQSPASLSASLGETVSIECLAS",
            "fwr2_aa": "LAWYQLKPGKSPQFLI",
            "fwr3_aa": "SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVYYC",
            "fwr4_aa": "FGSGTKLKIK",
            "cdr1_aa": "EGISNS",
            "cdr2_aa": "YATS",
            "cdr3_aa": "QQGYKFPLT",
            "v_gene": "IGKV12S11",
            "species": "rat",
            "j_gene": "IGKJ5*01",
            "locus": "IGK"
            "leader":"",
            "tail":""
        }
        """
        dictionary = json.loads(super().get_json())
        dictionary["locus"] = "IGK"
        return json.dumps(dictionary, indent=indent)

    @staticmethod
    def read_json(file: str) -> "KappaChainAA":
        """
        Read json file into KappaChainAA object

        Parameters
        ----------
        file : str path
            file path string

        Returns
        -------
        KappaChainAA
           KappaChainAA object
        """
        return KappaChainAA(**json.load(open(file)))

    @staticmethod
    def from_json(json_object: json) -> "KappaChainAA":
        """
        Read json string and return KappaChainAA object

        Returns
        -------
        KappaChainAA
          KappaChainAA
        """
        return KappaChainAA(**json.loads(json_object))


class LambdaChainAA(AntibodyChainAA):
    """lambda chain object of amino acid sequences

    Examples
    --------

    >>> lambda_chain_aa = antibody.LambdaChainAA(
    ... name="fezakinumab",
    ... fwr1_aa="QAVLTQPPSVSGAPGQRVTISCTGS",
    ... cdr1_aa="SSNIGAGYG",
    ... fwr2_aa="VHWYQQLPGTAPKLLIY",
    ... cdr2_aa="GDS",
    ... fwr3_aa="NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC",
    ... cdr3_aa="QSYDNSLSGYV",
    ... fwr4_aa="FGGGTQLTVL",
    ... v_gene="IGLV1-40*01",
    ... j_gene="IGLJ7*01",
    ... species="human",
    )
    Show segmented VDJ recombination of amino acids

    >>> lambda_chain_aa.get_segmented_vdj_aa()
    'QAVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYG VHWYQQLPGTAPKLLIY GDS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDNSLSGYV FGGGTQLTVL'

    Get segmented alignments between gemrline and mature

    >>> lambda_chain_aa.get_segmented_alignment_aa()
    fezakinumab QSVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYD VHWYQQLPGTAPKLLIY GNS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDSSLSGAV FGGGTQLTVL
    IGLV1-40*01 .A....................... ........G ................. .D. .................................... ....N....Y. .......... IGLJ7*01

    Get mature or germline of specifc segments

    >>> lambda_chain_aa.fwr1_aa
    FrameWork1AA 1-26: QAVLTQPPSVSGAPGQRVTISCTGSk

    >>> lambda_chain_aa.fwr1_aa_germline
    QSVLTQPPSVSGAPGQRVTISCTGS

    >>> lambda_chain_aa.fwr2_aa
    FrameWork2AA 35-52: VHWYQQLPGTAPKLLIY

    >>> lamba_chain_aa.fwr2_aa_germline
    VHWYQQLPGTAPKLLIY
    """

    def __init__(
        self,
        name="lambda_chain",
        fwr1_aa="",
        cdr1_aa="",
        fwr2_aa="",
        cdr2_aa="",
        fwr3_aa="",
        cdr3_aa="",
        fwr4_aa="",
        v_gene="",
        j_gene="",
        species="",
        locus="IGL",
        leader="",
        tail="",
        region_def="imgt",
    ):
        """Constructor for LambdaChainAA object

        Parameters
        ----------
        name : str, optional
           name of object, by default "lambda_chain"
        fwr1_aa : str
           framework 1 amino acid sequence
        cdr1_aa : str
            cdr1 amino acid sequence
        fwr2_aa : str
            framework 2 amino acid sequence
        cdr2_aa : str
            cdr 2 amino acid sequence
        fwr3_aa : str
            framework 3 amino acid sequence
        cdr3_aa : str
            cdr3 amino acid sequence
        fwr4_aa : str
            framework 4 amino acid sequence
        v_gene : str
           v gene name, ex. "IGLV10-67"
        j_gene : str
           J gene name, ex. "IGLJ9*01"
        species : str
           species name, ex "human"
        locus : str
            locus name, defaults "IGL"
        leader : str
            leading amino acid sequence, defaults ""
        tail : str
            tailing amino acid sequence, defaults ""
        """
        super().__init__(
            name,
            fwr1_aa,
            cdr1_aa,
            fwr2_aa,
            cdr2_aa,
            fwr3_aa,
            cdr3_aa,
            fwr4_aa,
            v_gene,
            j_gene,
            species,
            leader,
            tail,
            region_def,
        )
        self.locus = locus

    def get_json(self, indent=4) -> json:
        """Return json serialized string object

        This is useful for serializing the object to pass around

        Parameters
        ----------
        indent : int, optional
            how much to indent the json string, by default 4

        Returns
        -------
        json
            json encoded object

        Example
        -------

        >>> print(lambda_chain_aa.get_json())
        {
            "name": "fezakinumab",
            "fwr1_aa": "QAVLTQPPSVSGAPGQRVTISCTGS",
            "fwr2_aa": "VHWYQQLPGTAPKLLIY",
            "fwr3_aa": "NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC",
            "fwr4_aa": "FGGGTQLTVL",
            "cdr1_aa": "SSNIGAGYG",
            "cdr2_aa": "GDS",
            "cdr3_aa": "QSYDNSLSGYV",
            "v_gene": "IGLV1-40*01",
            "species": "human",
            "j_gene": "IGLJ7*01",
            "locus": "IGL"
            "leader":"",
            "tail":""
        }
        """
        dictionary = json.loads(super().get_json())
        dictionary["locus"] = "IGL"
        return json.dumps(dictionary, indent=indent)

    @staticmethod
    def read_json(file: str) -> "LambdaChainAA":
        """
        Read json file into LambdaChainAA object

        Parameters
        ----------
        file : str path
            file path string

        Returns
        -------
        LambdaChainAA
           LambdaChainAA object
        """
        return LambdaChainAA(**json.load(open(file)))

    @staticmethod
    def from_json(json_object: json) -> "LambdaChainAA":
        """
        Read json string and return LambdaChainAA object

        Returns
        -------
        LambdaChainAA
          LambdaChainAA
        """
        return LambdaChainAA(**json.loads(json_object))


class HeavyChainNT(AntibodyChainNT):
    """heavy chain object of nucleotide sequences


    Examples
    --------
    >>> chain_nt = antibody.HeavyChainNT(
    ... fwr1_nt="CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA",
    ... cdr1_nt="GGCCTATCGCTCACCAGCAACTCC",
    ... fwr2_nt="GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG",
    ... cdr2_nt="ATTTGGTCCAACGGTGGCACC",
    ... fwr3_nt="GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC",
    ... cdr3_nt="GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC",
    ... fwr4_nt="TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC",
    ... v_gene="IGHV2-47",
    ... j_gene="IGHJ3",
    ... species="rat",
    ... locus="IGH")

    Show segmented VDJ recombination of nucleotides

    >>> heavy_chain_nt.get_segmented_vdj_nt()
    CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA GGCCTATCGCTCACCAGCAACTCC GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG ATTTGGTCCAACGGTGGCACC GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC

    Show segmented VDJ recombination of aa

    >>> heavy_chain_nt.get_segmented_vdj_aa()
    QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS


    Get segmented alignments between gemrline and mature nucleotides

    >>> heavy_chain_nt.get_segmented_alignment_nt()
    antibodychain CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT GGGTTATCATTAACCAGCAATAGT GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA ATATGGAGTAATGGAGGCACA GATTATAATTCAGCTATCAAATCCCGACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTTCTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT GCCAGAAA--------------------------ACAATTGGTTTGCTTAC TGGGGCCAAGGCACTCTGGTCACTGTCTCTTCAG
    IGHV2-47      ..G.....G..G......AGC..C......T....G...........A..T..T.....G..A.....C..G..A ..CC....GC.C........CTCC ..C........A..T.....G.....C..A.................T..G ..T...TCC..C..T.....C ..C..C..C..C......G.GAG...CT..TCT....A.C.C........T..ATC...G........G.........TCG..T...C....G..T..G..T........T..C ...TCC.TTTATTACTATGACGCTGACTACCTCC..TGG.AC..C.A..T. .......CC......A....G..C..GAGC..-C IGHJ3

    Get segmented alignments between gemrline and mature aa

    >>> heavy_chain_nt.get_segmented_alignment_aa()
    antibodychain QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC AR--------YYWYFDF WGPGTMVTVSS
    IGHV2-47      ......................... ........ ................. ....... ......E.....N................P........ .SIYYYDADYLH..... ........... IGHJ1"

    Get specifc segments for amino acids or nucletodies for mature or germline segments

    >>> heavy_chain_nt.fwr1_aa
    FrameWork1AA 1-26: QVQLKESGPGLVQPSQTLSLTCTVSe
    QVQLKESGPGLVQPSQTLSLTCTVS

    >>> heavy_chain_nt.fwr1_nt
    FrameWork1NT 1-76: CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA

    >>> heavy_chain_nt.fwr1_nt_germline
    CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT

    >>> heavy_chain_nt.fwr1_aa_germline
    QVQLKESGPGLVQPSQTLSLTCTVS
    """

    def __init__(
        self,
        name="heavy_chain",
        fwr1_nt="",
        cdr1_nt="",
        fwr2_nt="",
        cdr2_nt="",
        fwr3_nt="",
        cdr3_nt="",
        fwr4_nt="",
        v_gene="",
        j_gene="",
        species="",
        locus="IGH",
        leader="",
        tail="",
    ):
        """Constructor for HeavyChainNT object

        Parameters
        ----------
        name : str, optional
           name of object, by default "heavy_chain"
        fwr1_nt : str
           framework 1 nucleotide sequence
        cdr1_nt : str
            cdr1 nucleotide sequence
        fwr2_nt : str
            framework 2 nucleotide sequence
        cdr2_nt : str
            cdr 2 nucleotide sequence
        fwr3_nt : str
            framework 3 nucleotide sequence
        cdr3_nt : str
            cdr3 nucleotide sequence
        fwr4_nt : str
            framework 4 nucleotide sequence
        v_gene : str
           v gene name, ex. "IGHV3-15"
        j_gene : str
           J gene name, ex. "IGHJ6*01"
        species : str
           species name, ex "human"
        locus : str
            locus name, defaults "IGH"
        leader : str
            leading nt sequence before antibody
        tail : str
            tail nt sequence before antibody
        """
        super().__init__(
            name,
            fwr1_nt,
            cdr1_nt,
            fwr2_nt,
            cdr2_nt,
            fwr3_nt,
            cdr3_nt,
            fwr4_nt,
            v_gene,
            j_gene,
            species,
            leader,
            tail,
        )
        self.locus = locus

    def get_json(self, indent=4) -> json:
        """Return json serialized string object


        Parameters
        ----------
        indent : int, optional
            how much to indent the json string, by default 4

        Returns
        -------
        json
            json encoded object

        Example
        -------
        >>> print(heavy_chain_nt.get_json())
        {
            "name": "heavy_chain",
            "fwr1_nt": "CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA",
            "fwr2_nt": "GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG",
            "fwr3_nt": "GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC",
            "fwr4_nt": "TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC",
            "cdr1_nt": "GGCCTATCGCTCACCAGCAACTCC",
            "cdr2_nt": "ATTTGGTCCAACGGTGGCACC",
            "cdr3_nt": "GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC",
            "v_gene": "IGHV2-47",
            "species": "rat",
            "j_gene": "IGHJ3",
            "locus": "IGH"
            "leader":"",
            "tail":""
        }
        """
        dictionary = json.loads(super().get_json())
        dictionary["locus"] = "IGH"
        return json.dumps(dictionary, indent=indent)

    @staticmethod
    def read_json(file: str) -> "HeavyChainNT":
        """
        Read json file into HeavyChainNT object

        Parameters
        ----------
        file : str path
            file path string

        Returns
        -------
        HeavyChainNT
           HeavyChainNT object
        """

        return HeavyChainNT(**json.load(open(file)))

    @staticmethod
    def from_json(json_object: json) -> "HeavyChainNT":
        """
        Read json string and return HeavyChainNT object

        Returns
        -------
        HeavyChainNT
          HeavyChainNT
        """
        return HeavyChainNT(**json.loads(json_object))


class KappaChainNT(AntibodyChainNT):
    """Kappa chain object of nucleotide sequences


    Examples
    --------
    >>> kappa_chain_nt = antibody.KappaChainNT(
    ... fwr1_nt="GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC",
    ... cdr1_nt="CACTCAATCTACCGGAAT",
    ... fwr2_nt="CTGGCCTGGTACCAACAGAAACTCGGTGAGGCTCCAAAACTACTCATCTAC",
    ... cdr2_nt="AACGCCAAC",
    ... fwr3_nt="TCTCTGCAGACAGGAATCCCGTCTAGATTTAGCGGATCCGGCTCCGGTACCGACTTCACCCTGACCATTAGCTCCCTGCAGCCCGAGGATGTGGCGACCTATTTCTGC",
    ... cdr3_nt="CAACAGTACTATCGAGGATGGACG",
    ... fwr4_nt="TGGACGTTCGGTGGAGGTACAAAGCTGGAGCTG",
    ... v_gene="IGKV22S4",
    ... j_gene="IGKJ1",
    ... species="rat",
    )

    Show segmented VDJ recombination of nucleotides

    >>> kappa_chain_nt.get_segmented_vdj_nt()
    GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC CACTCAATCTACCGGAAT CTGGCCTGGTACCAACAGAAACTCGGTGAGGCTCCAAAACTACTCATCTAC AACGCCAAC TCTCTGCAGACAGGAATCCCGTCTAGATTTAGCGGATCCGGCTCCGGTACCGACTTCACCCTGACCATTAGCTCCCTGCAGCCCGAGGATGTGGCGACCTATTTCTGC CAACAGTACTATCGAGGATGGACG TGGACGTTCGGTGGAGGTACAAAGCTGGAGCTG

    Show segmented VDJ recombination of aa

    >>> kappa_chain_nt.get_segmented_vdj_aa()
    DIQMTHSPSLLSASVGDRVSLNCKAS HSIYRN LAWYQQKLGEAPKLLIY NAN SLQTGIPSRFSGSGSGTDFTLTISSLQPEDVATYFC QQYYRGWT WTFGGGTKLEL


    Get segmented alignments between gemrline and mature nucleotides

    >>> kappa_chain_nt.get_segmented_alignment_nt()
    kappa_chain GACATCCAGATGACCCAGTCTCCTTCATTCCTGTCTGCATCTGTGGGAGACAGAGTCACTATCAACTGCAAAGCAAGT CAGAATATTAACAGGTAC TTAAACTGGTACCAGCAAAAGCTTGGAGAAGCTCCCAAACTCCTGATATAT AATGCAAAC AGTTTGCAAACGGGCATCCCATCAAGGTTCAGTGGCAGTGGATCTGGTACTGATTTCACACTCACCATCAGCAGCCTGCAGCCTGAAGATGTTGCCACATATTTCTGC TTGCAGCATAATAGTTGGCCGGTGGACG TT--CGGTGGAGGCACCAAGCTGGAATTGAAAC
    IGKV22S4    ........A.....A..T..G........G...AG...G........T...C.C....G.C.G........G..CTCC ..CTCA..CT..C..A.T C.GGC.........A..G..A..C..T..G.....A.....A..C..C..C ..C..C... TC.C....G..A..A.....G..T..A..T..C..ATCC..C..C.....C..C.....C..G.....T...TC.........C..G.....G..G..C......... CAA...--..C..--.C.AG.A...... .GGA..T.C.GT.G.GGT.CAAA.CTGGAGCTG IGKJ1

    Get segmented alignments between gemrline and mature aa

    >>> kappa_chain_nt.get_segmented_alignment_aa()
    kappa_chain DIQMTQSPSFLSASVGDRVTINCKAS QNINRY LNWYQQKLGEAPKLLIY NAN SLQTGIPSRFSGSGSGTDFTLTISSLQPEDVATYFC LQHNSWPWT FGGGTKLELK-
    IGKV22S4    .....H...L.........SL..... HS.Y.N .A............... ... .................................... -.QYYRG.. WTF.GGTK.EL IGKJ1

    Get specifc segments for amino acids or nucletodies for mature or germline segments

    >>> kappa_chain_nt.fwr1_aa
    FrameWork1AA 1-27: DIQMTHSPSLLSASVGDRVSLNCKAS

    >>> kappa_chain_nt.fwr1_nt
    FrameWork1NT 1-79: GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC

    >>> kappa_chain_nt.fwr1_nt_germline
    GACATCCAGATGACCCAGTCTCCTTCATTCCTGTCTGCATCTGTGGGAGACAGAGTCACTATCAACTGCAAAGCAAGT

    >>> kappa_chain_nt.fwr1_aa_germline
    DIQMTQSPSFLSASVGDRVTINCKAS
    """

    def __init__(
        self,
        name="kappa_chain",
        fwr1_nt="",
        cdr1_nt="",
        fwr2_nt="",
        cdr2_nt="",
        fwr3_nt="",
        cdr3_nt="",
        fwr4_nt="",
        v_gene="",
        j_gene="",
        species="",
        locus="IGK",
        leader="",
        tail="",
    ):
        """Constructor for KappaChainNT object

        Parameters
        ----------
        name : str, optional
           name of object, by default "kappa_chain"
        fwr1_nt : str
           framework 1 nucleotide sequence
        cdr1_nt : str
            cdr1 nucleotide sequence
        fwr2_nt : str
            framework 2 nucleotide sequence
        cdr2_nt : str
            cdr 2 nucleotide sequence
        fwr3_nt : str
            framework 3 nucleotide sequence
        cdr3_nt : str
            cdr3 nucleotide sequence
        fwr4_nt : str
            framework 4 nucleotide sequence
        v_gene : str
           v gene name, ex. "IGKV22S4"
        j_gene : str
           J gene name, ex. "IGKJ1*01"
        species : str
           species name, ex "human"
        locus : str
            locus name, defaults "IGK"
        leader: str
            leadint nt sequeence before kappa,
        tail: str
            trailing nt sequeence before kappa,


        """
        super().__init__(
            name,
            fwr1_nt,
            cdr1_nt,
            fwr2_nt,
            cdr2_nt,
            fwr3_nt,
            cdr3_nt,
            fwr4_nt,
            v_gene,
            j_gene,
            species,
        )
        self.locus = locus

    def get_json(self, indent=4) -> json:
        """Return json serialized string object


        Parameters
        ----------
        indent : int, optional
            how much to indent the json string, by default 4

        Returns
        -------
        json
            json encoded object

        Example
        -------

        >>> print(kappa_chain_nt.get_json())
        {
        "name": "kappa_chain",
        "fwr1_nt": "GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC",
        "fwr2_nt": "CTGGCCTGGTACCAACAGAAACTCGGTGAGGCTCCAAAACTACTCATCTAC",
        "fwr3_nt": "TCTCTGCAGACAGGAATCCCGTCTAGATTTAGCGGATCCGGCTCCGGTACCGACTTCACCCTGACCATTAGCTCCCTGCAGCCCGAGGATGTGGCGACCTATTTCTGC",
        "fwr4_nt": "TGGACGTTCGGTGGAGGTACAAAGCTGGAGCTG",
        "cdr1_nt": "CACTCAATCTACCGGAAT",
        "cdr2_nt": "AACGCCAAC",
        "cdr3_nt": "CAACAGTACTATCGAGGATGGACG",
        "v_gene": "IGKV22S4",
        "species": "rat",
        "j_gene": "IGKJ1",
        "locus": "IGK"
        "leader": ""
        "tail": ""
        }
        """
        dictionary = json.loads(super().get_json())
        dictionary["locus"] = "IGK"
        return json.dumps(dictionary, indent=indent)

    @staticmethod
    def read_json(file: str) -> "KappaChainNT":
        """
        Read json file into KappaChainNT object

        Parameters
        ----------
        file : str path
            file path string

        Returns
        -------
        KappaChainNT
           KappaChainNT object
        """

        return KappaChainNT(**json.load(open(file)))

    @staticmethod
    def from_json(json_object: json) -> "KappaChainNT":
        """
        Read json string and return KappaChainNT object

        Returns
        -------
        KappaChainNT
          KappaChainNT
        """
        return KappaChainNT(**json.loads(json_object))


class LambdaChainNT(AntibodyChainNT):
    """Lambda chain object of nucleotide sequences


    Examples
    --------
    >>> lambda_chain_nt = antibody.LambdaChainNT(
    ... name="fezakinumab",
    ... fwr1_nt="CAGGCGGTGCTCACCCAGCCACCTAGTGTGAGCGGTGCACCTGGGCAGCGTGTGACCATCTCTTGCACTGGGTCC",
    ... cdr1_nt="TCTTCCAACATCGGCGCCGGTTACGGC",
    ... fwr2_nt="GTGCACTGGTACCAACAGCTTCCGGGCACCGCCCCCAAGCTGCTCATCTAC",
    ... cdr2_nt="GGCGACAGC",
    ... fwr3_nt="AATCGTCCATCAGGGGTTCCGGATCGCTTTAGCGGGTCTAAGTCAGGGACCTCAGCCTCCCTGGCGATCACTGGGCTGCAGGCGGAGGACGAGGCAGACTATTACTGC",
    ... cdr3_nt="CAGTCTTATGACAATTCCTTGAGTGGC",
    ... fwr4_nt="TTCGGGGGAGGGACCCAGTTGACTGTTCTT",
    ... v_gene="IGLV1-40*01",
    ... j_gene="IGLJ7*01",
    ... species="human",
    )

    Show segmented VDJ recombination of nucleotides

    >>> lambda_chain_nt.get_segmented_vdj_nt()
    CAGGCGGTGCTCACCCAGCCACCTAGTGTGAGCGGTGCACCTGGGCAGCGTGTGACCATCTCTTGCACTGGGTCC TCTTCCAACATCGGCGCCGGTTACGGC GTGCACTGGTACCAACAGCTTCCGGGCACCGCCCCCAAGCTGCTCATCTAC GGCGACAGC AATCGTCCATCAGGGGTTCCGGATCGCTTTAGCGGGTCTAAGTCAGGGACCTCAGCCTCCCTGGCGATCACTGGGCTGCAGGCGGAGGACGAGGCAGACTATTACTGC CAGTCTTATGACAATTCCTTGAGTGGC TTCGGGGGAGGGACCCAGTTGACTGTTCTT

    Show segmented VDJ recombination of aa

    >>> lambda_chain_nt.get_segmented_vdj_aa()
    QAVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYG VHWYQQLPGTAPKLLIY GDS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDNSLSG FGGGTQLTVL


    Get segmented alignments between gemrline and mature nucleotides

    >>> lambda_chain_nt.get_segmented_alignment_nt()
    fezakinumab CAGTCTGTGCTGACGCAGCCGCCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACCATCTCCTGCACTGGGAGC AGCTCCAACATCGGGGCAGGTTATGAT GTACACTGGTACCAGCAGCTTCCAGGAACAGCCCCCAAACTCCTCATCTAT GGTAACAGC AATCGGCCCTCAGGGGTCCCTGACCGATTCTCTGGCTCCAAGTCTGGCACCTCAGCCTCCCTGGCCATCACTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACTGC CAGTCCTATGACAGCAGCCTGAGTGGTTCTGCTGTG TTCGGAGGAGGCACCCAGCTGACCGTCCTCG
    IGLV1-40*01 ...G.G.....C..C.....A..TAGT...AGC..T..A..T......C.T..G........T.........TC. TCT...........C..C.....C.GC ..G...........A........G..C..C........G..G........C ..CG..... .....T..A........T..G..T..C..TAGC..G..T.....A..G.................G...........G.....G.....C.....A..C......... .....T.......ATTC.T.....-.----..---- .....G.....G......T....T..T..-T IGLJ7*01

    Get segmented alignments between gemrline and mature aa

    >>> lambad_chain_nt.get_segmented_alignment_aa()
    fezakinumab QSVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYD VHWYQQLPGTAPKLLIY GNS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDSSLSGAV FGGGTQLTVL
    IGLV1-40*01 .A....................... ........G ................. .D. .................................... ....N....-- .......... IGLJ7*01'0


    Get specifc segments for amino acids or nucletodies for mature or germline segments

    >>> lambda_chain_nt.fwr1_aa
    FrameWork1AA 1-26: QAVLTQPPSVSGAPGQRVTISCTGS

    >>> lambda_chain_nt.fwr1_nt
    FrameWork1NT 1-76: CAGGCGGTGCTCACCCAGCCACCTAGTGTGAGCGGTGCACCTGGGCAGCGTGTGACCATCTCTTGCACTGGGTCC

    >>> kappa_chain_nt.fwr1_nt_germline
    CAGTCTGTGCTGACGCAGCCGCCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACCATCTCCTGCACTGGGAGC

    >>> kappa_chain_nt.fwr1_aa_germline
    QSVLTQPPSVSGAPGQRVTISCTGS
    """

    def __init__(
        self,
        name="lambda_chain",
        fwr1_nt="",
        cdr1_nt="",
        fwr2_nt="",
        cdr2_nt="",
        fwr3_nt="",
        cdr3_nt="",
        fwr4_nt="",
        v_gene="",
        j_gene="",
        species="",
        locus="IGL",
        leader="",
        tail="",
    ):
        """Constructor for LambdaChainNT object

        Parameters
        ----------
        name : str, optional
           name of object, by default "lambda_chain"
        fwr1_nt : str
           framework 1 nucleotide sequence
        cdr1_nt : str
            cdr1 nucleotide sequence
        fwr2_nt : str
            framework 2 nucleotide sequence
        cdr2_nt : str
            cdr 2 nucleotide sequence
        fwr3_nt : str
            framework 3 nucleotide sequence
        cdr3_nt : str
            cdr3 nucleotide sequence
        fwr4_nt : str
            framework 4 nucleotide sequence
        v_gene : str
           v gene name, ex. "IGLV1-40*01"
        j_gene : str
           J gene name, ex. "IGLJ7*01"
        species : str
           species name, ex "human"
        locus : str
            locus name, defaults "IGL"
        leader: str
            leadint nt sequeence before lambda,
        tail: str
            trailing nt sequeence before lambda,
        """
        super().__init__(
            name,
            fwr1_nt,
            cdr1_nt,
            fwr2_nt,
            cdr2_nt,
            fwr3_nt,
            cdr3_nt,
            fwr4_nt,
            v_gene,
            j_gene,
            species,
            leader,
            tail,
        )
        self.locus = locus

    def get_json(self, indent=4):
        """Return json serialized string object


        Parameters
        ----------
        indent : int, optional
            how much to indent the json string, by default 4

        Returns
        -------
        json
            json encoded object

        Example
        -------

        >>> print(lambda_chain_nt.get_json())
        {
            "name": "fezakinumab",
            "fwr1_nt": "CAGGCGGTGCTCACCCAGCCACCTAGTGTGAGCGGTGCACCTGGGCAGCGTGTGACCATCTCTTGCACTGGGTCC",
            "fwr2_nt": "GTGCACTGGTACCAACAGCTTCCGGGCACCGCCCCCAAGCTGCTCATCTAC",
            "fwr3_nt": "AATCGTCCATCAGGGGTTCCGGATCGCTTTAGCGGGTCTAAGTCAGGGACCTCAGCCTCCCTGGCGATCACTGGGCTGCAGGCGGAGGACGAGGCAGACTATTACTGC",
            "fwr4_nt": "TTCGGGGGAGGGACCCAGTTGACTGTTCTT",
            "cdr1_nt": "TCTTCCAACATCGGCGCCGGTTACGGC",
            "cdr2_nt": "GGCGACAGC",
            "cdr3_nt": "CAGTCTTATGACAATTCCTTGAGTGGC",
            "v_gene": "IGLV1-40*01",
            "species": "human",
            "j_gene": "IGLJ7*01",
            "locus": "IGL"
        }
        """
        dictionary = json.loads(super().get_json())
        dictionary["locus"] = "IGL"
        return json.dumps(dictionary, indent=indent)

    @staticmethod
    def read_json(file: str) -> "LambdaChainNT":
        """
        Read json file into LambdaChainNT object

        Parameters
        ----------
        file : str path
            file path string

        Returns
        -------
        LambdaChainNT
           LambdaChainNT object
        """

        return LambdaChainNT(**json.load(open(file)))

    @staticmethod
    def from_json(json_object):
        """
        Read json string and return LambdaChainNT object

        Returns
        -------
        LambdaChainNT
          LambdaChainNT
        """
        return LambdaChainNT(**json.loads(json_object))
