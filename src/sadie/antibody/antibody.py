import json
from typing import Union

from .chain import (
    HeavyChainAA,
    HeavyChainNT,
    KappaChainAA,
    KappaChainNT,
    LambdaChainAA,
    LambdaChainNT,
)
from .exception import HeavyChainException, LightChainException


class AntibodyAA:
    """Antibody Class with heavy and light chain pair of amino acid sequecnes

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
        ... species="human")
    >>> heavy_chain_aa = antibody.HeavyChainAA(
        ... name="2165H",
        ... fwr1_aa="QVQLKESGPGLVQPSQTLSLTCTVS",
        ... cdr1_aa="GLSLTSNS",
        ... fwr2_aa="VSWIRQPPGKGLEWMGV",
        ... cdr2_aa="IWSNGGT",
        ... fwr3_aa="DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
        ... cdr3_aa="ASIYYYDADYLHWYFDF",
        ... fwr4_aa="WGPGTMVTVSS",
        ... v_gene="IGHV2-47",
        ... j_gene="IGHJ3",
        ... species="rat")

    >>> antibody_aa = antibody.AntibodyAA(heavy_chain_aa, lambda_chain_aa)
    >>> antibody_aa
    Heavy
    2165H
    FrameWork1AA 1-25:QVQLKESGPGLVQPSQTLSLTCTVS
    CDR1AA 26-33:GLSLTSNS
    FrameWork2AA 34-50:VSWIRQPPGKGLEWMGV
    CDR2AA 51-57:IWSNGGT
    FrameWork3AA 58-95:DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC
    CDR3AA 96-112:ASIYYYDADYLHWYFDF
    FrameWork4AA 113-123:WGPGTMVTVSS

    Light
    fezakinumab
    FrameWork1AA 1-25:QAVLTQPPSVSGAPGQRVTISCTGS
    CDR1AA 26-34:SSNIGAGYG
    FrameWork2AA 35-51:VHWYQQLPGTAPKLLIY
    CDR2AA 52-54:GDS
    FrameWork3AA 55-90:NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC
    CDR3AA 91-101:QSYDNSLSGYV
    FrameWork4AA 102-111:FGGGTQLTVL
    """

    def __init__(self, heavy_chain: HeavyChainAA, light_chain: Union[KappaChainAA, LambdaChainAA]):
        """constructor object for heavy-light pair

        Parameters
        ----------
        heavy_chain : HeavyChainAA
            HeavyChainAA object
        light_chain : Union[KappaChainAA, LambdaChainAA]
            either a Kappa or Lambda chain object
        """
        self.heavy_chain = heavy_chain
        self.light_chain = light_chain

    @property
    def heavy_chain(self) -> HeavyChainAA:
        """heavy chain object getter

        Returns
        -------
        HeavyChainAA
            HeavyChainAA object
        """
        return self._heavy_chain

    @heavy_chain.setter
    def heavy_chain(self, hc: HeavyChainAA):
        """heavy chain setter

        Parameters
        ----------
        hc : HeavyChainAA
            heavy chain amino acid object

        Raises
        ------
        HeavyChainException
            if not a heavy chain AA
        """
        if not isinstance(hc, HeavyChainAA):
            raise HeavyChainException(f"{hc} not heavy chain")
        self._heavy_chain = hc

    @property
    def light_chain(self) -> Union[LambdaChainAA, KappaChainAA]:
        """light chain getter

        Returns
        -------
        Union[LambdaChainAA, KappaChainAA]
            returns a lambda chain or kappa chain amino acid object
        """
        return self._light_chain

    @light_chain.setter
    def light_chain(self, lc: Union[LambdaChainAA, KappaChainAA]):
        """light chain setter

        Parameters
        ----------
        lc : Union[LambdaChainAA, KappaChainAA]
            lambda chain or kappa chain amino acid object

        Raises
        ------
        LightChainException
            if not a lamba or kappa chain object
        """
        if not isinstance(lc, (KappaChainAA, LambdaChainAA)):
            raise LightChainException(f"{lc} not lambda or Kappa chain")
        self._light_chain = lc

    def get_segmented_vdj_aa(self) -> str:
        """Get segmented VDJ amino acid object

        Returns
        -------
        str
            vdj string segmentet

        Examples
        --------

        >>> antibody_aa.get_segmented_vdj_aa()
        QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS

        QAVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYG VHWYQQLPGTAPKLLIY GDS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDNSLSGYV FGGGTQLTVL
        """
        return "{}\n\n{}".format(
            self.heavy_chain.get_segmented_vdj_aa(),
            self.light_chain.get_segmented_vdj_aa(),
        )

    def get_segmented_alignment_aa(self, line_length=80) -> str:
        """Get segmented alignment of HL pair

        Parameters
        ----------
        line_length : int, optional
            the line length to break, by default 80

        Returns
        -------
        str
            segmented formatted string

        Examples
        --------

        >>> object.get_segmented_alignment_aa()

        2165H           QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKS
        IGHV2-47|IGHJ3  ......................... ........ ................. ....... ......E.....N......

        2165H           QVFLKMNSLQTEDTAMYFC AR----------NWFAY WGQGTLVTVSS
        IGHV2-47|IGHJ3  ..........P........ .SIYYYDADYLHWY.DF ..P..M.....



        fezakinumab           QSVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYD VHWYQQLPGTAPKLLIY GNS NRPSGVPDRFSGSKSGTSASLA
        IGLV1-40*01|IGLJ7*01  .A....................... ........G ................. .D. ......................

        fezakinumab           ITGLQAEDEADYYC QSYDSSLSGAV FGGGTQLTVL
        IGLV1-40*01|IGLJ7*01  .............. ....N....Y. ..........
        """
        return "{}\n\n{}".format(
            self.heavy_chain.get_segmented_alignment_aa(line_length=line_length),
            self.light_chain.get_segmented_alignment_aa(line_length=line_length),
        )

    def get_json(self, indent=4) -> json:
        """return json string serialization of object

        Parameters
        ----------
        indent : int, optional
            indentation of json, by default 4

        Returns
        -------
        json
            AntibodyAA json represention

        Examples
        --------

        >>> object.get_json()
        {
        "heavy": {
            "name": "2165H",
            "fwr1_aa": "QVQLKESGPGLVQPSQTLSLTCTVS",
            "fwr2_aa": "VSWIRQPPGKGLEWMGV",
            "fwr3_aa": "DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
            "fwr4_aa": "WGPGTMVTVSS",
            "cdr1_aa": "GLSLTSNS",
            "cdr2_aa": "IWSNGGT",
            "cdr3_aa": "ASIYYYDADYLHWYFDF",
            "v_gene": "IGHV2-47",
            "species": "rat",
            "j_gene": "IGHJ3",
            "locus": "IGH"
        },
        "light": {
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
            }
        }
        """
        dictionary = {
            "heavy": json.loads(self.heavy_chain.get_json()),
            "light": json.loads(self.light_chain.get_json()),
        }
        return json.dumps(dictionary, indent=indent)

    @staticmethod
    def from_json(json_object: json) -> "AntibodyAA":
        """take in json object serialized and return AntibodyAA

        Returns
        -------
        AntibodyAA
            AntiboydAA object

        """
        dictionary = json.loads(json_object)
        heavy_ = HeavyChainAA(**dictionary["heavy"])
        if dictionary["light"]["locus"] == "IGK":
            light_ = KappaChainAA(**dictionary["light"])
        elif dictionary["light"]["locus"] == "IGL":
            light_ = LambdaChainAA(**dictionary["light"])
        else:
            raise KeyError("IGK or IGL not in json file")
        return AntibodyAA(heavy_, light_)

    @staticmethod
    def read_json(file: str) -> "AntibodyAA":
        """Read a json file and deserialize to AntibodyAAObject

        Returns
        -------
        AntibodyAA
            AntibodyAA Object
        """
        dictionary = json.load(open(file))
        heavy_ = HeavyChainAA(**dictionary["heavy"])
        if dictionary["light"]["locus"] == "IGK":
            light_ = KappaChainAA(**dictionary["light"])
        elif dictionary["light"]["locus"] == "IGL":
            light_ = LambdaChainAA(**dictionary["light"])
        else:
            raise KeyError("IGK or IGL not in json file")
        return AntibodyAA(heavy_, light_)

    def to_json(self, file: str):
        """serialize AntibodyAA to json file

        Parameters
        ----------
        file : str
            file path
        """
        data = json.loads(self.get_json())
        with open(file, "w") as outfile:
            json.dump(data, outfile)

    def __repr__(self):
        return f"Heavy\n{self.heavy_chain.__repr__()}\n\nLight\n{self.light_chain.__repr__()}"

    def __eq__(self, other):
        return all(
            [
                self.heavy_chain == other.heavy_chain,
                self.light_chain == other.light_chain,
            ]
        )


class AntibodyNT:
    """Antibody Class with heavy and light chain pair of nucleotide sequecnes

    Examples
    --------
    >>> kappa_chain_nt = antibody.KappaChainNT(
        ... name='2165Kappa',
        ... fwr1_nt="GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC",
        ... cdr1_nt="CACTCAATCTACCGGAAT",
        ... fwr2_nt="CTGGCCTGGTACCAACAGAAACTCGGTGAGGCTCCAAAACTACTCATCTAC",
        ... cdr2_nt="AACGCCAAC",
        ... fwr3_nt="TCTCTGCAGACAGGAATCCCGTCTAGATTTAGCGGATCCGGCTCCGGTACCGACTTCACCCTGACCATTAGCTCCCTGCAGCCCGAGGATGTGGCGACCTATTTCTGC",
        ... cdr3_nt="CAACAGTACTATCGAGGATGGACG",
        ... fwr4_nt="TGGACGTTCGGTGGAGGTACAAAGCTGGAGCTG",
        ... v_gene="IGKV22S4",
        ... j_gene="IGKJ1",
        ... species="rat")

    >>> heavy_chain_nt = antibody.HeavyChainNT(
        ... name='2165Heavy',
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

    >>> antibody_nt = antibody.AntibodyNT(heavy_chain_nt,lambda_chain_nt)
    >>> antibody_nt
    Heavy
    2165Heavy
    FrameWork1NT 1-75:CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA
    CDR1NT 76-99:GGCCTATCGCTCACCAGCAACTCC
    FrameWork2NT 100-150:GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG
    CDR2NT 151-171:ATTTGGTCCAACGGTGGCACC
    FrameWork3NT 172-285:GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC
    CDR3NT 286-336:GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC
    FrameWork4NT 337-369:TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC

    Light
    fezakinumab
    FrameWork1NT 1-75:CAGGCGGTGCTCACCCAGCCACCTAGTGTGAGCGGTGCACCTGGGCAGCGTGTGACCATCTCTTGCACTGGGTCC
    CDR1NT 76-102:TCTTCCAACATCGGCGCCGGTTACGGC
    FrameWork2NT 103-153:GTGCACTGGTACCAACAGCTTCCGGGCACCGCCCCCAAGCTGCTCATCTAC
    CDR2NT 154-162:GGCGACAGC
    FrameWork3NT 163-270:AATCGTCCATCAGGGGTTCCGGATCGCTTTAGCGGGTCTAAGTCAGGGACCTCAGCCTCCCTGGCGATCACTGGGCTGCAGGCGGAGGACGAGGCAGACTATTACTGC
    CDR3NT 271-297:CAGTCTTATGACAATTCCTTGAGTGGC
    FrameWork4NT 298-327:TTCGGGGGAGGGACCCAGTTGACTGTTCTT

    """

    def __init__(self, heavy_chain: HeavyChainNT, light_chain: Union[KappaChainNT, LambdaChainNT]):
        """constructor object for heavy-light pair

        Parameters
        ----------
        heavy_chain : HeavyChainNT
            HeavyChainNT object
        light_chain : Union[KappaChainNT, LambdaChainNT]
            either a Kappa or Lambda chain object
        """
        self.heavy_chain = heavy_chain
        self.light_chain = light_chain

    @property
    def heavy_chain(self) -> HeavyChainNT:
        """heavy chain getter

        Returns
        -------
        HeavyChainNT
            HeavyChainNT object
        """
        return self._heavy_chain

    @heavy_chain.setter
    def heavy_chain(self, hc: HeavyChainNT):
        """heavy chain setter

        Parameters
        ----------
        hc : HeavyChainNT
            heavy chain nucleotide object

        Raises
        ------
        HeavyChainException
            if not a heavy chain NT
        """
        if not isinstance(hc, HeavyChainNT):
            raise HeavyChainException(f"{hc} not heavy chain")
        self._heavy_chain = hc

    @property
    def light_chain(self) -> Union[LambdaChainNT, KappaChainNT]:
        """get light chain

        Returns
        -------
        Union[LambdaChainNT,KappaChainNT]
            Either the lambda or kappa chain nucleotide object
        """
        return self._light_chain

    @light_chain.setter
    def light_chain(self, lc: Union[KappaChainNT, LambdaChainNT]):
        """set light chain

        Parameters
        ----------
        lc : Union[KappaChainNT, LambdaChainNT]
            Either the lambda or kappa chain nucleotide object

        Raises
        ------
        LightChainException
            if not a kappa or lambda chain NT
        """
        if not isinstance(lc, (KappaChainNT, LambdaChainNT)):
            raise LightChainException(f"{lc} not lambda or Kappa chain")
        self._light_chain = lc

    def get_segmented_vdj_nt(self):
        """Get segmented VDJ amino acid object

        Returns
        -------
        str
            vdj string segmentet

        Examples
        --------

        >>> antibody_aa.get_segmented_vdj_aa()
        QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS

        QAVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYG VHWYQQLPGTAPKLLIY GDS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDNSLSGYV FGGGTQLTVL
        """
        return "{}\n\n{}".format(
            self.heavy_chain.get_segmented_vdj_nt(),
            self.light_chain.get_segmented_vdj_nt(),
        )

    def get_segmented_vdj_aa(self):
        """Get segmented VDJ nucleotide object

        Returns
        -------
        str
            vdj string segment

        Examples
        --------

        >>> antibody_nt.get_segmented_vdj_nt()
        CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA GGCCTATCGCTCACCAGCAACTCC GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG ATTTGGTCCAACGGTGGCACC GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC

        CAGGCGGTGCTCACCCAGCCACCTAGTGTGAGCGGTGCACCTGGGCAGCGTGTGACCATCTCTTGCACTGGGTCC TCTTCCAACATCGGCGCCGGTTACGGC GTGCACTGGTACCAACAGCTTCCGGGCACCGCCCCCAAGCTGCTCATCTAC GGCGACAGC AATCGTCCATCAGGGGTTCCGGATCGCTTTAGCGGGTCTAAGTCAGGGACCTCAGCCTCCCTGGCGATCACTGGGCTGCAGGCGGAGGACGAGGCAGACTATTACTGC CAGTCTTATGACAATTCCTTGAGTGGC TTCGGGGGAGGGACCCAGTTGACTGTTCTT
        """
        return "{}\n\n{}".format(
            self.heavy_chain.get_segmented_vdj_aa(),
            self.light_chain.get_segmented_vdj_aa(),
        )

    def get_segmented_alignment_nt(self, line_length=80) -> str:
        """Get segmented alignment of HL pair

        Parameters
        ----------
        line_length : int, optional
            the line length to break, by default 80

        Returns
        -------
        str
            segmented formatted string

        Examples
        --------

        >>> object.get_segmented_alignment_nt()
        2165Heavy       CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT GGGT
        IGHV2-47|IGHJ3  ..G.....G..G......AGC..C......T....G...........A..T..T.....G..A.....C..G..A ..CC

        2165Heavy       TATCATTAACCAGCAATAGT GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA ATATGGA
        IGHV2-47|IGHJ3  ....GC.C........CTCC ..C........A..T.....G.....C..A.................T..G ..T...T

        2165Heavy       GTAATGGAGGCACA GATTATAATTCAGCTATCAAATCCCGACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTT
        IGHV2-47|IGHJ3  CC..C..T.....C ..C..C..C..C......G.GAG...CT..TCT....A.C.C........T..ATC...G.....

        2165Heavy       CTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT GCCAGAAA----------------------
        IGHV2-47|IGHJ3  ...G.........TCG..T...C....G..T..G..T........T..C ...TCC.TTTATTACTATGACGCTGACTAC

        2165Heavy       ----ACAATTGGTTTGCTTAC TGGGGCCAAGGCACTCTGGTCACTGTCTCTTCAG
        IGHV2-47|IGHJ3  CTCC..TGG.AC..C.A..T. .......CC......A....G..C..GAGC..-C



        fezakinumab           CAGTCTGTGCTGACGCAGCCGCCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACCATCTCCTGCACTGGGAGC AGCT
        IGLV1-40*01|IGLJ7*01  ...G.G.....C..C.....A..TAGT...AGC..T..A..T......C.T..G........T.........TC. TCT.

        fezakinumab           CCAACATCGGGGCAGGTTATGAT GTACACTGGTACCAGCAGCTTCCAGGAACAGCCCCCAAACTCCTCATCTAT GGTA
        IGLV1-40*01|IGLJ7*01  ..........C..C.....C.GC ..G...........A........G..C..C........G..G........C ..CG

        fezakinumab           ACAGC AATCGGCCCTCAGGGGTCCCTGACCGATTCTCTGGCTCCAAGTCTGGCACCTCAGCCTCCCTGGCCATCACTGG
        IGLV1-40*01|IGLJ7*01  ..... .....T..A........T..G..T..C..TAGC..G..T.....A..G.................G........

        fezakinumab           GCTCCAGGCTGAGGATGAGGCTGATTATTACTGC CAGTCCTATGACAGCAGCCTGAGTGGTTCTGCTGTG TTCGGAGG
        IGLV1-40*01|IGLJ7*01  ...G.....G.....C.....A..C......... .....T.......ATTC.T.......---------C .....G..

        fezakinumab           AGGCACCCAGCTGACCGTCCTCG
        IGLV1-40*01|IGLJ7*01  ...G......T....T..T..-T
        """
        return "{}\n\n{}".format(
            self.heavy_chain.get_segmented_alignment_nt(line_length=line_length),
            self.light_chain.get_segmented_alignment_nt(line_length=line_length),
        )

    def get_segmented_alignment_aa(self, line_length=80) -> str:
        """Get segmented alignment of HL pair

        Parameters
        ----------
        line_length : int, optional
            the line length to break, by default 80

        Returns
        -------
        str
            segmented formatted string

        Examples
        --------

        >>> object.get_segmented_alignment_aa()

        2165H           QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKS
        IGHV2-47|IGHJ3  ......................... ........ ................. ....... ......E.....N......

        2165H           QVFLKMNSLQTEDTAMYFC AR----------NWFAY WGQGTLVTVSS
        IGHV2-47|IGHJ3  ..........P........ .SIYYYDADYLHWY.DF ..P..M.....



        fezakinumab           QSVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYD VHWYQQLPGTAPKLLIY GNS NRPSGVPDRFSGSKSGTSASLA
        IGLV1-40*01|IGLJ7*01  .A....................... ........G ................. .D. ......................

        fezakinumab           ITGLQAEDEADYYC QSYDSSLSGAV FGGGTQLTVL
        IGLV1-40*01|IGLJ7*01  .............. ....N....Y. ..........
        """
        return "{}\n\n{}".format(
            self.heavy_chain.get_segmented_alignment_aa(line_length=line_length),
            self.light_chain.get_segmented_alignment_aa(line_length=line_length),
        )

    def get_json(self, indent=4) -> json:
        """return json string serialization of object

        Parameters
        ----------
        indent : int, optional
            indentation of json, by default 4

        Returns
        -------
        json
            AntibodyNT json represention

        Examples
        --------

        >>> object.get_json()
        {
            "heavy": {
                "name": "2165Heavy",
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
            },
            "light": {
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
        }
        """
        dictionary = {
            "heavy": json.loads(self.heavy_chain.get_json()),
            "light": json.loads(self.light_chain.get_json()),
        }
        return json.dumps(dictionary, indent=indent)

    @staticmethod
    def from_json(json_object: json) -> "AntibodyNT":
        """take in json object serialized and return AntibodyNT

        Returns
        -------
        AntibodyNT
            AntiboydNT object

        """
        dictionary = json.loads(json_object)
        heavy_ = HeavyChainNT(**dictionary["heavy"])
        if dictionary["light"]["locus"] == "IGK":
            light_ = KappaChainNT(**dictionary["light"])
        elif dictionary["light"]["locus"] == "IGL":
            light_ = LambdaChainNT(**dictionary["light"])
        else:
            raise KeyError("IGK or IGL not in json file")
        return AntibodyNT(heavy_, light_)

    @staticmethod
    def read_json(file: str) -> "AntibodyNT":
        """Read a json file and deserialize to AntibodyNT object

        Returns
        -------
        AntibodyNT
            AntibodyNT Object
        """
        dictionary = json.load(open(file))
        heavy_ = HeavyChainNT(**dictionary["heavy"])
        if dictionary["light"]["locus"] == "IGK":
            light_ = KappaChainNT(**dictionary["light"])
        elif dictionary["light"]["locus"] == "IGL":
            light_ = LambdaChainNT(**dictionary["light"])
        else:
            raise KeyError("IGK or IGL not in json file")
        return AntibodyNT(heavy_, light_)

    def to_json(self, file: str):
        """serialize AntibodyNT to json file

        Parameters
        ----------
        file : str
            file path
        """
        data = json.loads(self.get_json())
        with open(file, "w") as outfile:
            json.dump(data, outfile)

    def __repr__(self):
        return f"Heavy\n{self.heavy_chain.__repr__()}\n\nLight\n{self.light_chain.__repr__()}"

    def __eq__(self, other):
        return all(
            [
                self.heavy_chain == other.heavy_chain,
                self.light_chain == other.light_chain,
            ]
        )
