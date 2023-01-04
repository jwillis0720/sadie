"""A collection of base classes for higher level antibody objects"""
import re
import warnings

from Bio.Align import MultipleSeqAlignment
from Bio.pairwise2 import align

# Third Party
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy import isnan

from sadie.utility.util import format_alignment

# Lib Level
from .exception import (
    BadAASequenceError,
    BadAASequenceWarning,
    BadNTSequenceError,
    BadNTSequenceWarning,
)

# Common amino acids
STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVYW"
EXTRA_AMINO_ACIDS = "X*"
STANDARD_NT = "ATCG"
EXTRA_NT = "NRYX"


class AntibodySegment:
    """
    Base class for defining antibody segments, either amino acids or nucletoide sequecnes
    """

    def __init__(self, sequence):
        """
        constructor for generic antibody segments


        Parameters
        ----------
        sequence : str
            amino acid or nucleotide sequence for the antibody segment


        Examples
        --------
        >>> segment = AntibodySegment("DIQMTQSPASLSASLGETVSIECLAS")
        """
        self.sequence = str(sequence)
        self._germline = ""
        self._alignment = ("", "", "", "")

    @property
    def sequence(self) -> str:
        return self._sequence

    @sequence.setter
    def sequence(self, seq: str):
        if seq:
            self._start_index = 0
            self._start = 1
            self._end_index = len(seq) - 1
            self._end = len(seq)
            self._sequence = seq
        elif isinstance(seq, str):
            self._sequence = ""
            self._start, self._end, self._start_index, self._end_index = None, None, None, None
        elif isnan(seq):
            self.sequence = ""
        else:
            raise TypeError(f"Setting sequence seq {seq} as {type(seq)} is currently not supported")

    @property
    def start_index(self) -> int:
        """return start index of antibody segment in 0 indexing


        Returns
        -------
        int
          start index
        """
        return self._start_index

    @start_index.setter
    def start_index(self, start):
        """set start index of segment in 0 indexing

        Parameters
        ----------
        start : int
            the start index
        """
        if not isinstance(start, int):
            raise ValueError(f"{start} must be integer")
        if start < 0:
            raise ValueError(f"{start} must be greater than 0")
        self._start_index = start
        self._start = start + 1
        self._end_index = start + len(self.sequence) - 1
        self._end = start + len(self.sequence)

    @property
    def start(self) -> int:
        """return the start of the antibody segment

        Returns
        -------
        int
            the start of the antibody segment
        """
        return self._start

    @start.setter
    def start(self, start):
        """set the start of the antibody segment

        Parameters
        ----------
        start : int
            the start of the antibody segment
        """
        # setting the start index will call start_index setter and set everything else
        self._start_index = start - 1

    @property
    def end_index(self) -> int:
        """return the 0-based index of the end

        Returns
        -------
        int
           0-based end index
        """
        return self._end_index

    @property
    def end(self) -> int:
        """return the end of the antibody segment

        Returns
        -------
        int
            the end  of the antibody segment
        """
        return self._end

    @property
    def germline(self) -> str:
        """get the germline sequence if set

        Returns
        -------
        str
            the germline sequence
        """
        return self._germline

    @germline.setter
    def germline(self, germ):
        """set the germline sequence.

        Will also set the alignment between the sequence and germline

        Parameters
        ----------
        germ : [str, Bio.Seq]
            the germline sequence

        Raises
        ------
        ValueError
            if the the germline sequence is not str or Bio.Seq
        """
        if not isinstance(germ, (type(self), str, Seq)):
            raise ValueError(f"{germ} must be instance of str or Bio.Seq")
        self._germline = str(germ)

        # If no seq set sequence as '-'
        if not self._sequence:
            _sequence = "-" * len(str(germ))
        else:
            _sequence = self._sequence

        # we can have it where germline segment is also empty
        if self._germline:
            _alignments = align.globalxs(self._germline, _sequence, -10, -1)
            self._alignment = _alignments[0]

    @property
    def alignment(self) -> tuple:
        """Get alignment tuple

        Example
        -------
        >>> object.alignment
        ('CAGGCGCAGCACG-------...---------C', 'CAGGTTCAGCTGGTGCAGTC...CAAGGCTTCT', -61.0, 0, 75)

        Returns
        -------
        tuple
            alignment tuple from biopython, eg (source, target, score, begin, end)
        """
        return self._alignment

    def get_formatted_alignment(
        self,
        source_name="germline",
        target_name="target",
        line_length=80,
        ljust=12,
        maxid=30,
    ) -> str:
        """Return a formatted alignment string using dots for matches

        Parameters
        ----------
        source_name : str, optional
            source name in alignment, by default "germline"
        target_name : str, optional
            target name in alignment string, by default "target"
        line_length : int, optional
            how long until break of line in alignment, by default 80
        ljust : int, optional
           buffer between alignment title and start of sequence, by default 12
        maxid : int, optional
            truncate id if longer than this line, by default 30

        Returns
        -------
        str
            formatted alignment str

        Example
        -------

        >>> print(object.get_formated_alignment())
        germline    CAGGCGCAGCACG-------------------------------------------------------------C
        target      ....TT....TG.TGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCT

        """
        if self._alignment == ("", "", "", ""):
            warnings.warn("No germline set", UserWarning)
            return self._alignment
        # grab germ and target
        _germ, _target = (self.alignment[0], self.alignment[1])
        formatted_germ, formatted_target = (_germ, "")

        # iterate through and replace same sequecne with dot
        for g, t in zip(_germ, _target):
            if g != t:
                formatted_target += t
            else:
                formatted_target += "."

        # Change to alignment object since its so much easier to deal
        alignment_object = MultipleSeqAlignment(
            [
                SeqRecord(Seq(formatted_germ), id=source_name),
                SeqRecord(Seq(formatted_target), id=target_name),
            ]
        )

        return format_alignment(alignment_object, line_length, ljust, maxid)

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return str(self.sequence) == other


class AntibodySegmentNT(AntibodySegment):
    """Antibody Segment base class of nucleotide sequences"""

    def __init__(self, nt):
        """Constructor for Antibody Segment

        Parameters
        ----------
        nt : str
            nucleotide sequence
        """
        super().__init__(nt)
        self.nt = nt

    @property
    def nt(self) -> str:
        """Return nuceltoide sequence

        Returns
        -------
        str
            nucleotide sequence
        """
        return self._nt

    @nt.setter
    def nt(self, nt_sequence):
        """set the nucletodie sequence

        Parameters
        ----------
        nt_sequence : str
            nucleotide sequecne

        Raises
        ------
        BadNTSequenceError
            If the nucleotide sequence contains bad characters not found in ATCGNRYX

        Warns
        -----
        BadNTSequenceWarning
            If an ambiguous nucleotide is passed NRYX
        """
        self._nt = self._validate_nt(nt_sequence)

    @property
    def aa(self) -> str:
        """return amino acid translation

        Returns
        -------
        str
            amino acid sequence
        """
        trim = len(self._nt) - (len(self._nt) % 3)
        return str(Seq(self._nt)[0:trim].translate())

    def _validate_nt(self, nt_seq):
        # If its false, let's set it to empty
        if not nt_seq:
            nt_seq = ""
        """Private Nucleotide validation method"""
        check_nt = re.match("^[{}]+$".format(STANDARD_NT.upper()), nt_seq.upper())
        if not bool(check_nt):
            accepted_nt = list(STANDARD_NT.upper()) + list(EXTRA_NT.upper())
            for position, nt in enumerate(nt_seq, start=1):
                if nt.upper() not in accepted_nt:
                    raise BadNTSequenceError(self.__class__.__name__, position, nt, "".join(accepted_nt))
                if nt.upper() in EXTRA_NT:
                    warnings.warn(
                        f"positon {position} in {self.__class__.__name__} is {nt}",
                        BadNTSequenceWarning,
                    )
        return nt_seq

    def __repr__(self):
        string = "{} {}-{}: {}".format(self.__class__.__name__, self.start, self.end, self.nt)
        return string

    def __add__(self, other):
        return self.nt + other.nt


class AntibodySegmentAA(AntibodySegment):
    """Antibody Segment base class of Amino acids"""

    def __init__(self, amino_acid):
        """
        constructor for antibody amino acid sequences


        Parameters
        ----------
        amino_acid : str
            amino acid sequence for the antibody segment
        """
        super().__init__(amino_acid)
        self.aa = amino_acid

    @property
    def aa(self) -> str:
        """returns amino acid sequence of antibody segmetn

        Returns
        -------
        str
            amino acid sequence for the antibody segment
        """
        return self._aa

    @aa.setter
    def aa(self, amino_acid):
        """set the amino acid sequence

        Parameters
        ----------
        amino_acid : str
            amino acid  sequecne

        Raises
        ------
        BadAASequenceError
           amino acid sequence contains bad characters not found in the ACDEFGHIKLMNPQRSTVYW*X

        Warns
        -----
        BadAASequenceWarning
            If an ambiguous amino acid X or stop codon * is passed
        """
        self._aa = self._validate_aa(amino_acid)

    def _validate_aa(self, amino_acid):
        """Private method for validating amino acid sequence"""
        check_aa = re.match("^[{}]+$".format(STANDARD_AMINO_ACIDS.upper()), amino_acid.upper())
        if not bool(check_aa):
            accepted_aa = list(STANDARD_AMINO_ACIDS.upper()) + list(EXTRA_AMINO_ACIDS.upper())
            for position, amino in enumerate(amino_acid, start=1):
                if amino.upper() not in accepted_aa:
                    raise BadAASequenceError(self.__class__.__name__, position, amino, "".join(accepted_aa))
                if amino.upper() in EXTRA_AMINO_ACIDS:
                    warnings.warn(
                        f"positon {position} in {self.__class__.__name__} is {amino}",
                        BadAASequenceWarning,
                    )
        return amino_acid

    def __repr__(self):
        string = "{} {}-{}: {}".format(self.__class__.__name__, self.start, self.end, self.aa)
        return string

    def __add__(self, other):
        return self.aa + other.aa


class FrameWork1NT(AntibodySegmentNT):
    """Framework 1 nucelotide segments

    Parameters
    ----------
    AntibodySegmentNT :
       AntibodySegment superclass
    """

    def __init__(self, seq):
        """constructur

        Parameters
        ----------
        seq : str
            nucleotide sequence

        Examples
        --------
        >>> fw1 = Framework1NT("GACGTTCAGCTGGTGGAAAGTGGGGGTGACCTCTTGAAGCCGGGGGGCTCTCTCAGGTTAACATGCGTGGCTCCC")
        """
        super().__init__(seq)


class FrameWork2NT(AntibodySegmentNT):
    """Framework 2 nucelotide segments

    Parameters
    ----------
    AntibodySegmentNT :
       AntibodySegment superclass
    """

    def __init__(self, seq):
        """constructur

        Parameters
        ----------
        seq : str
            nucleotide sequence

        Examples
        --------
        >>> fw2 = Framework2NT("ATGAACTGGGTCTGTCAGGCACCTGGTAAGGGACTCCAGTGGGTAGCCTAT")
        """
        super().__init__(seq)


class FrameWork3NT(AntibodySegmentNT):
    """Framework 3 nucelotide segments

    Parameters
    ----------
    AntibodySegmentNT :
       AntibodySegment superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            nucleotide sequence

        Examples
        --------
        >>> fw3 = Framework3NT("TACTACGCCGATAGCGTAAAGGGGAGATTTACCATCTCACGAGATAACGCCAAGAACACCTTGTACCTGCAGATGAACTCCTTGAAGGCCGAGGACACAGCCACACACTACTGT")
        """
        super().__init__(seq)


class FrameWork4NT(AntibodySegmentNT):
    """Framework 4 nucelotide segments

    Parameters
    ----------
    AntibodySegmentNT :
       AntibodySegment nucleotide superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            nucleotide sequence

        Examples
        --------
        >>> fw4 = Framework4NT("TACTACGCCGATAGCGTAAAGGGGAGATTTACCATCTCACGAGATAACGCCAAGAACACCTTGTACCTGCAGATGAACTCCTTGAAGGCCGAGGACACAGCCACACACTACTGT")
        """
        super().__init__(seq)


class CDR1NT(AntibodySegmentNT):
    """CDR1 Nucleotide Sequence Class

    Parameters
    ----------
    AntibodySegmentNT :
        Antibody nucletodie segment superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            nucleotide sequence

        Examples
        --------
        >>> cdr1_nt = CDR1NT("GGCCTGTCCCTGACTTCTGGGTCT")
        """
        super().__init__(seq)


class CDR2NT(AntibodySegmentNT):
    """CDR2 Nucleotide Sequence Class

    Parameters
    ----------
    AntibodySegmentNT :
        Antibody nucletodie segment superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            nucleotide sequence

        Examples
        --------
        >>> cdr2_nt = CDR2NT("GGCCTGTCCCTGACTTCTGGGTCT")
        """
        super().__init__(seq)


class CDR3NT(AntibodySegmentNT):
    """CDR3 Nucleotide Sequence Class

    Parameters
    ----------
    AntibodySegmentNT :
        Antibody nucletodie segment superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            nucleotide sequence

        Examples
        --------
        >>> cdr2_nt = CDR3NT("GGCCTGTCCCTGACTTCTGGGTCT")
        """
        super().__init__(seq)
        self.germline = ""

    @property
    def germline(self) -> str:
        """get the germline sequence if set

        Returns
        -------
        str
            the germline sequence
        """
        return self._germline

    @germline.setter
    def germline(self, germ):
        """set the germline sequence of the HCDR3 based on the v and j gene.


        Parameters
        ----------
        germ : (v_gene, j_gene)
            a tuple of v gene and jgene strings representing the CDR3 parts

        Examples
        --------
        >>>cdr3_nt = CDR3NT("GGCCTGTCCCTGACTTCTGGGTCT")
        >>>cdr3_nt.germline = ("GGC","GAAACT")
        >>>cdr3_nt.germline
        >>>"GGC------GAACT"
        Raises
        ------
        ValueError
            if the the germline sequence is not a tuple of strings

        """
        if not isinstance(germ, (tuple, str)):
            raise ValueError(f"{germ} must be instance of tuple of strings")
        if isinstance(germ, str):
            self._germline = germ
        elif isinstance(germ, tuple):
            _v, _j = germ
            _dashes = len(self.sequence) - len(_v) - len(_j)
            self._germline = str(_v) + "-" * _dashes + str(_j)
            self._v_portion = self.sequence[0 : len(_v)]
        if self._germline:
            if not self._sequence:
                _sequence = "-" * len(self._germline)
            else:
                _sequence = self._sequence

            _alignments = align.globalxs(self._germline, _sequence, -10, -1)
            self._alignment = _alignments[0]

    @property
    def v_portion(self) -> str:
        """The V region of the nucleotide"""
        return self._v_portion

    @v_portion.setter
    def v_portion(self, v_portion: str):
        self._v_portion = v_portion


class FrameWork1AA(AntibodySegmentAA):
    """Framework 1 amino acid segments

    Parameters
    ----------
    AntibodySegmentAA :
       AntibodySegment AA superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            amino acid sequence

        Examples
        --------
        >>> fw1_aa = FrameWork1AA("DIQMTQSPASLSASLGETVSIECLAS")
        """
        super().__init__(seq)


class FrameWork2AA(AntibodySegmentAA):
    """Framework 2 amino acid segments

    Parameters
    ----------
    AntibodySegmentAA :
       AntibodySegment AA superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            amino acid sequence

        Examples
        --------
        >>> fw2_aa = FrameWork2AA("LAWYQLKPGKSPQFLI")
        """
        super().__init__(seq)


class FrameWork3AA(AntibodySegmentAA):
    """Framework 3 amino acid segments

    Parameters
    ----------
    AntibodySegmentAA :
       AntibodySegment AA superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            amino acid sequence

        Examples
        --------
        >>> fw2_aa = FrameWork2AA("SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVYYC")
        """
        super().__init__(seq)


class FrameWork4AA(AntibodySegmentAA):
    """Framework 4 amino acid segments

    Parameters
    ----------
    AntibodySegmentAA :
       AntibodySegment AA superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            amino acid sequence

        Examples
        --------
        >>> fw4_aa = FrameWork4AA("FGSGTKLKIK")
        """
        super().__init__(seq)


class CDR1AA(AntibodySegmentAA):
    """CDR1 Amino Acid Sequence Class

    Parameters
    ----------
    AntibodySegmentAA :
        Antibody AA segment superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            amino acid sequence

        Examples
        --------
        >>> cdr1_aa = CDR1AA("DIQMTQSPASLSS")
        """
        super().__init__(seq)


class CDR2AA(AntibodySegmentAA):
    """CDR2 Amino Acid Sequence Class

    Parameters
    ----------
    AntibodySegmentAA :
        Antibody AA segment superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            amino acid sequence

        Examples
        --------
        >>> cdr2_aa = CDR2AA("DIQMTQSPASLSS")
        """
        super().__init__(seq)


class CDR3AA(AntibodySegmentAA):
    """CDR3 Amino Acid Sequence Class

    Parameters
    ----------
    AntibodySegmentAA :
        Antibody amino acid segment superclass
    """

    def __init__(self, seq):
        """constructor

        Parameters
        ----------
        seq : str
            amino acid sequence

        Examples
        --------
        >>> cdr3_aa = CDR3AA("VRQSPASV")
        """
        super().__init__(seq)
        self.germline = ""

    @property
    def germline(self) -> str:
        """get the germline sequence if set

        Returns
        -------
        str
            the germline sequence
        """
        return self._germline

    @germline.setter
    def germline(self, germ):
        """set the germline sequence of the HCDR3 based on the v and j gene.


        Parameters
        ----------
        germ : (v_gene, j_gene)
            a tuple of v gene and jgene strings representing the CDR3 parts

        Examples
        --------
        >>>cdr3_aa = CDR3AA("VRQSPASV")
        >>>cdr3_aa.germline = ("V","YV")
        >>>cdr3_aa.germline
        >>>"V----YV"
        Raises
        ------
        ValueError
            if the the germline sequence is not a tuple of strings

        """
        if not isinstance(germ, (tuple, str)):
            raise ValueError(f"{germ} must be instance of tuple of strings")
        if isinstance(germ, str):
            self._germline = germ
        elif isinstance(germ, tuple):
            _v, _j = germ
            _dashes = len(self.sequence) - len(_v) - len(_j)
            self._germline = str(_v) + "-" * _dashes + str(_j)

        if self._germline:
            if not self._sequence:
                _sequence = "-" * len(self._germline)
            else:
                _sequence = self._sequence
            _alignments = align.globalxs(self._germline, _sequence, -10, -1)
            self._alignment = _alignments[0]
