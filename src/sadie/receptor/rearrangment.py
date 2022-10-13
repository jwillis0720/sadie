"""
This file contains objects to represent the rearrangemnt scheme
https://docs.airr-community.org/en/stable/datarep/rearrangements.html
"""
import re
from functools import lru_cache
from typing import Any, List, Optional, Set, Union
from uuid import UUID, uuid4

from Bio.Seq import Seq
from pandas._libs.missing import NAType
from pydantic import BaseModel, validator


@lru_cache(maxsize=1)
def get_nt_validator_regex() -> Any:
    return re.compile(r"^[ACNTGYRWSKMDVHBacntgyrwskmdvhb]+\Z")


@lru_cache(maxsize=1)
def get_aa_validator_regex() -> Any:
    return re.compile(r"^[ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy]+\Z")


class RearrargmentCategory(BaseModel):
    """
    The category of rearrangement annotation object

    Atttributes
    -----------

    category: str
        The category of rearrangement objectect:

    Categories
    ----------
        Input:
            The input sequence to the V(D)J assignment process.

        Identifiers:
            Primary and foreign key identifiers for linking AIRR data across files and databases.

        Primary Annotations:
            The primary outputs of the V(D)J assignment process, which includes the gene locus, V, D, J, and C gene calls, various flags, V(D)J junction sequence, copy number (duplicate_count), and the number of reads contributing to a consensus input sequence (consensus_count).

        Alignment Annotations:
            Detailed alignment annotations including the input and germline sequences used in the alignment; score, identity, statistical support (E-value, likelihood, etc); and the alignment itself through CIGAR strings for each aligned gene.

        Alignment Positions:
            The start/end positions for genes in both the input and germline sequences.

        Region Sequence:
            Sequence annotations for the framework regions (FWRs) and complementarity-determining regions (CDRs).

        Region Positions:
            Positional annotations for the framework regions (FWRs) and complementarity-determining regions (CDRs).

        Junction Lengths:
            Lengths for junction sub-regions associated with aspects of the V(D)J recombination process.
    """

    category: str

    @validator("category")
    @classmethod
    def validate_category(cls, v: str) -> str:
        valid_categories: Set[str] = {
            "input",
            "identifiers",
            "primary_annotations",
            "alignment_annotations",
            "alignment_positions",
            "region_sequence_annotations",
            "region_positions",
            "junction_lengths",
        }
        if v not in valid_categories:
            raise ValueError(f"{v} is not a valid category, use {valid_categories}")
        return v


class InputSequence(BaseModel):
    """
    These required VDJ Sequences are taken from https://docs.airr-community.org/en/stable/datarep/rearrangements.html

    Attributes - AIRR 1.3
    ---------------------
    sequence_id: Optional[Union[str,UUID]]:
        Unique query sequence identifier for the Rearrangment. Most often this will be the input sequence header or a substring thereof, but may also be a custom identifier defined by the tool in cases where query sequences have been combined in some fashion prior to alignment. If not, given explicitly, will be randomly generated a UUID

    sequence: Union[Seq,str]
        The sequence vdj_c sequence or the rearrangment that starts at the first nt of the alignment. It is reverse complimented if necessary.

    sequence_aa: Optional[Union[Seq,str]]:
        Amino acid translation of the nucleotide sequence (not the raw sequence).

    Attributes - non-AIRR
    ---------------------
    raw_sequence: Optional[Union[Seq,str]]:
        The unmodified query sequence
    category: Optional[RearrargmentCategory]:
        The category of rearrangement objectect
    """

    sequence_id: Optional[Union[str, UUID]]
    sequence: Union[Seq, str]
    raw_sequence: Optional[Union[Seq, str]] = None  # non airr
    sequence_aa: Optional[Union[Seq, str]] = None
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="input")

    @validator("sequence_id", always=True)
    @classmethod
    def validate_sequence_id(cls, v: Optional[Union[str, UUID]]) -> str:
        """If no sequence id is provided, use a unique UUID"""
        if not v:
            v = uuid4()
        return str(v)

    @validator("sequence")
    @classmethod
    def validate_sequence(cls, v: Union[Seq, str]) -> Seq:
        """Check for valid nt sequence"""
        nt_validator: re.Pattern[Any] = get_nt_validator_regex()
        if isinstance(v, Seq):
            v = str(v)
        if not nt_validator.match(v):
            raise ValueError(f"{v} is not a valid nucleotide sequence, must only contain ACGTN")
        else:
            return Seq(v)

    @validator("raw_sequence")
    @classmethod
    def validate_raw_sequence(cls, v: Union[Seq, str]) -> Seq:
        """Check for valid nt sequence"""
        nt_validator: re.Pattern[Any] = get_nt_validator_regex()
        if isinstance(v, Seq):
            v = str(v)
        if not nt_validator.match(v):
            raise ValueError(f"{v} is not a valid nucleotide sequence, must only contain ACGTN")
        else:
            return Seq(v)

    @validator("sequence_aa")
    @classmethod
    def validate_sequence_aa(cls, v: Union[Seq, str]) -> Union[None, Seq]:
        aa_validator: re.Pattern[Any] = get_aa_validator_regex()
        if isinstance(v, Seq):
            v = str(v)
        if not aa_validator.match(v):
            raise ValueError(f"{v} is not a valid amino acid sequence, must only contain ACDEFGHIKLMNPQRSTVWXY")
        return Seq(v)

    @staticmethod
    def get_airr_fields() -> List[str]:
        return ["sequence_id", "sequence"]

    class Config:
        arbitrary_types_allowed = True


class PrimaryAnnotations(BaseModel):
    """
    The primary outputs of the V(D)J assignment process, which includes the gene locus, V, D, J, and C gene calls, various flags, V(D)J junction sequence, copy number (duplicate_count), and the number of reads contributing to a consensus input sequence (consensus_count). Taken from https://docs.airr-community.org/en/stable/datarep/rearrangements.html


    Attributes - AIRR 1.3
    ---------------------
    rev_comp: bool
        True if the alignment is on the opposite strand (reverse complemented) with respect to the query sequence. If True then all output data, such as alignment coordinates and sequences, are based on the reverse complement of ‘sequence’.
    productive: bool
        True if the V(D)J sequence is predicted to be productive:
            1. Coding region has an open reading frame
            2. No defect in the start codon, splicing sites or regulatory elements.
            3. No internal stop codons.
            4. An in-frame junction region.
    vj_in_frame bool:
        True if the V and J gene alignments are in-frame
    stop_codon: Optional[bool]
        True if the aligned sequence contains a stop codon.
    complete_vdj: Optional[bool]
        True if the sequence alignment spans the entire V(D)J region. Meaning, sequence_alignment includes both the first V gene codon that encodes the mature polypeptide chain (i.e., after the leader sequence) and the last complete codon of the J gene (i.e., before the J:C splice site). This does not require an absence of deletions within the internal FWR and CDR regions of the alignment.
    locus: Optional[str]
        Gene locus (chain type). Note that this field uses a controlled vocabulary that is meant to provide a generic classification of the locus, not necessarily the correct designation according to a specific nomenclature. here are the loci IGH, IGK, IGL, TRA, TRB, TRD, or TRG.
    v_call: Union[str, List[str]]
        V gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHV4-59*01 if using IMGT/GENE-DB).
    d_call: Optional[str]
        First or only D gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHD3-10*01 if using IMGT/GENE-DB).
    d2_call: Optional[str]
        Second D gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHD3-10*01 if using IMGT/GENE-DB).
    j_call: Union[str, List[str]]
        J gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHJ4*02 if using IMGT/GENE-DB).
    c_call: Optional[str]
        Constant region gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHG1*01 if using IMGT/GENE-DB).

    Attributes - Non-AIRR
    ---------------------
    v_call_top: Optional[str]
        The top V gene call in a comma seperated list of V gene calls.
    v_call_top_gene: Optional[str]
        The top V gene call without the allele. (ex. IGHV1-69)
    v_call_top_allele: Optional[str]
        The top V gene call alleld. (ex. *01)
    d_call_top: Optional[str]
        The top D gene call in a comma seperated list of D gene calls.
    d_call_gene: Optional[str]
        The top D gene call without the allele. (ex. IGHD3-3)
    d_call_allele: Optional[str]
        The top D gene call alleld. (ex. *01)
    j_call_top: Optional[str]
        The top J gene call in a comma seperated list of J gene calls.
    j_call_top_gene: Optional[str]  # Non Airr
        The top J gene call without the allele. (ex. IGHJ4)
    j_call_top_allele: Optional[str]  # Non Airr
        The top J gene call alleld. (ex. *02)
    c_call_allele: Optional[str]  # Non Airr
        The top C gene call alleld. (ex. *01)
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="primary_annotations")
        The category of the rearrangement object
    """

    # Airr fields
    rev_comp: bool
    productive: bool
    vj_in_frame: Optional[bool] = None
    stop_codon: Optional[bool] = None
    complete_vdj: Optional[bool] = None
    locus: Optional[str] = None
    v_call: Union[str, List[str]]
    d_call: Optional[Union[str, List[str]]] = None
    d2_call: Optional[str] = None
    j_call: Union[str, List[str]]
    c_call: Optional[str] = None

    # Non Airr fields
    v_call_top: Optional[str] = None
    v_call_top_gene: Optional[str] = None
    v_call_top_allele: Optional[str] = None
    d_call_top: Optional[str] = None
    d_call_gene: Optional[str] = None
    d_call_allele: Optional[str] = None
    j_call_top: Optional[str] = None
    j_call_top_gene: Optional[str] = None
    j_call_top_allele: Optional[str] = None
    c_call_allele: Optional[str] = None
    reference_name: Optional[str] = None
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="primary_annotations")

    @staticmethod
    def get_airr_fields() -> List[str]:
        return [
            "rev_comp",
            "productive",
            "vj_in_frame",
            "stop_codon",
            "complete_vdj",
            "locus",
            "reference_name",
            "v_call",
            "d_call",
            "j_call",
            "v_call_top",
            "d_call_top",
            "j_call_top",
        ]

    class Config:
        arbitrary_types_allowed = True


class AlignmentAnnotations(BaseModel):
    """
    Detailed alignment annotations including the input and germline sequences used in the alignment; score, identity, statistical support (E-value, likelihood, etc); and the alignment itself through CIGAR strings for each aligned gene.

    Attributes - AIRR 1.3:
    ---------------------
    sequecne_alignment: Union[str, Seq]
       Aligned portion of query sequence, including any indel corrections or numbering spacers, such as IMGT-gaps. Typically, this will include only the V(D)J region, but that is not a requirement.
    seqeunce_alignment_aa: Optional[Union[str, Seq]]
        Amino acid translation of the aligned query sequence.
    germline_alignement: Union[str, Seq]
        Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and spacers (if any)
    germline_alignment_aa: Optional[Union[str, Seq]]
        Amino acid translation of the aligned germline sequence.
    v_score: Optional[float]
        Alignment score for the V gene.
    v_identity: Optional[float]
        Fractional identity of the V gene alignment
    v_support: Optional[float]
        V gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the V gene assignment as defined by the alignment tool.
    v_cigar: str
        CIGAR string for the V gene alignment.
    d_score: Optional[float]
        Alignment score for the D gene.
    d_identity: Optional[float]
        Fractional identity of the D gene alignment
    d_support: Optional[float]
        D gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the D gene assignment as defined by the alignment tool.
    d_cigar: str
        CIGAR string for the D gene alignment.
    d2_score: Optional[float]
        Alignment score for the second D gene.
    d2_identity: Optional[float]
        Fractional identity of the second D gene alignment
    d2_support: Optional[float]
        Second D gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the second D gene assignment as defined by the alignment tool.
    d2_cigar: Optional[str]
        CIGAR string for the second D gene alignment.
    j_score: Optional[float]
        Alignment score for the J gene.
    j_identity: Optional[float]
        Fractional identity of the J gene alignment
    j_support: Optional[float]
        J gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the J gene assignment as defined by the alignment tool.
    j_cigar: str
        CIGAR string for the J gene alignment.
    junction: Union[str, Seq]
        Junction region nucleotide sequence, where the junction is defined as the CDR3 plus the two flanking conserved codons.
    junction_aa: Optional[Union[str, Seq]]
        Amino acid translation of the junction.
    np1: Optional[Union[str, Seq]]
        Nucleotide sequence of the combined N/P region between the V gene and first D gene alignment or between the V gene and J gene alignments.
    np1_aa: Optional[Union[str, Seq]]
        Amino acid translation of the np1 field.
    np2: Optional[Union[str, Seq]]
        Nucleotide sequence of the combined N/P region between either the first D gene and J gene alignments or the first D gene and second D gene alignments.
    np2_aa: Optional[Union[str, Seq]]
        Amino acid translation of the np2 field.
    np3: Optional[Union[str, Seq]]
        Nucleotide sequence of the combined N/P region between the second D gene and J gene alignments.
    np3_aa: Optional[Union[str, Seq]]
        Amino acid translation of the np3 field.
    c_score: Optional[float]
        Alignment score for the C gene alignment.
    c_identity: Optional[float]
        Fractional identity of the C gene alignment
    c_support: Optional[float]
        C gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the C gene assignment as defined by the alignment tool.
    c_cigar: Optional[str]
        CIGAR string for the C gene alignment.

    Attributes - Non AIRR:
    ---------------------
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="alignment_annotations")
        The category of the rearrangement object
    """

    sequence_alignment: Union[str, Seq]
    sequence_alignment_aa: Optional[Union[str, Seq]] = None
    germline_alignment: Union[str, Seq]
    germline_alignment_aa: Optional[Union[str, Seq]] = None
    v_score: Optional[float] = None
    v_identity: Optional[float] = None
    v_support: Optional[float] = None
    v_cigar: str
    d_score: Optional[float] = None
    d_identity: Optional[float] = None
    d_support: Optional[float] = None
    d_cigar: str
    d2_score: Optional[float] = None
    d2_identity: Optional[float] = None
    d2_support: Optional[float] = None
    d2_cigar: Optional[str] = None
    j_score: Optional[float] = None
    j_identity: Optional[float] = None
    j_support: Optional[float] = None
    j_cigar: str
    junction: Union[str, Seq]
    junction_aa: Optional[Union[str, Seq]] = None
    np1: Optional[Union[str, Seq]] = None
    np1_aa: Optional[Union[str, Seq]] = None
    np2: Optional[Union[str, Seq]] = None
    np2_aa: Optional[Union[str, Seq]] = None
    np3: Optional[Union[str, Seq]] = None
    np3_aa: Optional[Union[str, Seq]] = None
    c_score: Optional[float] = None
    c_identity: Optional[float] = None
    c_support: Optional[float] = None
    c_cigar: Optional[str] = None

    # Non Airr
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="alignment_annotations")

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def get_airr_fields() -> List[str]:
        return [
            "sequence_alignment",
            "sequence_alignment_aa",
            "germline_alignment",
            "germline_alignment_aa",
            "v_score",
            "v_identity",
            "v_support",
            "v_cigar",
            "d_score",
            "d_identity",
            "d_support",
            "d_cigar",
            "j_score",
            "j_identity",
            "j_support",
            "j_cigar",
            "junction",
            "junction_aa",
            "np1",
            "np2",
        ]


class AlignmentPositions(BaseModel):
    """
    The start/end positions for genes in both the input and germline sequences.


    Attributes:
    ----------
    v_sequence_start: Optional[int]
        Start position of the V gene in the query sequence (1-based closed interval).

    v_sequence_end: Optional[int]
        End position of the V gene in the query sequence (1-based closed interval).

    v_germline_start: Optional[int]
        Alignment start position in the V gene reference sequence (1-based closed interval).

    v_germline_end: Optional[int]
        Alignment end position in the V gene reference sequence (1-based closed interval).

    v_alignment_start: Optional[int]:
        Start position of the V gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).

    v_alignment_end: Optional[int]:
        End position of the V gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).

    d_sequence_start: Optional[int]:
        Start position of the first or only D gene in the query sequence. (1-based closed interval).

    d_sequence_end: Optional[int]
        End position of the first or only D gene in the query sequence. (1-based closed interval).

    d_germline_start: Optional[int]
        Alignment start position in the D gene reference sequence for the first or only D gene (1-based closed interval).

    d_germline_end: Optional[int]
        Alignment end position in the D gene reference sequence for the first or only D gene (1-based closed interval).

    d_alignment_start: Optional[int]
        Start position of the first or only D gene in both the sequence_alignment and germline_alignment fields (1-based closed interval).

    d_alignment_end: Optional[int]
        End position of the first or only D gene in both the sequence_alignment and germline_alignment fields (1-based closed interval).

    d2_sequence_start: Optional[int]
        Start position of the second D gene in the query sequence (1-based closed interval).

    d2_sequence_end: Optional[int]
        End position of the second D gene in the query sequence (1-based closed interval).

    d2_germline_start: Optional[int]
        Alignment start position in the second D gene reference sequence (1-based closed interval).

    d2_germline_end: Optinal[int]
        Alignment end position in the second D gene reference sequence (1-based closed interval).

    d2_alignment_start: Optinal[int]
        Start position of the second D gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).

    d2_alignment_end: Optional[int]
        End position of the second D gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).

    j_sequence_start: Optional[int]
        Start position of the J gene in the query sequence (1-based closed interval).

    j_sequence_end: Opitional[int]
        End position of the J gene in the query sequence (1-based closed interval).

    j_germline_start: Optional[int]
        Alignment start position in the J gene reference sequence (1-based closed interval).

    j_germline_end: Optional[int]
        Alignment end position in the J gene reference sequence (1-based closed interval).

    j_alignment_start: Optional[int]
        Start position of the J gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).

    j_alignment_end: Optional[int]
        End position of the J gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).
    """

    v_sequence_start: Optional[Union[int, NAType]] = None
    v_sequence_end: Optional[Union[int, NAType]] = None
    v_germline_start: Optional[Union[int, NAType]] = None
    v_germline_end: Optional[Union[int, NAType]] = None
    v_alignment_start: Optional[Union[int, NAType]] = None
    v_alignment_end: Optional[Union[int, NAType]] = None
    d_sequence_start: Optional[Union[int, NAType]] = None
    d_sequence_end: Optional[Union[int, NAType]] = None
    d_germline_start: Optional[Union[int, NAType]] = None
    d_germline_end: Optional[Union[int, NAType]] = None
    d_alignment_start: Optional[Union[int, NAType]] = None
    d_alignment_end: Optional[Union[int, NAType]] = None
    d2_sequence_start: Optional[Union[int, NAType]] = None
    d2_sequence_end: Optional[Union[int, NAType]] = None
    d2_germline_start: Optional[Union[int, NAType]] = None
    d2_germline_end: Optional[Union[int, NAType]] = None
    d2_alignment_start: Optional[Union[int, NAType]] = None
    d2_alignment_end: Optional[Union[int, NAType]] = None
    j_sequence_start: Optional[int] = None
    j_sequence_end: Optional[int] = None
    j_germline_start: Optional[int] = None
    j_germline_end: Optional[int] = None
    j_alignment_start: Optional[int] = None
    j_alignment_end: Optional[int] = None
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="alignment_positions")

    @validator(
        "v_sequence_start",
        "v_sequence_end",
        "v_germline_start",
        "v_germline_end",
        "v_alignment_start",
        "v_alignment_end",
        "d_sequence_start",
        "d_sequence_end",
        "d_germline_start",
        "d_germline_end",
        "d_alignment_start",
        "d_alignment_end",
        "d2_sequence_start",
        "d2_sequence_end",
        "d2_germline_start",
        "d2_germline_end",
        "d2_alignment_start",
        "d2_alignment_end",
        "j_sequence_start",
        "j_sequence_end",
        "j_germline_start",
        "j_germline_end",
        "j_alignment_start",
        "j_alignment_end",
    )
    @classmethod
    def validate_with_na(cls, v: Union[int, NAType]) -> Union[int, None]:
        if v is None or isinstance(v, int):
            return v
        if isinstance(v, NAType):
            return None
        raise ValueError(f"Invalid value for alignment_positions: {v}")

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def get_airr_fields() -> List[str]:
        return [
            "v_sequence_start",
            "v_sequence_end",
            "v_germline_start",
            "v_germline_end",
            "v_alignment_start",
            "v_alignment_end",
            "d_sequence_start",
            "d_sequence_end",
            "d_germline_start",
            "d_germline_end",
            "d_alignment_start",
            "d_alignment_end",
            "j_sequence_start",
            "j_sequence_end",
            "j_germline_start",
            "j_germline_end",
            "j_alignment_start",
            "j_alignment_end",
        ]


class RegionSequences(BaseModel):
    """
    Sequence annotations for the framework regions (FWRs) and complementarity-determining regions (CDRs).

    Attributes
    ----------
    fwr1: Optional[Union[str,Seq]]
        Nucleotide sequence of the aligned FWR1 region.
    fwr1_aa: Optional[Union[str,Seq]]
        Amino acid translation of the fwr1 field.
    cdr1: Optional[Union[str,Seq]]
        Nucleotide sequence of the aligned CDR1 region.
    cdr1_aa: Optional[Union[str,Seq]]
        Amino acid translation of the cdr1 field.
    fwr2: Optional[Union[str,Seq]]
        Nucleotide sequence of the aligned FWR2 region.
    fwr2_aa: Optional[Union[str,Seq]]
        Amino acid translation of the fwr2 field.
    cdr2: Optional[Union[str,Seq]]
        Nucleotide sequence of the aligned CDR2 region.
    cdr2_aa: Optional[Union[str,Seq]]
        Amino acid translation of the cdr2 field.
    fwr3: Optional[Union[str,Seq]]
        Nucleotide sequence of the aligned FWR3 region.
    fwr3_aa: Optional[Union[str,Seq]]
        Amino acid translation of the fwr3 field.
    cdr3: Optional[Union[str,Seq]]
        Nucleotide sequence of the aligned CDR3 region.
    cdr3_aa: Optional[Union[str,Seq]]
        Amino acid translation of the cdr3 field.
    fwr4: Optional[Union[str,Seq]]
        Nucleotide sequence of the aligned FWR4 region.
    fwr4_aa: Optional[Union[str,Seq]]
        Amino acid translation of the fwr4 field.
    """

    fwr1: Optional[Union[str, Seq]] = None
    fwr1_aa: Optional[Union[str, Seq]] = None
    cdr1: Optional[Union[str, Seq]] = None
    cdr1_aa: Optional[Union[str, Seq]] = None
    fwr2: Optional[Union[str, Seq]] = None
    fwr2_aa: Optional[Union[str, Seq]] = None
    cdr2: Optional[Union[str, Seq]] = None
    cdr2_aa: Optional[Union[str, Seq]] = None
    fwr3: Optional[Union[str, Seq]] = None
    fwr3_aa: Optional[Union[str, Seq]] = None
    cdr3: Optional[Union[str, Seq]] = None
    cdr3_aa: Optional[Union[str, Seq]] = None
    fwr4: Optional[Union[str, Seq]] = None
    fwr4_aa: Optional[Union[str, Seq]] = None
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="region_sequence_annotations")

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def get_airr_fields() -> List[str]:
        return [
            "fwr1",
            "fwr1_aa",
            "cdr1",
            "cdr1_aa",
            "fwr2",
            "fwr2_aa",
            "cdr2",
            "cdr2_aa",
            "fwr3",
            "fwr3_aa",
            "cdr3",
            "cdr3_aa",
            "fwr4",
            "fwr4_aa",
        ]


class RegionPositions(BaseModel):
    """
    Positional annotations for the framework regions (FWRs) and complementarity-determining regions (CDRs).

    Attributes
    ----------
    fwr1_start: Optinal[int]
        FWR1 start position in the query sequence (1-based closed interval).
    fwr1_end: Optional[int]
        FWR1 end position in the query sequence (1-based closed interval).
    cdr1_start : Optional[int]
        CDR1 start position in the query sequence (1-based closed interval).
    cdr1_end: Optional[int]
        CDR1 end position in the query sequence (1-based closed interval).
    fwr2_start: Optinal[int]
        FWR2 start position in the query sequence (1-based closed interval).
    fwr2_end: Optional[int]
        FWR2 end position in the query sequence (1-based closed interval).
    cdr2_start : Optional[int]
        CDR2 start position in the query sequence (1-based closed interval).
    cdr2_end: Optional[int]
        CDR2 end position in the query sequence (1-based closed interval).
    fwr3_start: Optinal[int]
        FWR3 start position in the query sequence (1-based closed interval).
    fwr3_end: Optional[int]
        FWR3 end position in the query sequence (1-based closed interval).
    cdr3_start : Optional[int]
        CDR3 start position in the query sequence (1-based closed interval).
    cdr3_end: Optional[int]
        CDR3 end position in the query sequence (1-based closed interval).
    fwr4_start: Optinal[int]
        FWR4 start position in the query sequence (1-based closed interval).
    fwr4_end: Optional[int]
        FWR4 end position in the query sequence (1-based closed interval).
    """

    fwr1_start: Optional[int] = None
    fwr1_end: Optional[int] = None
    cdr1_start: Optional[int] = None
    cdr1_end: Optional[int] = None
    fwr2_start: Optional[int] = None
    fwr2_end: Optional[int] = None
    cdr2_start: Optional[int] = None
    cdr2_end: Optional[int] = None
    fwr3_start: Optional[int] = None
    fwr3_end: Optional[int] = None
    cdr3_start: Optional[int] = None
    cdr3_end: Optional[int] = None
    fwr4_start: Optional[int] = None
    fwr4_end: Optional[int] = None
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="region_positions")

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def get_airr_fields() -> List[str]:
        return [
            "fwr1_start",
            "fwr1_end",
            "cdr1_start",
            "cdr1_end",
            "fwr2_start",
            "fwr2_end",
            "cdr2_start",
            "cdr2_end",
            "fwr3_start",
            "fwr3_end",
            "cdr3_start",
            "cdr3_end",
            "fwr4_start",
            "fwr4_end",
        ]


class JunctionLengths(BaseModel):
    """
    Lengths for junction sub-regions associated with aspects of the V(D)J recombination process.

    Attributes
    ----------
    junction_length: Optional[int]
        Number of nucleotides in the junction sequence.
    junction_aa_length: Optional[int]
        Number of amino acids in the junction sequence.
    np1_length: Optinal[int]
        Number of nucleotides between the V gene and first D gene alignments or between the V gene and J gene alignments.
    np2_length: Optinal[int]
        Number of nucleotides between either the first D gene and J gene alignments or the first D gene and second D gene alignments.
    np3_length: Optinal[int]
        Number of nucleotides between the second D gene and J gene alignments.
    n1_length: Optinal[int]
        Number of untemplated nucleotides 5’ of the first or only D gene alignment.
    n2_length: Optinal[int]
        Number of untemplated nucleotides 3’ of the first or only D gene alignment.
    n3_length: Optinal[int]
        Number of untemplated nucleotides 3’ of the second D gene alignment.
    p3v_length: Optinal[int]
        Number of palindromic nucleotides 3’ of the V gene alignment.
    p5d_length: Optinal[int]
        Number of palindromic nucleotides 5’ of the first or only D gene alignment.
    p3d_length: Optinal[int]
        Number of palindromic nucleotides 3’ of the first or only D gene alignment.
    p5d2_length: Optinal[int]
        Number of palindromic nucleotides 5’ of the second D gene alignment.
    p3d2_length: Optinal[int]
        Number of palindromic nucleotides 3’ of the second D gene alignment.
    p5j_length: Optinal[int]
        Number of palindromic nucleotides 5’ of the J gene alignment.
    """

    junction_length: Optional[int] = None
    junction_aa_length: Optional[int] = None
    np1_length: Optional[int] = None
    np2_length: Optional[int] = None
    np3_length: Optional[int] = None
    n1_length: Optional[int] = None
    n2_length: Optional[int] = None
    n3_length: Optional[int] = None
    p3v_length: Optional[int] = None
    p5d_length: Optional[int] = None
    p3d_length: Optional[int] = None
    p5d2_length: Optional[int] = None
    p3d2_length: Optional[int] = None
    p5j_length: Optional[int] = None
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="junction_lengths")

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def get_airr_fields() -> List[str]:
        return [
            "junction_length",
            "junction_aa_length",
        ]


class ReceptorChain(BaseModel):
    input_sequence: InputSequence
    primary_annotations: PrimaryAnnotations
    alignment_annotations: AlignmentAnnotations
    alignment_positions: Optional[AlignmentPositions] = None
    region_sequences: Optional[RegionSequences] = None
    region_positions: Optional[RegionPositions] = None
    junction_lengths: Optional[JunctionLengths] = None

    @staticmethod
    def from_single(
        sequence_id: str,
        sequence: Union[str, Seq],
        reference_name: str = "human",
        database: str = "imgt",
    ) -> "ReceptorChain":
        """
        Create a receptor chain from a single sequence.

        Parameters
        ----------
        sequence_id: str
            Identifier for the sequence.
        sequence: Union[str,Seq]
            Sequence data.

        Returns
        -------
        receptor_chain: ReceptorChain
            Receptor chain.
        """
        from sadie.airr import Airr
        from sadie.airr.airrtable import AirrSeries, AirrTable

        airr_api = Airr(reference_name)
        result: AirrTable = airr_api.run_single(sequence_id, str(sequence))
        result_sliced: AirrSeries = result.iloc[0]  # type: ignore

        # my py won't stop complaining unless I pass the sliced object back through itself
        return ReceptorChain(**result_sliced.to_receptor_chain_object().__dict__)

    def __str__(self) -> str:
        printable: List[str] = []
        for key in self.__fields__.keys():
            sub_obj = self.__dict__[key].__dict__
            printable.append(key)
            printable.append("-" * len(key))
            for sub_key in sub_obj.keys():
                printable.append(f"{sub_key} : {sub_obj[sub_key]}")
            printable.append("\n")
        return "\n".join(printable)


class Antibody(BaseModel):
    heavy_chain: ReceptorChain
    light_chain: ReceptorChain


class Antibodies(BaseModel):
    antibodies: List[Antibody]


class TCR(BaseModel):
    heavy_chain: ReceptorChain
    light_chain: ReceptorChain
