"""
This file contains objects to represent the rearrangemnt scheme
https://docs.airr-community.org/en/stable/datarep/rearrangements.html
"""
from functools import lru_cache
from pydantic import BaseModel, validator
from typing import Any, Dict, List, Optional, Union, Set
from Bio.Seq import Seq
from uuid import UUID, uuid4
import re


@lru_cache
def get_nt_validator_regex() -> re.Pattern[Any]:
    return re.compile(r"^[ACNTGacntg]+\Z")


@lru_cache
def get_aa_validator_regex() -> re.Pattern[Any]:
    return re.compile(r"^[ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy]+\Z")


class RearrargmentCategory(BaseModel):
    """
    The category of rearrangement

    Input: The input sequence to the V(D)J assignment process.

    Identifiers: Primary and foreign key identifiers for linking AIRR data across files and databases.

    Primary Annotations : The primary outputs of the V(D)J assignment process, which includes the gene locus, V, D, J, and C gene calls, various flags, V(D)J junction sequence, copy number (duplicate_count), and the number of reads contributing to a consensus input sequence (consensus_count).

    Alignment Annotations: Detailed alignment annotations including the input and germline sequences used in the alignment; score, identity, statistical support (E-value, likelihood, etc); and the alignment itself through CIGAR strings for each aligned gene.

    Alignment Positions: The start/end positions for genes in both the input and germline sequences.

    Region Sequence Sequence annotations for the framework regions (FWRs) and complementarity-determining regions (CDRs).

    Region Positions:Positional annotations for the framework regions (FWRs) and complementarity-determining regions (CDRs).

    Junction Lengths:Lengths for junction sub-regions associated with aspects of the V(D)J recombination process.
    """

    category: str

    @validator("category")
    @classmethod
    def validate_category(cls, v: str):
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


class VDJSequence(BaseModel):
    """
    These required VDJ Sequences are taken from https://docs.airr-community.org/en/stable/datarep/rearrangements.html

    Attributes
    ----------
    sequence_id: Optional[Union[str,UUID]]:
        Unique query sequence identifier for the Rearrangment. Most often this will be the input sequence header or a substring thereof, but may also be a custom identifier defined by the tool in cases where query sequences have been combined in some fashion prior to alignment. If not, given explicitly, will be randomly generated a UUID

    sequence: Union[Seq,str]
        The sequence vdj_c sequence or the rearrangment that starts at the first nt of the alignment. It is reverse complimented if necessary.

    raw_sequence: Optional[Union[Seq,str]]:
        The unmodified query sequence
    """

    sequence_id: Optional[Union[str, UUID]]
    sequence: Union[Seq, str]
    raw_sequence: Optional[Union[Seq, str]] = None
    sequnce_aa: Optional[Union[Seq, str]] = None
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="input")

    @validator("sequence_id", always=True)
    @classmethod
    def validate_sequence_id(cls, v: Optional[Union[str, UUID]]) -> str:
        """If no sequence id is provided, use a unique UUID"""
        if not v:
            v = uuid4()
        return str(v)

    @validator("sequence", "raw_sequece")
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

    @validator("sequnce_aa")
    @classmethod
    def validate_sequence_aa(cls, v: Union[Seq, str], values: Dict[Any, Any]) -> Union[None, Seq]:
        aa_validator: re.Pattern[Any] = get_aa_validator_regex()
        if isinstance(v, Seq):
            v = str(v)
        if not aa_validator.match(v):
            raise ValueError(f"{v} is not a valid amino acid sequence, must only contain ACDEFGHIKLMNPQRSTVWXY")
        return Seq(v)

    class Config:
        arbitrary_types_allowed = True


class PrimaryAnnotations(BaseModel):
    """
    The primary outputs of the V(D)J assignment process, which includes the gene locus, V, D, J, and C gene calls, various flags, V(D)J junction sequence, copy number (duplicate_count), and the number of reads contributing to a consensus input sequence (consensus_count). Taken from https://docs.airr-community.org/en/stable/datarep/rearrangements.html

    rev_comp  - True if the alignment is on the opposite strand (reverse complemented) with respect to the query sequence. If True then all output data, such as alignment coordinates and sequences, are based on the reverse complement of ‘sequence’.

    productive : True if the V(D)J sequence is predicted to be productive.

    vj_in_frame : True if the V and J gene alignments are in-frame.

    stop_codon :  True if the aligned sequence contains a stop codon.

    complete_vdj : True if the sequence alignment spans the entire V(D)J region. Meaning, sequence_alignment includes both the first V gene codon that encodes the mature polypeptide chain (i.e., after the leader sequence) and the last complete codon of the J gene (i.e., before the J:C splice site). This does not require an absence of deletions within the internal FWR and CDR regions of the alignment.

    locus :  Gene locus (chain type). Note that this field uses a controlled vocabulary that is meant to provide a generic classification of the locus, not necessarily the correct designation according to a specific nomenclature.

    v_call : V gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHV4-59*01 if using IMGT/GENE-DB).

    d_call :  First or only D gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHD3-10*01 if using IMGT/GENE-DB).

    d2_call : Second D gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHD3-10*01 if using IMGT/GENE-DB).

    j_call : J gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHJ4*02 if using IMGT/GENE-DB).

    c_call : Constant region gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHG1*01 if using IMGT/GENE-DB).
    """

    # Airr fields
    rev_comp: bool
    productive: bool
    vj_in_frame: Optional[bool]
    stop_codon: Optional[bool]
    complete_vdj: Optional[bool]
    locus: Optional[str]
    v_call: Union[str, List[str]]
    d_call: Optional[str]
    d2_call: Optional[str]
    j_call: Union[str, List[str]]
    c_call: Optional[str]

    # Non Airr fields
    v_call_top: Optional[str]  # Non Airr
    v_call_top_gene: Optional[str]  # Non Airr
    v_call_top_allele: Optional[str]  # Non Airr
    d_call_gene: Optional[str]  # Non Airr
    d_call_allele: Optional[str]  # Non Airr
    j_call_top: Optional[str]  # Non Airr
    j_call_top_gene: Optional[str]  # Non Airr
    j_call_top_allele: Optional[str]  # Non Airr
    c_call_allele: Optional[str]  # Non Airr
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="primary_annotations")

    class Config:
        arbitrary_types_allowed = True


class AlignmentAnnotations(BaseModel):
    """
    Detailed alignment annotations including the input and germline sequences used in the alignment; score, identity, statistical support (E-value, likelihood, etc); and the alignment itself through CIGAR strings for each aligned gene.

    Attributes:
    -----------

    sequence_alignment: Aligned portion of query sequence, including any indel corrections or numbering spacers, such as IMGT-gaps. Typically, this will include only the V(D)J region, but that is not a requirement.

    sequence_alignment_aa :Amino acid translation of the aligned query sequence.

    germline_alignment: Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and spacers (if any).

    germline_alignment_aa : Amino acid translation of the assembled germline sequence.

    v_score: Alignment score for the V gene.

    v_identity: Fractional identity for the V gene alignment.

    v_support: V gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the V gene assignment as defined by the alignment tool.

    v_cigar :CIGAR string for the V gene alignment.

    d_score: Alignment score for the first or only D gene alignment.

    d_identity: Fractional identity for the first or only D gene alignment.

    d_support : D gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the first or only D gene as defined by the alignment tool.

    d_cigar: CIGAR string for the first or only D gene alignment.

    d2_score: Alignment score for the second D gene alignment.

    d2_identity: number: Fractional identity for the second D gene alignment.

    d2_support: D gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the second D gene as defined by the alignment tool.

    d2_cigar: CIGAR string for the second D gene alignment.

    j_score:Alignment score for the J gene alignment.

    j_identity: Fractional identity for the J gene alignment.

    j_support: J gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the J gene assignment as defined by the alignment tool.

    j_cigar: CIGAR string for the J gene alignment.

    c_score: Alignment score for the C gene alignment.

    c_identity: Fractional identity for the C gene alignment.

    c_support: C gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the C gene assignment as defined by the alignment tool.

    c_cigar: CIGAR string for the C gene alignment.
    """

    sequecne_alignment: Union[str, Seq]
    seqeunce_alignment_aa: Optional[Union[str, Seq]]
    germline_alignement: Union[str, Seq]
    germline_alignment_aa: Optional[Union[str, Seq]]
    v_score: Optional[float]
    v_identity: Optional[float]
    v_support: Optional[float]
    v_cigar: str
    d_score: Optional[float]
    d_identity: Optional[float]
    d_support: Optional[float]
    d_cigar: str
    d2_score: Optional[float]
    d2_identity: Optional[float]
    d2_support: Optional[float]
    d2_cigar: Optional[str]
    j_score: Optional[float]
    j_identity: Optional[float]
    j_support: Optional[float]
    j_cigar: str
    junction: Union[str, Seq]
    junction_aa: Optional[Union[str, Seq]]
    np1: Optional[Union[str, Seq]]
    np1_aa: Optional[Union[str, Seq]]
    np2: Optional[Union[str, Seq]]
    np2_aa: Optional[Union[str, Seq]]
    np3: Optional[Union[str, Seq]]
    np3_aa: Optional[Union[str, Seq]]
    c_score: Optional[float]
    c_identity: Optional[float]
    c_support: Optional[float]
    c_cigar: Optional[str]
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="alignment_annotations")

    class Config:
        arbitrary_types_allowed = True


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

    v_sequence_start: Optional[int] = None
    v_sequence_end: Optional[int] = None
    v_germline_start: Optional[int] = None
    v_germline_end: Optional[int] = None
    v_alignment_start: Optional[int] = None
    v_alignment_end: Optional[int] = None
    d_sequence_start: Optional[int] = None
    d_sequence_end: Optional[int] = None
    d_germline_start: Optional[int] = None
    d_germline_end: Optional[int] = None
    d_alignment_start: Optional[int] = None
    d_alignment_end: Optional[int] = None
    d2_sequence_start: Optional[int] = None
    d2_sequence_end: Optional[int] = None
    d2_germline_start: Optional[int] = None
    d2_germline_end: Optional[int] = None
    d2_alignment_start: Optional[int] = None
    d2_alignment_end: Optional[int] = None
    j_sequence_start: Optional[int] = None
    j_sequence_end: Optional[int] = None
    j_germline_start: Optional[int] = None
    j_germline_end: Optional[int] = None
    j_alignment_start: Optional[int] = None
    j_alignment_end: Optional[int] = None
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="alignment_positions")

    class Config:
        arbitrary_types_allowed = True


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
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="region_sequences")

    class Config:
        arbitrary_types_allowed = True


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


class JunctionLenghts(BaseModel):
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


class ConstantChain(BaseModel):
    constant_chain: Union[Seq, str]
    isotype: str
    species: str
    category: RearrargmentCategory

    class Config:
        arbitrary_types_allowed = True


class Chain(BaseModel):
    vdj_seq: VDJSequence
    constant_seq: Optional[ConstantChain]


class Chains(BaseModel):
    chains: List[Chain]


class Antibody(BaseModel):
    heavy_chain: Chain
    light_chain: Chain


class Antibodies(BaseModel):
    antibodies: List[Antibody]


class TCR(BaseModel):
    heavy_chain: Chain
    light_chain: Chain
