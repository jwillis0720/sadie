from typing import Any, Dict, Optional

import pandas as pd
from pydantic import BaseModel, Field, root_validator


# TODO: https://docs.airr-community.org/en/stable/news.html
# confirm match 1.4.1 schema
class AirrSeriesModel(BaseModel):
    sequence_id: Optional[str] = Field(
        description="Unique query sequence identifier for the Rearrangement.",
    )  # Should be required
    c_alignment_end: Optional[int] = Field(
        description="End position of the C gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).",
    )  # New in Version 1.4.1
    c_alignment_start: Optional[int] = Field(
        description="Start position of the C gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).",
    )  # New in Version 1.4.1
    c_call: Optional[str] = Field(description="Constant region gene with allele.")  # New in Version 1.4.1
    c_germline_end: Optional[int] = Field(
        description="Alignment end position in the C gene reference sequence (1-based closed interval).",
    )  # New in Version 1.4.1
    c_germline_start: Optional[int] = Field(
        description="Alignment start position in the C gene reference sequence (1-based closed interval).",
    )  # New in Version 1.4.1
    c_sequence_end: Optional[int] = Field(
        description="End position of the C gene in the query sequence (1-based closed interval).",
    )  # New in Version 1.4.1
    c_sequence_start: Optional[int] = Field(
        description="Start position of the C gene in the query sequence (1-based closed interval).",
    )  # New in Version 1.4.1
    cdr1: Optional[str] = Field(description="Nucleotide sequence of the aligned CDR1 region.")
    cdr1_aa: Optional[str] = Field(description="Amino acid translation of the cdr1 field.")
    cdr1_end: Optional[int] = Field(description="CDR1 end position in the query sequence (1-based closed interval).")
    cdr1_start: Optional[int] = Field(
        description="CDR1 start position in the query sequence (1-based closed interval).",
    )
    cdr2: Optional[str] = Field(description="Nucleotide sequence of the aligned CDR2 region.")
    cdr2_aa: Optional[str] = Field(description="Amino acid translation of the cdr2 field.")
    cdr2_end: Optional[int] = Field(description="CDR2 end position in the query sequence (1-based closed interval).")
    cdr2_start: Optional[int] = Field(
        description="CDR2 start position in the query sequence (1-based closed interval).",
    )
    cdr3: Optional[str] = Field(description="Nucleotide sequence of the aligned CDR3 region.")
    cdr3_aa: Optional[str] = Field(description="Amino acid translation of the cdr3 field.")
    # cdr3_aa_length: Optional[int] = Field(description="CDR1 end position in the query sequence (1-based closed interval).")
    cdr3_end: Optional[int] = Field(description="CDR13 end position in the query sequence (1-based closed interval).")
    cdr3_start: Optional[int] = Field(
        description="CDR3 start position in the query sequence (1-based closed interval).",
    )
    complete_vdj: Optional[bool] = Field(
        description="True if the sequence alignment spans the entire V(D)J region. Can contain deletions within the internal FWR and CDR regions of the alignment.",
    )
    d_alignment_end: Optional[int] = Field(
        description="End position of the first or only D gene in both the sequence_alignment and germline_alignment fields (1-based closed interval).",
    )
    d_alignment_start: Optional[int] = Field(
        description="Start position of the first or only D gene in both the sequence_alignment and germline_alignment fields (1-based closed interval).",
    )
    d_call: Optional[str] = Field(
        description="First or only D gene with allele. Use relevant gene/allele nomenclature should be followed (e.g., IGHD3-10*01 if using IMGT/GENE-DB) for known reference DB.",
    )
    d_cigar: Optional[str] = Field(description="CIGAR string for the first or only D gene alignment.")
    # d_family: Optional[str] = Field(description="") CK: Can't find this on https://docs.airr-community.org/en/stable/datarep/rearrangements.html
    d_frame: Optional[bool] = Field(
        description="Numerical reading frame (1, 2, 3) of the first or only D gene in the query nucleotide sequence, where frame 1 is relative to the first codon of D gene reference sequence.",
    )  # Added in 1.4.1
    d_germline_alignment: Optional[str] = Field(
        description="Aligned D gene germline sequence spanning the same region as the d_sequence_alignment field and including the same set of corrections and spacers (if any).",
    )
    d_germline_alignment_aa: Optional[str] = Field(
        description="Amino acid translation of the d_germline_alignment field.",
    )
    d_germline_end: Optional[int] = Field(
        description="Alignment end position in the D gene reference sequence for the first or only D gene (1-based closed interval).",
    )
    d_germline_start: Optional[int] = Field(
        description="Alignment start position in the D gene reference sequence for the first or only D gene (1-based closed interval).",
    )
    d_identity: Optional[float] = Field(description="Fractional identity for the first or only D gene alignment.")
    d_score: Optional[float] = Field(description="Alignment score for the first or only D gene alignment.")
    d_sequence_alignment: Optional[str] = Field(
        description="Aligned portion of query sequence assigned to the first or only D gene, including any indel corrections or numbering spacers.",
    )
    d_sequence_alignment_aa: Optional[str] = Field(
        description="Amino acid translation of the d_sequence_alignment field.",
    )
    d_sequence_end: Optional[int] = Field(
        description="End position of the first or only D gene in the query sequence. (1-based closed interval).",
    )
    d_sequence_start: Optional[int] = Field(
        description="Start position of the first or only D gene in the query sequence. (1-based closed interval).",
    )
    d_support: Optional[str] = Field(
        description="D gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the first or only D gene as defined by the alignment tool.",
    )
    fwr1: Optional[str] = Field(description="Nucleotide sequence of the aligned FWR1 region.")
    fwr1_aa: Optional[str] = Field(description="Amino acid translation of the fwr1 field.")
    fwr1_end: Optional[int] = Field(description="FWR1 end position in the query sequence (1-based closed interval).")
    fwr1_start: Optional[int] = Field(
        description="FWR1 start position in the query sequence (1-based closed interval)."
    )
    fwr2: Optional[str] = Field(description="Nucleotide sequence of the aligned FWR2 region.")
    fwr2_aa: Optional[str] = Field(description="Amino acid translation of the fwr2 field.")
    fwr2_end: Optional[int] = Field(description="FWR2 end position in the query sequence (1-based closed interval).")
    fwr2_start: Optional[int] = Field(
        description="FWR2 start position in the query sequence (1-based closed interval).",
    )
    fwr3: Optional[str] = Field(description="Nucleotide sequence of the aligned FWR3 region.")
    fwr3_aa: Optional[str] = Field(description="Amino acid translation of the fwr3 field.")
    fwr3_end: Optional[int] = Field(description="FWR3 end position in the query sequence (1-based closed interval).")
    fwr3_start: Optional[int] = Field(
        description="FWR3 start position in the query sequence (1-based closed interval).",
    )
    fwr4: Optional[str] = Field(description="Nucleotide sequence of the aligned FWR4 region.")
    fwr4_aa: Optional[str] = Field(description="Amino acid translation of the fwr4 field.")
    fwr4_end: Optional[int] = Field(description="FWR4 end position in the query sequence (1-based closed interval).")
    fwr4_start: Optional[int] = Field(
        description="FWR4 start position in the query sequence (1-based closed interval).",
    )
    germline_alignment: Optional[str] = Field(
        description="Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and spacers (if any).",
    )
    germline_alignment_aa: Optional[str] = Field(
        description="Amino acid translation of the assembled germline sequence.",
    )
    j_alignment_end: Optional[int] = Field(
        description="End position of the J gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).",
    )
    j_alignment_start: Optional[int] = Field(
        description="Start position of the J gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).",
    )
    j_call: Optional[str] = Field(description="J gene with allele.")
    j_cigar: Optional[str] = Field(description="CIGAR string for the J gene alignment.")
    # j_family: Optional[str] = Field(description=" ") #CK Not found
    j_frameshift: Optional[str] = Field(
        description="True if the J gene in the query nucleotide sequence contains a translational frameshift relative to the frame of the J gene reference sequence.",
    )  # added in 1.4
    j_germline_alignment: Optional[str] = Field(
        description="Aligned J gene germline sequence spanning the same region as the j_sequence_alignment field and including the same set of corrections and spacers (if any).",
    )
    j_germline_alignment_aa: Optional[str] = Field(
        description="Amino acid translation of the j_germline_alignment field.",
    )
    j_germline_end: Optional[int] = Field(
        description="Alignment end position in the J gene reference sequence (1-based closed interval).",
    )
    j_germline_start: Optional[int] = Field(
        description="Alignment start position in the J gene reference sequence (1-based closed interval)."
    )
    j_identity: Optional[float] = Field(description="Fractional identity for the J gene alignment.")
    j_score: Optional[float] = Field(description="Alignment score for the J gene alignment.")
    j_sequence_alignment: Optional[str] = Field(
        description="Aligned portion of query sequence assigned to the J gene, including any indel corrections or numbering spacers.",
    )
    j_sequence_alignment_aa: Optional[str] = Field(
        description="Amino acid translation of the j_sequence_alignment field.",
    )
    j_sequence_end: Optional[int] = Field(
        description="End position of the J gene in the query sequence (1-based closed interval).",
    )
    j_sequence_start: Optional[int] = Field(
        description="Start position of the J gene in the query sequence (1-based closed interval).",
    )
    j_support: Optional[str] = Field(
        description="J gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the J gene assignment as defined by the alignment tool.",
    )
    junction: Optional[str] = Field(
        description="Junction region nucleotide sequence, where the junction is defined as the CDR3 plus the two flanking conserved codons.",
    )
    junction_aa: Optional[str] = Field(description="Amino acid translation of the junction.")
    junction_aa_length: Optional[int] = Field(description="Number of amino acids in the junction sequence.")
    junction_length: Optional[int] = Field(description="Number of nucleotides in the junction sequence.")
    locus: Optional[str] = Field(
        description="Gene locus (chain type). Note that this field uses a controlled vocabulary that is meant to provide a generic classification of the locus, not necessarily the correct designation according to a specific nomenclature.",
    )
    np1: Optional[str] = Field(
        description="Nucleotide sequence of the combined N/P region between the V gene and first D gene alignment or between the V gene and J gene alignments.",
    )
    np1_length: Optional[int] = Field(
        description="Number of nucleotides between the V gene and first D gene alignments or between the V gene and J gene alignments.",
    )
    np2: Optional[str] = Field(
        description="Nucleotide sequence of the combined N/P region between either the first D gene and J gene alignments or the first D gene and second D gene alignments.",
    )
    np2_length: Optional[int] = Field(
        description="Number of nucleotides between either the first D gene and J gene alignments or the first D gene and second D gene alignments.",
    )
    productive: Optional[bool] = Field(description="True if the V(D)J sequence is predicted to be productive.")
    quality: Optional[str] = Field(
        description="The Sanger/Phred quality scores for assessment of sequence quality.",
    )  # CK Added in 1.4.1
    quality_alignment: Optional[str] = Field(
        description="Sanger/Phred quality scores for assessment of sequence_alignment quality.",
    )  # CK Added in 1.4.1
    rev_comp: Optional[bool] = Field(
        description="True if the alignment is on the opposite strand (reverse complemented) with respect to the query sequence.",
    )
    sequence: Optional[str] = Field(description="The query nucleotide sequence.")
    sequence_alignment: Optional[str] = Field(
        description="Aligned portion of query sequence, including any indel corrections or numbering spacers, such as IMGT-gaps. Typically, this will include only the V(D)J region, but that is not a requirement.",
    )
    sequence_alignment_aa: Optional[str] = Field(description="Amino acid translation of the aligned query sequence.")
    stop_codon: Optional[bool] = Field(description="True if the aligned sequence contains a stop codon.")
    v_alignment_end: Optional[int] = Field(
        description="End position of the V gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).",
    )
    v_alignment_start: Optional[int] = Field(
        description="Start position of the V gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).",
    )
    v_call: Optional[str] = Field(
        description="V gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHV4-59*01 if using IMGT/GENE-DB).",
    )
    v_cigar: Optional[str] = Field(description="CIGAR string for the V gene alignment.")
    # v_family: Optional[str] = Field(description=" ")
    v_frameshift: Optional[str] = Field(
        description="True if the V gene in the query nucleotide sequence contains a translational frameshift relative to the frame of the V gene reference sequence.",
    )  # added in 1.4.1
    v_germline_alignment: Optional[str] = Field(
        description="Aligned V gene germline sequence spanning the same region as the v_sequence_alignment field and including the same set of corrections and spacers (if any).",
    )
    v_germline_alignment_aa: Optional[str] = Field(
        description="Amino acid translation of the v_germline_alignment field.",
    )
    v_germline_end: Optional[int] = Field(
        description="Alignment end position in the V gene reference sequence (1-based closed interval).",
    )
    v_germline_start: Optional[int] = Field(
        description="Alignment start position in the V gene reference sequence (1-based closed interval).",
    )
    v_identity: Optional[float] = Field(description="Fractional identity for the V gene alignment.")
    v_score: Optional[float] = Field(description="Alignment score for the V gene.")
    v_sequence_alignment: Optional[str] = Field(
        description="Aligned portion of query sequence assigned to the V gene, including any indel corrections or numbering spacers.",
    )
    v_sequence_alignment_aa: Optional[str] = Field(
        description="Amino acid translation of the v_sequence_alignment field.",
    )
    v_sequence_end: Optional[int] = Field(
        description="End position of the V gene in the query sequence (1-based closed interval).",
    )
    v_sequence_start: Optional[int] = Field(
        description="Start position of the V gene in the query sequence (1-based closed interval).",
    )
    v_support: Optional[str] = Field(
        description="V gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the V gene assignment as defined by the alignment tool.",
    )
    vj_in_frame: Optional[bool] = Field(description="True if the V and J gene alignments are in-frame.")
    # chain: Optional[str] = Field(description=" ")
    species: Optional[str] = Field(
        description="Species"
    )  # CK missing on https://docs.airr-community.org/en/stable/datarep/rearrangements.html

    @root_validator(pre=True)
    def fix_dependent_attrs(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        """Fixes dependent attributes that are not the proper type"""
        cleaned_values = {}
        # Remove all null values
        for k, v in values.items():
            # any null values as string by mistake
            if isinstance(v, str):
                if v.lower() in ["<na>", "na", "nan", "", "none"]:
                    v = None
            # nan to None: we have hybrid types in the data so nan types cannot be leveraged & will break some logic
            # elif isinstance(v, pd._libs.missing.NAType):
            #     v = None
            cleaned_values[k] = v
        return cleaned_values
