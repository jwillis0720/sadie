"""
This file contains objects to represent the rearrangemnt scheme
https://docs.airr-community.org/en/stable/datarep/rearrangements.html
"""
import re
import warnings
from functools import lru_cache
from typing import Any, Dict, List, Optional, Set, Union
from uuid import UUID, uuid4

from Bio.Seq import Seq
from pydantic import BaseModel, validator


@lru_cache(maxsize=1)
def get_nt_validator_regex() -> re.Pattern[Any]:
    return re.compile(r"^[ACNTGacntg]+\Z")


@lru_cache(maxsize=1)
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

    sequence_id - Unique query sequence identifier for the Rearrangment. Most often this will be the input sequence header or a substring thereof, but may also be a custom identifier defined by the tool in cases where query sequences have been combined in some fashion prior to alignment. When downloaded from an AIRR Data Commons repository, this will usually be a universally unique record locator for linking with other objects in the AIRR Data Model.

    sequence - The query nucleotide sequence. Usually, this is the unmodified input sequence, which may be reverse complemented if necessary. In some cases, this field may contain consensus sequences or other types of collapsed input sequences if these steps are performed prior to alignment.
    """

    sequence_id: Optional[Union[str, UUID]]
    sequence: Union[Seq, str]
    sequnce_aa: Optional[Union[Seq, str]]
    species: str
    category: Optional[RearrargmentCategory] = RearrargmentCategory(category="input")

    @validator("sequence_id", always=True)
    @classmethod
    def validate_sequence_id(cls, v: Optional[Union[str, UUID]]):
        """If no sequence id is provided, use a unique UUID"""
        if not v:
            v = uuid4()
        return str(v)

    @validator("sequence")
    @classmethod
    def validate_sequence(cls, v: Union[Seq, str]) -> Seq:
        """Check for valid nt sequence"""
        if isinstance(v, Seq):
            v = str(v)
        nt_validator: re.Pattern[Any] = get_nt_validator_regex()
        if not nt_validator.match(v):
            raise ValueError(f"{v} is not a valid nucleotide sequence, must only contain ACGTN")
        else:
            return Seq(v)

    @validator("sequnce_aa", always=True)
    @classmethod
    def validate_sequence_aa(cls, v: Optional[Union[Seq, str]], values: Dict[Any, Any]) -> Seq:
        sequence_from_nt: Seq = Seq(values["sequence"])
        if not v:
            """If no sequence aa is provided, use the nt sequence to translate"""
            v = sequence_from_nt.translate().__str__()
        else:
            if v.__str__() != sequence_from_nt.translate().__str__():
                warnings.warn(
                    f"{v.__str__()} is not the same as the translated seq: {sequence_from_nt.translate().__str__()}",
                    UserWarning,
                )
        aa_validator: re.Pattern[Any] = get_aa_validator_regex()
        if not aa_validator.match(v):
            raise ValueError(f"{v} is not a valid amino acid sequence, must only contain ACDEFGHIKLMNPQRSTVWXY")
        return Seq(v)

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
