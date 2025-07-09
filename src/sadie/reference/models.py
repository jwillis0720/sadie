from typing import Dict, List

from pydantic import BaseModel, field_validator


class Species(BaseModel):
    """What species to retrieve"""

    species: str


class Source(BaseModel):
    """What source to retrieve"""

    source: str


class GeneEntry(BaseModel):
    """V,D or J Gene Entry with validation"""

    species: str
    gene: str
    source: str

    # values: a dict containing the name-to-value mapping of any previously-validated fields
    @field_validator("species")
    @classmethod
    def check_species(cls, v: str) -> str:
        return Species(**{"species": v}).species

    @field_validator("gene")
    @classmethod
    def check_vgene(cls, v: str) -> str:
        if v[3] not in ["V", "D", "J"]:
            raise ValueError(f"gene must contain V,D or J at 3rd index, current have {v[3]} in {v} ")
        return v

    @field_validator("source")
    @classmethod
    def check_source(cls, v: str) -> str:
        if v not in ["imgt", "custom"]:
            raise ValueError(f"{v} is not a valid source, chocies are 'imgt' or 'custom'")
        return v


class GeneEntries(BaseModel):
    """V,D or J Gene Entry with validation"""

    species: str
    genes: List[str]
    source: str

    @field_validator("species")
    @classmethod
    def check_species(cls, v: str) -> str:
        return Species(**{"species": v}).species

    @field_validator("genes")
    @classmethod
    def check_vgene(cls, v: List[str]) -> List[str]:
        for gene in v:
            if gene[3] not in ["V", "D", "J", "C", "A", "G", "M", "E"]:
                raise ValueError(f"gene must contain V,D,J or C at 3rd index, current have {gene[3]} in {gene} ")
        return v

    @field_validator("source")
    @classmethod
    def check_source(cls, v: str) -> str:
        if v not in ["imgt", "custom"]:
            raise ValueError(f"{v} is not a valid source, chocies are 'imgt' or 'custom'")
        return v
