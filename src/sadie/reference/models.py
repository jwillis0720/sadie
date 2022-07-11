from pydantic import BaseModel, validator
from typing import Dict, List


class Species(BaseModel):
    """What species to retrieve"""

    species: str


class Source(BaseModel):
    """What database to retrieve"""

    source: str


class GeneEntry(BaseModel):
    """V,D or J Gene Entry with validation"""

    species: str
    gene: str
    source: str

    # values: a dict containing the name-to-value mapping of any previously-validated fields
    @validator("species")
    def check_species(cls, v: str) -> str:
        # pylint: disable=no-self-argument
        return Species(**{"species": v}).species

    @validator("gene")
    def check_vgene(cls, v: str, values: Dict[str, str]) -> str:
        # pylint: disable=no-self-argument
        if v[3] not in ["V", "D", "J"]:
            raise ValueError(f"gene must contain V,D or J at 3rd index, current have {v[3]} in {v} ")
        return v

    @validator("source")
    def check_source(cls, v: str) -> str:
        # pylint: disable=no-self-argument
        if v not in ["imgt", "custom"]:
            raise ValueError(f"{v} is not a valid source, chocies are 'imgt' or 'custom'")
        return v


class GeneEntries(BaseModel):
    """V,D or J Gene Entry with validation"""

    species: str
    genes: List[str]
    source: str

    @validator("species")
    def check_species(cls, v: str) -> str:
        # pylint: disable=no-self-argument
        return Species(**{"species": v}).species

    @validator("genes", each_item=True)
    def check_vgene(cls, v: str, values: Dict[str, str]) -> str:
        # pylint: disable=no-self-argument
        if v[3] not in ["V", "D", "J", "C", "A", "G", "M", "E"]:
            raise ValueError(f"gene must contain V,D,J or C at 3rd index, current have {v[3]} in {v} ")
        return v

    @validator("source")
    def check_source(cls, v: str) -> str:
        # pylint: disable=no-self-argument
        if v not in ["imgt", "custom"]:
            raise ValueError(f"{v} is not a valid source, chocies are 'imgt' or 'custom'")
        return v
