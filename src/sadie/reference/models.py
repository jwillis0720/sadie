from pydantic import BaseModel, validator
from typing import Dict, Optional, List


class Species(BaseModel):
    species: str


class GeneEntry(BaseModel):
    """V,D or J Gene Entry with validation"""

    species: str
    sub_species: Optional[str] = ""
    gene: str
    database: str

    # values: a dict containing the name-to-value mapping of any previously-validated fields
    @validator("species")
    def check_species(cls, v: str) -> str:
        # pylint: disable=no-self-argument
        return Species(**{"species": v}).species

    @validator("sub_species", always=True)
    def check_sub_species(cls, v: str, values: Dict[str, str]) -> str:
        # pylint: disable=no-self-argument
        if not v:
            return values["species"]
        return Species(**{"species": v}).species

    @validator("gene")
    def check_vgene(cls, v: str, values: Dict[str, str]) -> str:
        # pylint: disable=no-self-argument
        if v[3] not in ["V", "D", "J"]:
            raise ValueError(f"gene must contain V,D or J at 3rd index, current have {v[3]} in {v} ")
        return v

    @validator("database")
    def check_database(cls, v: str) -> str:
        # pylint: disable=no-self-argument
        if v not in ["imgt", "custom"]:
            raise ValueError(f"{v} is not a valid database, chocies are 'imgt' or 'custom'")
        return v


class GeneEntries(BaseModel):
    """V,D or J Gene Entry with validation"""

    species: str
    sub_species: Optional[str] = ""
    gene: List[str]
    database: str

    @validator("species")
    def check_species(cls, v: str) -> str:
        # pylint: disable=no-self-argument
        return Species(**{"species": v}).species

    @validator("sub_species", always=True)
    def check_sub_species(cls, v: str, values: Dict[str, str]) -> str:
        # pylint: disable=no-self-argument
        if not v:
            return values["species"]
        return Species(**{"species": v}).species

    @validator("gene", each_item=True)
    def check_vgene(cls, v: str, values: Dict[str, str]) -> str:
        # pylint: disable=no-self-argument
        if v[3] not in ["V", "D", "J"]:
            raise ValueError(f"gene must contain V,D or J at 3rd index, current have {v[3]} in {v} ")
        return v

    @validator("database")
    def check_database(cls, v: str) -> str:
        # pylint: disable=no-self-argument
        if v not in ["imgt", "custom"]:
            raise ValueError(f"{v} is not a valid database, chocies are 'imgt' or 'custom'")
        return v
