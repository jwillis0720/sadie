from collections import UserString
from typing import Callable, Generator

from pydantic.fields import ModelField

CHAINS = ["L", "H", "K", "A", "B", "G", "D"]


class Chain(UserString):

    chains = CHAINS

    @classmethod
    def __get_validators__(cls) -> Generator[Callable[[str, ModelField], str], None, None]:
        yield cls.validate

    @classmethod
    def validate(cls, value: str, field: ModelField) -> str:
        if not isinstance(value, str):
            raise ValueError(f"{field} [{value}] must be a string")
        value = value.strip().upper()
        if value not in CHAINS:
            # print('CHAINS', value)
            raise ValueError(f"{field} [{value}] must be a string")

        return value
