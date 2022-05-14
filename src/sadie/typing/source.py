from collections import UserString
from typing import Callable, Generator

from pydantic.fields import ModelField

SOURCES = ["imgt", "custom"]


class Source(UserString):

    sources = SOURCES

    @classmethod
    def __get_validators__(cls) -> Generator[Callable[[str, ModelField], str], None, None]:
        yield cls.validate

    @classmethod
    def validate(cls, value: str, field: ModelField) -> str:
        if not isinstance(value, str):
            raise ValueError(f"{field} [{value}] must be a string")
        value = value.strip().lower()
        if value not in SOURCES:
            raise ValueError(f"{field} [{value}] must be in {SOURCES}")

        return value
