from collections import UserString
from typing import Any

from pydantic_core import core_schema

SOURCES = ["imgt", "custom"]


class Source(UserString):
    sources = SOURCES

    @classmethod
    def __get_pydantic_core_schema__(cls, source_type: Any, handler: Any) -> core_schema.CoreSchema:
        return core_schema.no_info_after_validator_function(
            cls.validate,
            core_schema.str_schema(),
        )

    @classmethod
    def validate(cls, value: str) -> str:
        if not isinstance(value, str):
            raise ValueError(f"value [{value}] must be a string")
        value = value.strip().lower()
        if value not in SOURCES:
            raise ValueError(f"value [{value}] must be in {SOURCES}")

        return value
