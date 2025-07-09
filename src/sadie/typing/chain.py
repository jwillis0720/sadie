from collections import UserString
from typing import Any

from pydantic_core import core_schema

CHAINS = ["L", "H", "K", "A", "B", "G", "D"]


class Chain(UserString):
    chains = CHAINS

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
        value = value.strip().upper()
        if value not in CHAINS:
            # print('CHAINS', value)
            raise ValueError(f"value [{value}] must be a string")

        return value
