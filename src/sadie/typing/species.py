from collections import UserString
from typing import Any

from pydantic import field_validator
from pydantic_core import core_schema

# TODO: go through and see which are viable to use; tests need to be fixed first in test_g3 to handle this
SPECIES = {
    "rhesus": "macaque",
    "homo_sapiens": "human",
    "mus": "mouse",
    "rattus_norvegicus": "rat",
    "oryctolagus_cuniculus": "rabbit",
    "macaca_mulatta": "rhesus",
    "sus_scrofa": "pig",
    "vicugna_pacos": "alpaca",
    "bos_taurus": "cow",
    "alpaca": "alpaca",
    "human": "human",
    "macaque": "macaque",
    "mouse": "mouse",
    "rabbit": "rabbit",
    "dog": "dog",
    "cat": "cat",
    "rat": "rat",
    "pig": "pig",
    # 'amberjack': 'amberjack',
    # 'bass': 'bass',
    # 'boar': 'boar',
    # 'bull_shark': 'bull_shark',
    # 'camel': 'camel',
    # 'carp': 'carp',
    # 'catfish': 'catfish',
    # 'char': 'char',
    # 'chinese_perch': 'chinese_perch',
    # 'clearnose_skate': 'clearnose_skate',
    # 'cod': 'cod',
    # 'crab_eating_macaque': 'crab_eating_macaque',
    # 'dolphin': 'dolphin',
    # 'ferret': 'ferret',
    # 'flounder': 'flounder',
    # 'goat': 'goat',
    # 'goldfish': 'goldfish',
    # 'horn_shark': 'horn_shark',
    # 'horse': 'horse',
    # 'icefish': 'icefish',
    # 'junglefowl': 'junglefowl',
    # 'ladyfish': 'ladyfish',
    # 'little_skate': 'little_skate',
    # 'night_monkey': 'night_monkey',
    # 'nurse_shark': 'nurse_shark',
    # 'platypus': 'platypus',
    # 'pufferfish': 'pufferfish',
    # 'ratfish': 'ratfish',
    # 'rockcod': 'rockcod',
    # 'salmon': 'salmon',
    # 'sandbar_shark': 'sandbar_shark',
    # 'shark': 'shark',
    # 'sheep': 'sheep',
    # 'spotted_wolffish': 'spotted_wolffish',
    # 'trout': 'trout',
    # 'tubot': 'tubot',
    # 'wobbegong': 'wobbegong',
    # 'zebrafish': 'zebrafish',
}


class Species(UserString):
    species = SPECIES

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
        value = value.strip().lower().replace(" ", "_")
        if value not in SPECIES:
            raise ValueError(f"value [{value}] must be in {SPECIES.keys()}")
        value = SPECIES[value]
        return value
