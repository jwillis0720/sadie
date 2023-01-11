from __future__ import annotations

from contextlib import AbstractContextManager
from contextlib import nullcontext as does_not_raise

import pytest
from pydantic import validate_arguments

from sadie.typing.chain import Chain
from sadie.typing.source import Source
from sadie.typing.species import SPECIES, Species


class TestTyping:
    @pytest.mark.parametrize(
        "source, expectation",
        [
            ("imgt", does_not_raise()),
            ("custom", does_not_raise()),
            ("", pytest.raises(ValueError)),
            ("test", pytest.raises(ValueError)),
            (1, pytest.raises(ValueError)),
        ],
    )
    def test_source(self, source: str, expectation: AbstractContextManager | pytest.raises) -> None:  # type: ignore
        @validate_arguments
        def dummy(source: Source) -> None:
            pass

        with expectation:
            dummy(source)

    @pytest.mark.parametrize(
        "chain, expectation",
        [(chain, does_not_raise()) for chain in ["L", "H", "K", "A", "B", "G", "D"]]  # type: ignore
        + [
            ("", pytest.raises(ValueError)),
            ("test", pytest.raises(ValueError)),
            (1, pytest.raises(ValueError)),
        ],
    )
    def test_chain(self, chain: str, expectation: AbstractContextManager | pytest.raises) -> None:  # type: ignore
        @validate_arguments
        def dummy(chain: Chain) -> None:
            pass

        with expectation:
            dummy(chain)

    @pytest.mark.parametrize(
        "species, expectation",
        [(species, does_not_raise()) for species in SPECIES]  # type: ignore
        + [
            ("", pytest.raises(ValueError)),
            ("test", pytest.raises(ValueError)),
            (1, pytest.raises(ValueError)),
        ],
    )
    def test_species(self, species: str, expectation: AbstractContextManager | pytest.raises) -> None:  # type: ignore
        @validate_arguments
        def dummy(species: Species) -> None:
            pass

        with expectation:
            dummy(species)
