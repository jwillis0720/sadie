from __future__ import annotations

from contextlib import AbstractContextManager
from contextlib import nullcontext as does_not_raise

import pytest
from pydantic import BaseModel, ValidationError

from sadie.typing.chain import Chain
from sadie.typing.source import Source
from sadie.typing.species import SPECIES, Species


class TestModels(BaseModel):
    source: Source
    species: Species
    chain: Chain


class TestSpeciesValidation:
    def test_species_validation_success(self):
        model = TestModels(source="imgt", species="human", chain="H")
        assert model.species == "human"
        assert model.source == "imgt"
        assert model.chain == "H"

    def test_species_validation_transformation(self):
        model = TestModels(source="imgt", species="homo_sapiens", chain="H")
        assert model.species == "human"

    def test_species_validation_failure(self):
        with pytest.raises(ValidationError):
            TestModels(source="imgt", species="invalid_species", chain="H")


class TestChainValidation:
    def test_chain_validation_success(self):
        model = TestModels(source="imgt", species="human", chain="H")
        assert model.chain == "H"

    def test_chain_validation_lowercase(self):
        model = TestModels(source="imgt", species="human", chain="h")
        assert model.chain == "H"

    def test_chain_validation_failure(self):
        with pytest.raises(ValidationError):
            TestModels(source="imgt", species="human", chain="X")


class TestSourceValidation:
    def test_source_validation_success(self):
        model = TestModels(source="imgt", species="human", chain="H")
        assert model.source == "imgt"

    def test_source_validation_custom(self):
        model = TestModels(source="custom", species="human", chain="H")
        assert model.source == "custom"

    def test_source_validation_failure(self):
        with pytest.raises(ValidationError):
            TestModels(source="invalid_source", species="human", chain="H")
