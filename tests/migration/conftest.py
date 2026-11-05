"""Session-scoped fixtures for G3-Germlines parity testing."""
from pathlib import Path

import pytest

from sadie.reference import References

# Use human-only reference for parity testing (avoids multi-species chimera issues)
YAML_PATH = Path(__file__).parent / "reference_parity_test.yml"
# Mixed-source reference: same alleles but 5 from OGRDB + 5 from VDJbase
YAML_PATH_MIX = Path(__file__).parent / "reference_parity_test_mix.yml"


@pytest.fixture(scope="session")
def germlines_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference_parity_test.yml using Germlines backend.

    Returns:
        Path to the built database directory.
    """
    refs = References.from_yaml(YAML_PATH, use_germlines=True)
    outpath = tmp_path_factory.mktemp("germlines_db")
    return refs.make_airr_database(outpath)


@pytest.fixture(scope="session")
def mixed_source_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference_parity_test_mix.yml using Germlines backend.

    This database has the same alleles as germlines_database but pulls:
    - 5 alleles from OGRDB: IGHV1-18*01, IGHV1-2*02, IGHD1-1*01, IGHJ1*01, IGHJ3*02
    - 5 alleles from VDJbase: IGHV3-30*01, IGHV3-21*01, IGHD2-2*01, IGHJ2*01, IGHJ4*02
    - Remaining alleles from IMGT

    All selected cross-source alleles have IDENTICAL sequences across providers.

    Returns:
        Path to the built database directory.
    """
    refs = References.from_yaml(YAML_PATH_MIX, use_germlines=True)
    outpath = tmp_path_factory.mktemp("mixed_db")
    return refs.make_airr_database(outpath)
