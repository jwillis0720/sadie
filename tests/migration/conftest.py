"""Session-scoped fixtures for G3-Germlines parity testing."""
import pytest
from pathlib import Path

from sadie.reference import References


# Use human-only reference for parity testing (avoids multi-species chimera issues)
YAML_PATH = Path(__file__).parent / "reference_parity_test.yml"


@pytest.fixture(scope="session")
def g3_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference.g3.yml using G3 backend.
    
    Returns:
        Path to the built database directory.
    """
    refs = References.from_yaml(YAML_PATH, use_germlines=False)
    outpath = tmp_path_factory.mktemp("g3_db")
    return refs.make_airr_database(outpath)


@pytest.fixture(scope="session")
def germlines_database(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Build database from reference.g3.yml using Germlines backend.
    
    Returns:
        Path to the built database directory.
    """
    refs = References.from_yaml(YAML_PATH, use_germlines=True)
    outpath = tmp_path_factory.mktemp("germlines_db")
    return refs.make_airr_database(outpath)
