"""Pytest fixtures for germlines integration tests.

Provides fixtures for testing offline operation and network isolation.
"""

import socket
from unittest.mock import patch

import pytest


@pytest.fixture
def network_disabled(monkeypatch):
    """Fixture to disable network access for testing offline operation.

    This fixture blocks all socket connections to simulate offline mode.
    Tests using this fixture verify that germlines module works without
    network access.

    Usage:
        def test_offline(network_disabled):
            # All socket connections will raise OSError
            pass
    """

    def block_socket(*args, **kwargs):
        raise OSError("Network is disabled for this test")

    # Block socket creation
    monkeypatch.setattr(socket, "socket", block_socket)

    # Also block socket.create_connection
    monkeypatch.setattr(socket, "create_connection", block_socket)

    yield


@pytest.fixture
def enable_germlines(monkeypatch):
    """Fixture to enable germlines module via feature flag.

    Sets SADIE_USE_GERMLINES_MODULE=true environment variable.
    """
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")
    yield


@pytest.fixture
def disable_germlines(monkeypatch):
    """Fixture to disable germlines module (use G3 API).

    Sets SADIE_USE_GERMLINES_MODULE=false environment variable.
    """
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "false")
    yield


@pytest.fixture
def mock_unpopulated_germlines(monkeypatch, tmp_path):
    """Fixture to simulate unpopulated germlines database.

    Creates a mock germlines directory without data files to test
    error handling for first-time setup scenarios.
    """
    # Create empty germlines directory structure
    empty_sources = tmp_path / "germlines" / "sources"
    empty_sources.mkdir(parents=True)

    # Patch the germlines base directory to point to empty location
    def mock_base_dir():
        return tmp_path / "germlines"

    monkeypatch.setattr(
        "sadie.germlines.get_germlines_base_dir",
        mock_base_dir
    )

    yield tmp_path / "germlines"
