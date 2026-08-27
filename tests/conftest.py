"""Shared pytest configuration for NGP Sovereign Synesis bounty test suites."""
import pathlib

import pytest


def pytest_configure(config):
    pathlib.Path("results").mkdir(exist_ok=True)


def pytest_collection_modifyitems(config, items):
    for item in items:
        if "docker" in item.name.lower():
            item.add_marker(pytest.mark.integration)
        if "determinism" in item.name.lower():
            item.add_marker(pytest.mark.slow)
