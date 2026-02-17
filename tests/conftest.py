#!/usr/bin/env python

"""Shared test fixtures for omnipath_metabo tests."""

import pytest


@pytest.fixture
def tmp_data_dir(tmp_path):
    """Create a temporary data directory for test files."""

    data_dir = tmp_path / 'data'
    data_dir.mkdir()
    return data_dir
