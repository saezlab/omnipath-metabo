#!/usr/bin/env python

#
# This file is part of the `omnipath_metabo` Python module
#
# Copyright 2026
# Heidelberg University Hospital
#
# File author(s): OmniPath Team (omnipathdb@gmail.com)
#
# Distributed under the BSD-3-Clause license
# See the file `LICENSE` or read a copy at
# https://opensource.org/license/bsd-3-clause
#

"""COSMOS-specific data files."""

from pathlib import Path


DATA_DIR = Path(__file__).parent


def data_path(filename: str) -> Path:
    """
    Path to a data file in the COSMOS data directory.

    Args:
        filename: Name of the data file.

    Returns:
        Path to the data file.
    """

    return DATA_DIR / filename
