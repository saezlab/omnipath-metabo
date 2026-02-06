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

"""
COSMOS prior-knowledge network builder.

This module contains all logic for building the COSMOS PKN from multiple
data sources including TCDB, SLC, BRENDA, MRCLinksDB, and Rhea.
"""

__all__ = [
    'build',
    'sources',
]

from . import sources
from ._build import build
