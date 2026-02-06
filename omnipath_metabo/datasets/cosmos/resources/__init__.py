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

"""COSMOS PKN data resource processors."""

__all__ = [
    'tcdb_interactions',
    'slc_interactions',
    'brenda_regulations',
    'mrclinksdb_interactions',
    'rhea_reactions',
]

from .tcdb import tcdb_interactions
from .slc import slc_interactions
from .brenda import brenda_regulations
from .mrclinksdb import mrclinksdb_interactions
from .rhea import rhea_reactions
