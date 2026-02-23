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
    'brenda_regulations',
    'gem_interactions',
    'mrclinksdb_interactions',
    'rhea_reactions',
    'slc_interactions',
    'stitch_interactions',
    'tcdb_interactions',
]

from .brenda import brenda_regulations
from .gem import gem_interactions
from .mrclinksdb import mrclinksdb_interactions
from .rhea import rhea_reactions
from .slc import slc_interactions
from .stitch import stitch_interactions
from .tcdb import tcdb_interactions
