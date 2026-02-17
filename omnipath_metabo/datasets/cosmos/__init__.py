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
data sources including TCDB, SLC, BRENDA, MRCLinksDB, STITCH, and Rhea.

Each resource processor is a generator that yields
:class:`_record.Interaction` named tuples in a uniform format.

Usage — test individual resources::

    from omnipath_metabo.datasets.cosmos.resources import (
        stitch_interactions,
        tcdb_interactions,
        slc_interactions,
        brenda_regulations,
        mrclinksdb_interactions,
    )

    # STITCH (fast, no location lookup)
    stitch = list(stitch_interactions())
    len(stitch)
    stitch[:3]

    # SLC (fast, human only)
    slc = list(slc_interactions())
    len(slc)
    slc[:3]

    # BRENDA
    brenda = list(brenda_regulations())
    len(brenda)
    brenda[:3]

    # TCDB (slower — needs UniProt location query)
    tcdb = list(tcdb_interactions())
    len(tcdb)
    tcdb[:3]

    # MRCLinksDB (slower — needs UniProt location query)
    mrc = list(mrclinksdb_interactions())
    len(mrc)
    mrc[:3]

    # Collect all into a DataFrame
    import pandas as pd
    from omnipath_metabo.datasets.cosmos._record import Interaction
    all_records = stitch + slc + brenda + tcdb + mrc
    df = pd.DataFrame(all_records, columns=Interaction._fields)
    df
"""

__all__ = [
    'build',
    'resources',
]

from . import resources
from ._build import build
