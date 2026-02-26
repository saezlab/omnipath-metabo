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

    # Or use build() to collect all into a single DataFrame
    from omnipath_metabo.datasets.cosmos import build
    df = build()
    df
"""

__all__ = [
    'build',
    'build_allosteric',
    'build_enzyme_metabolite',
    'build_receptors',
    'build_transporters',
    'config',
    'default_config',
    'format_pkn',
    'resources',
]

from . import resources
from ._build import build, build_allosteric, build_enzyme_metabolite, build_receptors, build_transporters
from ._config import config, default_config
from ._format import format_pkn
