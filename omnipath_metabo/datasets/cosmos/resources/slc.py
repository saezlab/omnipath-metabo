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
SLC transporter processing for COSMOS PKN.

The SLC (Solute Carrier) table provides human transporter-metabolite
interactions with subcellular localization information.
"""

from __future__ import annotations

__all__ = ['slc_interactions']

import re
from collections.abc import Generator

from .._record import Interaction


def slc_interactions() -> Generator[Interaction, None, None]:
    """
    Yield SLC transporter-substrate interactions as uniform records.

    Compound location strings (e.g. ``"ER; Plasma membrane"``) are split
    into individual locations and each is mapped to a compartment code.
    The resulting record contains all unique codes as a tuple.

    Note: SLC data is only available for human (NCBI taxonomy ID 9606).

    Yields:
        :class:`Interaction` records with *source_type*
        ``'small_molecule'`` and *target_type* ``'protein'``.
    """

    from pypath.inputs import slc

    from ..location import slc_locations

    location_mapping = slc_locations()

    for r in slc.slc_interactions():

        chebi = r.substrate.chebi
        localization = r.localization

        if not chebi or not localization or localization in ('', 'Unknown', 'None'):
            continue

        parts = re.split(r';\s*', localization)
        codes = {
            location_mapping[part]
            for part in parts
            if part in location_mapping
        }

        if not codes:
            continue

        yield Interaction(
            source=chebi,
            target=r.transporter.uniprot,
            source_type='small_molecule',
            target_type='protein',
            id_type_a='chebi',
            id_type_b='uniprot',
            interaction_type='transport',
            resource='SLC',
            mor=0,
            locations=tuple(sorted(codes)),
        )
