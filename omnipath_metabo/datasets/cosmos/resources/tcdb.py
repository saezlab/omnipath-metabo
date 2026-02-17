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
TCDB transporter processing for COSMOS PKN.

TCDB (Transporter Classification Database) is a multi-species database
with transporter substrate information. This module extracts species-specific
transporter-metabolite interactions.
"""

from __future__ import annotations

__all__ = ['tcdb_interactions']

from collections.abc import Generator

from .._record import Interaction


def tcdb_interactions(
    ncbi_tax_id: int = 9606,
) -> Generator[Interaction, None, None]:
    """
    Yield TCDB transporter-substrate interactions as uniform records.

    TCDB is a multi-species database where entries from different species
    are mixed together. This function extracts species-specific interactions
    based on protein UniProt IDs and queries UniProt for their subcellular
    locations.

    Args:
        ncbi_tax_id:
            NCBI taxonomy ID (default: 9606 for human).

    Yields:
        :class:`Interaction` records with *source_type*
        ``'small_molecule'`` and *target_type* ``'protein'``.
    """

    from pypath.inputs.tcdb import _substrates as tcdb_substrates
    from pypath.utils import reflists

    from ..location import (
        locations_to_abbreviations,
        tcdb_locations,
        uniprot_locations,
    )

    location_mapping = tcdb_locations()
    species_proteins = set(
        reflists.get_reflist('uniprot', ncbi_tax_id=ncbi_tax_id)
    )
    all_locations = uniprot_locations(organism=ncbi_tax_id, reviewed=True)

    for r in tcdb_substrates.tcdb_substrate():

        if r.transporter_uniprot not in species_proteins:
            continue

        if r.transporter_uniprot not in all_locations:
            continue

        abbreviations = locations_to_abbreviations(
            all_locations[r.transporter_uniprot],
            location_mapping,
        )

        if not abbreviations:
            continue

        yield Interaction(
            source=r.substrate_id,
            target=r.transporter_uniprot,
            source_type='small_molecule',
            target_type='protein',
            id_type_a='tcdb',
            id_type_b='uniprot',
            interaction_type='transport',
            resource='TCDB',
            mor=0,
            locations=tuple(sorted(abbreviations)),
        )
