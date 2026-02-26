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
MRCLinksDB receptor-metabolite processing for COSMOS PKN.

MRCLinksDB provides receptor-metabolite interaction data with
subcellular localization information.
"""

from __future__ import annotations

__all__ = ['mrclinksdb_interactions']

from collections.abc import Generator

from .._record import Interaction

ORGANISM_NAMES = {
    9606: 'human',
    10090: 'mouse',
}


def mrclinksdb_interactions(
    organism: int = 9606,
) -> Generator[Interaction, None, None]:
    """
    Yield MRCLinksDB receptor-metabolite interactions as uniform records.

    Args:
        organism:
            NCBI taxonomy ID (default: 9606 for human).

    Yields:
        :class:`Interaction` records with *source_type*
        ``'small_molecule'`` and *target_type* ``'protein'``.
    """

    from pypath.inputs.mrclinksdb import _interactions

    from ..location import (
        locations_to_abbreviations,
        tcdb_locations,
        uniprot_locations,
    )

    organism_name = ORGANISM_NAMES.get(organism, str(organism))
    location_mapping = tcdb_locations()
    all_locations = uniprot_locations(organism=organism, reviewed=True)

    for rec in _interactions.mrclinksdb_interaction(organism=organism_name):

        receptor = str(rec.receptor_uniprot)
        pubchem_raw = str(rec.pubchem)

        # PubChem IDs arrive with a 'CID:' prefix (e.g. 'CID:13712').
        # Records that do not match this pattern (including the CSV header
        # row that pypath passes through) are skipped.
        if not pubchem_raw.startswith('CID:'):
            continue

        pubchem = pubchem_raw[4:]  # strip 'CID:' → plain numeric string

        if not receptor or not pubchem:
            continue

        if receptor not in all_locations:
            continue

        abbreviations = locations_to_abbreviations(
            all_locations[receptor],
            location_mapping,
        )

        if not abbreviations:
            continue

        yield Interaction(
            source=pubchem,
            target=receptor,
            source_type='small_molecule',
            target_type='protein',
            id_type_a='pubchem',
            id_type_b='uniprot',
            interaction_type='ligand_receptor',
            resource='MRCLinksDB',
            mor=1,
            locations=tuple(sorted(abbreviations)),
        )
