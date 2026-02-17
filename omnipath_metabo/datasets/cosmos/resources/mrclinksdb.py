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


def mrclinksdb_interactions(
    organism: str = 'human',
) -> Generator[Interaction, None, None]:
    """
    Yield MRCLinksDB receptor-metabolite interactions as uniform records.

    Args:
        organism:
            Organism name (e.g. ``'human'``, ``'mouse'``).

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

    location_mapping = tcdb_locations()

    ncbi_tax_id = 9606 if organism == 'human' else None
    all_locations = (
        uniprot_locations(organism=ncbi_tax_id, reviewed=True)
        if ncbi_tax_id
        else {}
    )

    for rec in _interactions.mrclinksdb_interaction(organism=organism):

        receptor = str(rec.receptor_uniprot)
        pubchem = rec.pubchem

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
            mor=0,
            locations=tuple(sorted(abbreviations)),
        )
