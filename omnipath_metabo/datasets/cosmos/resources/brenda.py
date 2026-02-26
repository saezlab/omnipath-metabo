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
BRENDA allosteric regulation processing for COSMOS PKN.

BRENDA provides enzyme allosteric regulation data, including
activators and inhibitors.
"""

from __future__ import annotations

__all__ = ['brenda_regulations']

from collections.abc import Generator

from .._record import Interaction

ORGANISM_NAMES = {
    9606: 'human',
    10090: 'mouse',
}


def brenda_regulations(
    organism: int = 9606,
) -> Generator[Interaction, None, None]:
    """
    Yield BRENDA allosteric regulation interactions as uniform records.

    Args:
        organism:
            NCBI taxonomy ID (default: 9606 for human).

    Yields:
        :class:`Interaction` records with *source_type*
        ``'small_molecule'`` and *target_type* ``'protein'``.
    """

    from pypath.inputs.brenda._main import allosteric_regulation

    organism_name = ORGANISM_NAMES.get(organism, str(organism))

    for record in allosteric_regulation(
        organisms=[organism_name],
        limit=None,
    ):

        if not record.protein:
            continue

        mor = (
            1 if record.action == 'activator' else
            -1 if record.action == 'inhibitor' else
            0
        )

        if mor == 0:
            continue

        for protein_id in record.protein:
            yield Interaction(
                source=record.compound,
                target=protein_id,
                source_type='small_molecule',
                target_type='protein',
                id_type_a='synonym',
                id_type_b=record.id_type,
                interaction_type='allosteric_regulation',
                resource='BRENDA',
                mor=mor,
            )
