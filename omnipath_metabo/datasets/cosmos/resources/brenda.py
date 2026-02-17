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
from typing import TYPE_CHECKING

from .._record import Interaction

if TYPE_CHECKING:
    from collections.abc import Sequence


def brenda_regulations(
    organisms: Sequence[str] | None = None,
) -> Generator[Interaction, None, None]:
    """
    Yield BRENDA allosteric regulation interactions as uniform records.

    Args:
        organisms:
            List of organism names (e.g. ``['human']``).
            Defaults to ``['human']``.

    Yields:
        :class:`Interaction` records with *source_type*
        ``'small_molecule'`` and *target_type* ``'protein'``.
    """

    from pypath.inputs.brenda._main import allosteric_regulation

    if organisms is None:
        organisms = ['human']

    for record in allosteric_regulation(organisms=organisms, limit=None):

        if not record.protein:
            continue

        mor = (
            1 if record.action == 'activating' else
            -1 if record.action == 'inhibiting' else
            0
        )

        for protein_id in record.protein:
            yield Interaction(
                source=record.compound,
                target=protein_id,
                source_type='small_molecule',
                target_type='protein',
                id_type_a=record.id_type,
                id_type_b='uniprot',
                interaction_type='regulation',
                resource='BRENDA',
                mor=mor,
            )
