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
STITCH chemical-protein interactions for COSMOS PKN.

STITCH integrates chemical-protein interaction data from multiple
evidence channels. This module filters by confidence score and
interaction mode, normalising orientation so the small molecule
is always the source.
"""

from __future__ import annotations

__all__ = ['stitch_interactions']

from collections.abc import Generator
from typing import TYPE_CHECKING

from .._record import Interaction

if TYPE_CHECKING:
    from collections.abc import Sequence


def stitch_interactions(
    organism: int = 9606,
    score_threshold: int = 700,
    mode: str | Sequence[str] | None = ('activation', 'inhibition'),
) -> Generator[Interaction, None, None]:
    """
    Yield STITCH chemical-protein interactions as uniform records.

    Args:
        organism:
            NCBI taxonomy ID (default: 9606 for human).
        score_threshold:
            Minimum ``final_score`` for interactions.
        mode:
            Interaction mode(s) to keep.  Pass a single string, a
            sequence of strings, or ``None`` to keep all modes.
            Available modes per record: ``'binding'``, ``'pred_bind'``,
            ``'expression'``, ``'activation'``, ``'inhibition'``,
            ``'reaction'``, ``'catalysis'``.

    Yields:
        :class:`Interaction` records with *source_type*
        ``'small_molecule'`` and *target_type* ``'protein'``.
    """

    from pypath.inputs.new_stitch import interactions

    allowed = None

    if mode is not None:
        allowed = (mode,) if isinstance(mode, str) else tuple(mode)

    for rec in interactions(ncbi_tax_id=organism):

        if rec.final_score < score_threshold:
            continue

        if allowed and rec.mode not in allowed:
            continue

        # Normalise orientation: small molecule as source
        if rec.source.type == 'small_molecule':
            chemical_id = rec.source.id
            protein_id = rec.target.id
        elif rec.target.type == 'small_molecule':
            chemical_id = rec.target.id
            protein_id = rec.source.id
        else:
            continue

        if rec.activation:
            mor = 1
        elif rec.inhibition:
            mor = -1
        else:
            mor = 0

        yield Interaction(
            source=chemical_id,
            target=protein_id,
            source_type='small_molecule',
            target_type='protein',
            id_type_a='pubchem',
            id_type_b='ensp',
            interaction_type='unknown',
            resource='STITCH',
            mor=mor,
        )
