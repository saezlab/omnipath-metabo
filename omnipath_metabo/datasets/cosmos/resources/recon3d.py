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
Recon3D transporter interactions for COSMOS PKN.

Thin wrapper around
:func:`pypath.inputs.recon3d.recon3d_transporter_network` that converts
``GemInteraction`` records into COSMOS ``Interaction`` records.

The transport detection algorithm (compartment-crossing filter) lives in
pypath.  This module handles COSMOS-specific concerns: the
``cell_surface_only`` filter, ``Interaction`` record building, and
``id_type`` / resource labelling.

For reference on the transport detection approach and the original R
implementations, see the docstring of
``pypath.inputs.recon3d._gem.recon3d_transporter_network``.
"""

from __future__ import annotations

__all__ = ['recon3d_transporter_interactions']

from collections import defaultdict
from collections.abc import Generator

from .._record import Interaction


def recon3d_transporter_interactions(
    organism: int = 9606,
    include_reverse: bool = True,
    include_orphans: bool = True,
    cell_surface_only: bool = False,
) -> Generator[Interaction, None, None]:
    """
    Yield Recon3D transporter interactions as uniform Interaction records.

    Delegates transport detection to
    :func:`pypath.inputs.recon3d.recon3d_transporter_network`.  The
    compartment-crossing filter (keeping only metabolites that physically
    cross a membrane) is applied in pypath.

    Two directed edges are generated per transported metabolite per enzyme:

    - ``met[in_comp] → enzyme``: the metabolite enters the transporter.
    - ``enzyme → met[out_comp]``: the metabolite exits on the other side.

    Args:
        organism:
            NCBI taxonomy ID.  Only human (9606) is supported; any other
            value causes the function to yield nothing.
        include_reverse:
            If ``True``, include reversed edges for reversible transport
            reactions (``attrs['reverse'] = True``).  Default: ``True``.
        include_orphans:
            If ``True``, include transport reactions with no gene rule,
            using the reaction ID as a pseudo-enzyme node with
            ``id_type = 'reaction_id'`` and ``attrs['orphan'] = True``.
            Default: ``True``.
        cell_surface_only:
            If ``True``, restrict to transport events where at least one
            compartment is extracellular (``'e'``).  Applied per reaction
            group so met→enzyme / enzyme→met pairs are always kept or
            dropped together.  Default: ``False``.

    Yields:
        :class:`~omnipath_metabo.datasets.cosmos._record.Interaction`
        records.  For metabolite → enzyme edges, ``source_type`` is
        ``'small_molecule'`` and ``id_type_a`` is ``'bigg'``; ``target_type``
        is ``'protein'`` and ``id_type_b`` is ``'entrez'``.  Roles are
        swapped for enzyme → metabolite edges.  ``interaction_type`` is
        ``'transport'`` and ``resource`` is ``'Recon3D'``.

    Warning:
        Reverse edges are already included in the output when
        ``include_reverse=True`` (the default).  A downstream COSMOS
        formatter must NOT re-generate reverse edges.  Use
        ``attrs['reverse']`` to identify and label them, not to create
        new edges from them.
    """

    if organism != 9606:
        return

    from pypath.inputs.recon3d._gem import recon3d_transporter_network

    # Group by reaction_id for cell_surface_only filter.
    groups: defaultdict[str, list] = defaultdict(list)

    for rec in recon3d_transporter_network(
        include_reverse=include_reverse,
        include_orphans=include_orphans,
    ):
        groups[rec.reaction_id].append(rec)

    for rxn_id, edges in groups.items():

        if cell_surface_only:
            if not any(
                e.source_compartment == 'e' or e.target_compartment == 'e'
                for e in edges
            ):
                continue

        for rec in edges:

            is_orphan = rec.source_type == 'reaction' or rec.target_type == 'reaction'

            if rec.source_type == 'metabolite':
                source = rec.source
                source_type = 'small_molecule'
                id_type_a = 'bigg'
                target = rec.target
                target_type = 'protein'
                id_type_b = 'reaction_id' if is_orphan else 'entrez'
                locations = (rec.source_compartment,) if rec.source_compartment else ()

            else:
                source = rec.source
                source_type = 'protein'
                id_type_a = 'reaction_id' if is_orphan else 'entrez'
                target = rec.target
                target_type = 'small_molecule'
                id_type_b = 'bigg'
                locations = (rec.target_compartment,) if rec.target_compartment else ()

            is_complex = '_' in (rec.source if source_type == 'protein' else rec.target) and not is_orphan

            attrs = {
                'reverse': rec.reverse,
                'reaction_id': rec.reaction_id,
                'enzyme_complex': is_complex,
                'transport_from': rec.source_compartment if rec.source_type == 'metabolite' else rec.target_compartment,
                'transport_to': rec.target_compartment if rec.target_type == 'metabolite' else rec.source_compartment,
            }

            if is_orphan:
                attrs['orphan'] = True

            yield Interaction(
                source=source,
                target=target,
                source_type=source_type,
                target_type=target_type,
                id_type_a=id_type_a,
                id_type_b=id_type_b,
                interaction_type='transport',
                resource='Recon3D',
                mor=1,
                locations=locations,
                attrs=attrs,
            )
