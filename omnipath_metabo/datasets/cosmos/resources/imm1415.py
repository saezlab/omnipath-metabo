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
iMM1415 transporter and metabolic interactions for COSMOS PKN (mouse).

Thin wrappers around
:func:`pypath.inputs.imm1415.imm1415_transporter_network` and
:func:`pypath.inputs.imm1415.imm1415_network` that convert
``GemInteraction`` records into COSMOS ``Interaction`` records.

Only Mus musculus (10090) is supported; human builds yield nothing.
"""

from __future__ import annotations

__all__ = ['imm1415_transporter_interactions', 'imm1415_metabolic_interactions']

from collections import defaultdict
from collections.abc import Generator

from .._record import Interaction


def imm1415_transporter_interactions(
    organism: int = 9606,
    include_reverse: bool = True,
    include_orphans: bool = False,
    cell_surface_only: bool = False,
    **_kwargs,
) -> Generator[Interaction, None, None]:
    """
    Yield iMM1415 transporter interactions as uniform Interaction records.

    Mouse-only counterpart of
    :func:`~omnipath_metabo.datasets.cosmos.resources.recon3d.recon3d_transporter_interactions`.
    Delegates transport detection to
    :func:`pypath.inputs.imm1415.imm1415_transporter_network`.

    Only Mus musculus (10090) is supported; other organisms yield nothing.

    Args:
        organism:
            NCBI taxonomy ID.  Only 10090 (mouse) is supported.
        include_reverse:
            If ``True``, include reversed edges for reversible transport
            reactions.  Default: ``True``.
        include_orphans:
            If ``True``, include transport reactions with no gene rule,
            using the reaction ID as a pseudo-enzyme node.  Default: ``False``.
        cell_surface_only:
            If ``True``, restrict to transport events where at least one
            compartment is extracellular (``'e'``).  Applied per reaction
            group.  Default: ``False``.

    Yields:
        :class:`~omnipath_metabo.datasets.cosmos._record.Interaction`
        records with ``interaction_type='transport'`` and
        ``resource='iMM1415'``.  Metabolite IDs use ``id_type='bigg'``;
        protein IDs use ``id_type='entrez'`` (mouse NCBI Entrez).

    Warning:
        Reverse edges are already included when ``include_reverse=True``.
        A downstream formatter must NOT re-generate reverse edges.  Use
        ``attrs['reverse']`` to identify and label them.
    """

    if organism != 10090:
        return

    from pypath.inputs.imm1415._gem import imm1415_transporter_network

    groups: defaultdict[str, list] = defaultdict(list)

    for rec in imm1415_transporter_network(
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

            is_complex = (
                '_' in (rec.source if source_type == 'protein' else rec.target)
                and not is_orphan
            )

            attrs = {
                'reverse': rec.reverse,
                'reaction_id': rec.reaction_id,
                'enzyme_complex': is_complex,
                'transport_from': (
                    rec.source_compartment if rec.source_type == 'metabolite'
                    else rec.target_compartment
                ),
                'transport_to': (
                    rec.target_compartment if rec.target_type == 'metabolite'
                    else rec.source_compartment
                ),
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
                resource='iMM1415',
                mor=1,
                locations=locations,
                attrs=attrs,
            )


def imm1415_metabolic_interactions(
    organism: int = 9606,
    include_reverse: bool = True,
    **_kwargs,
) -> Generator[Interaction, None, None]:
    """
    Yield iMM1415 stoichiometric enzyme-metabolite interactions.

    Mouse-only counterpart of
    :func:`~omnipath_metabo.datasets.cosmos.resources.recon3d.recon3d_metabolic_interactions`.
    Delegates to :func:`pypath.inputs.imm1415.imm1415_network`.

    Only Mus musculus (10090) is supported; other organisms yield nothing.

    Args:
        organism:
            NCBI taxonomy ID.  Only 10090 (mouse) is supported.
        include_reverse:
            If ``True``, include reversed edges for reversible reactions.
            Default: ``True``.

    Yields:
        :class:`~omnipath_metabo.datasets.cosmos._record.Interaction`
        records with ``interaction_type='catalysis'`` and
        ``resource='GEM:iMM1415'``.  Metabolite IDs use ``id_type='bigg'``;
        protein IDs use ``id_type='entrez'`` (mouse NCBI Entrez).

    Warning:
        Reverse edges are already included when ``include_reverse=True``.
        A downstream formatter must NOT re-generate reverse edges.
    """

    if organism != 10090:
        return

    from pypath.inputs.imm1415._gem import imm1415_network

    for rec in imm1415_network(include_reverse=include_reverse):

        if rec.source_type == 'metabolite':
            source = rec.source
            source_type = 'small_molecule'
            id_type_a = 'bigg'
            target = rec.target
            target_type = 'protein'
            id_type_b = 'entrez'
            locations = (rec.source_compartment,) if rec.source_compartment else ()
        else:
            source = rec.source
            source_type = 'protein'
            id_type_a = 'entrez'
            target = rec.target
            target_type = 'small_molecule'
            id_type_b = 'bigg'
            locations = (rec.target_compartment,) if rec.target_compartment else ()

        yield Interaction(
            source=source,
            target=target,
            source_type=source_type,
            target_type=target_type,
            id_type_a=id_type_a,
            id_type_b=id_type_b,
            interaction_type='catalysis',
            resource='GEM:iMM1415',
            mor=1,
            locations=locations,
            attrs={
                'reverse': rec.reverse,
                'reaction_id': rec.reaction_id,
            },
        )
