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
Genome-scale metabolic model (GEM) processing for COSMOS PKN.

Thin wrapper around ``pypath.inputs.metatlas`` parsers that converts
``GemInteraction`` records into COSMOS ``Interaction`` records.

Metabolic edges come from
:func:`pypath.inputs.metatlas.metatlas_gem_network`; transport edges come
from :func:`pypath.inputs.metatlas.metatlas_gem_transport_network`, which
already applies the compartment-crossing filter.

Transport reactions (subsystem ``'Transport reactions'``) are labelled
``resource='GEM_transporter:<gem_name>'``; all other reactions use
``resource='GEM:<gem_name>'``.  This allows the two classes to be
filtered, weighted, or disabled independently in downstream analyses.
"""

from __future__ import annotations

__all__ = ['gem_interactions']

from collections import defaultdict
from collections.abc import Generator

from .._record import Interaction


def gem_interactions(
    gem: str | list[str] = 'Human-GEM',
    organism: int = 9606,
    include_reverse: bool = True,
    include_orphans: bool = True,
) -> Generator[Interaction, None, None]:
    """
    Yield GEM metabolite-enzyme interactions as uniform records.

    Metabolic and transport reactions are sourced separately via pypath:
    :func:`~pypath.inputs.metatlas.metatlas_gem_network` provides metabolic
    edges (all reactions except transport); the transport edges come from
    :func:`~pypath.inputs.metatlas.metatlas_gem_transport_network`, which
    already applies the compartment-crossing filter.

    Transport reactions are labelled
    ``resource='GEM_transporter:<gem_name>'``; metabolic reactions use
    ``resource='GEM:<gem_name>'``.

    Args:
        gem:
            GEM name or list of GEM names to load and merge.  Names must
            match entries in the Metabolic Atlas validation index (e.g.
            ``'Human-GEM'``, ``'Mouse-GEM'``).
        organism:
            NCBI taxonomy ID.  Accepted for API consistency; GEMs are
            inherently organism-specific so this parameter is not used
            for filtering.
        include_reverse:
            If ``True``, include reversed edges for reversible reactions
            (``attrs['reverse'] = True``).  Default: ``True``.
        include_orphans:
            If ``True``, include reactions with no gene rule, using the
            reaction ID as a pseudo-enzyme node with
            ``id_type = 'reaction_id'`` and ``attrs['orphan'] = True``.
            Default: ``True``.

    Yields:
        :class:`Interaction` records.  For metabolite → enzyme edges,
        ``source_type`` is ``'small_molecule'`` and ``target_type`` is
        ``'protein'``.  For enzyme → metabolite edges the roles are
        swapped.  ``id_type_a`` / ``id_type_b`` are the GEM's native
        gene ID type (``'ensembl'`` or ``'genesymbol'``) for enzymes,
        ``'metatlas'`` for metabolites, and ``'reaction_id'`` for orphan
        pseudo-enzyme nodes.

        Each yielded record carries ``attrs['gems']``: a sorted list of
        GEM names that contributed the edge.

    Warning:
        Reverse edges are already included in the output when
        ``include_reverse=True`` (the default).  A downstream COSMOS
        formatter must NOT re-generate reverse edges.  Use
        ``attrs['reverse']`` to identify and label them, not to create
        new edges from them.
    """

    from pypath.inputs.metatlas._gem import (
        metatlas_gem_network,
        metatlas_gem_transport_network,
        metatlas_gem_transport_ids,
        metatlas_gem_detect_gene_id_type,
        _strip_compartment,
    )

    gems = [gem] if isinstance(gem, str) else list(gem)

    # Pre-detect gene ID type per GEM.
    gene_id_types: dict[str, str] = {
        gem_name: metatlas_gem_detect_gene_id_type(gem_name)
        for gem_name in gems
    }

    # Collect edges in two labelled streams:
    #   - metabolic_raw: list of (GemInteraction, gem_name, resource_prefix)
    #   - transport_groups: (rxn_id, gem_name) → list[GemInteraction]
    metabolic_raw: list[tuple] = []
    transport_groups: dict[tuple, list] = defaultdict(list)

    for gem_name in gems:

        t_ids = metatlas_gem_transport_ids(gem_name)

        # Metabolic edges: all reactions except transport.
        for rec in metatlas_gem_network(
            gem=gem_name,
            include_orphans=include_orphans,
        ):
            if not include_reverse and rec.reverse:
                continue

            if rec.reaction_id in t_ids:
                continue

            metabolic_raw.append((rec, gem_name))

        # Transport edges: grouped by reaction for cell_surface_only filter.
        for rec in metatlas_gem_transport_network(
            gem=gem_name,
            include_reverse=include_reverse,
            include_orphans=include_orphans,
        ):
            transport_groups[(rec.reaction_id, gem_name)].append(rec)

    transport_raw: list[tuple] = [
        (rec, gem_name)
        for (_, gem_name), edges in transport_groups.items()
        for rec in edges
    ]

    # Build Interaction objects from both streams.
    # Tag each record with its resource prefix before merging.
    all_records: list[tuple] = (
        [(rec, gem_name, 'GEM') for rec, gem_name in metabolic_raw] +
        [(rec, gem_name, 'GEM_transporter') for rec, gem_name in transport_raw]
    )

    all_interactions: list[tuple] = []

    for rec, gem_name, resource_prefix in all_records:

        is_orphan = rec.source_type == 'reaction' or rec.target_type == 'reaction'
        gene_id_type = gene_id_types[gem_name]

        if rec.source_type == 'metabolite':
            met_id = _strip_compartment(rec.source, rec.source_compartment)
            enzyme_id = rec.target
            compartment = rec.source_compartment
            source = met_id
            target = enzyme_id
            source_type = 'small_molecule'
            target_type = 'protein'
            id_type_a = 'metatlas'
            id_type_b = 'reaction_id' if is_orphan else gene_id_type

        else:
            enzyme_id = rec.source
            met_id = _strip_compartment(rec.target, rec.target_compartment)
            compartment = rec.target_compartment
            source = enzyme_id
            target = met_id
            source_type = 'protein'
            target_type = 'small_molecule'
            id_type_a = 'reaction_id' if is_orphan else gene_id_type
            id_type_b = 'metatlas'

        locations = (compartment,) if compartment else ()
        is_complex = '_' in enzyme_id and not is_orphan

        attrs = {
            'reverse': rec.reverse,
            'reaction_id': rec.reaction_id,
            'enzyme_complex': is_complex,
        }

        if is_orphan:
            attrs['orphan'] = True

        all_interactions.append((
            Interaction(
                source=source,
                target=target,
                source_type=source_type,
                target_type=target_type,
                id_type_a=id_type_a,
                id_type_b=id_type_b,
                interaction_type='catalysis',
                resource=f'{resource_prefix}:{gem_name}',
                mor=1,
                locations=locations,
                attrs=attrs,
            ),
            gem_name,
        ))

    # Deduplicate across GEMs: merge edges with same (source, target,
    # reaction_id, reverse) key, collecting all GEM names in attrs['gems'].
    seen: dict[tuple, Interaction] = {}
    gem_provenance: dict[tuple, list[str]] = defaultdict(list)

    for interaction, gem_name in all_interactions:
        key = (
            interaction.source,
            interaction.target,
            interaction.attrs.get('reaction_id', ''),
            interaction.attrs.get('reverse', False),
        )

        if key not in seen:
            seen[key] = interaction

        gem_provenance[key].append(gem_name)

    for key, interaction in seen.items():
        gem_names = sorted(set(gem_provenance[key]))
        attrs = dict(interaction.attrs)
        attrs['gems'] = gem_names
        yield interaction._replace(attrs=attrs)
