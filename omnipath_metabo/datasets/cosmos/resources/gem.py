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

Fetches reaction networks from standard GEM repositories via
``pypath.inputs.metatlas`` and converts them into binary directed edges:
metabolite → enzyme and enzyme → metabolite.  Reversible reactions produce
additional edges with ``attrs['reverse'] = True``.

High-degree metabolites (cofactors such as ATP, NADH, H2O) are removed
using a degree threshold applied after all edges are collected.

Metabolite IDs are in the MetAtlas format (e.g. ``MAM01039``); enzyme IDs
are Ensembl gene IDs (e.g. ``ENSG00000115275``) or compound complex IDs
(e.g. ``ENSG00000115275_ENSG00000136881`` for AND-rule subunits).

Transport reactions (Human-GEM subsystem ``'Transport reactions'``) are
labelled ``resource='GEM_transporter:<gem_name>'``; all other reactions use
``resource='GEM:<gem_name>'``.  This allows the two classes to be filtered,
weighted, or disabled independently in downstream analyses.
"""

from __future__ import annotations

__all__ = ['gem_interactions']

from collections import Counter
from collections.abc import Generator
from typing import TYPE_CHECKING

from .._record import Interaction

if TYPE_CHECKING:
    pass


def _transport_reaction_ids(gem_name: str) -> frozenset[str]:
    """
    Return the set of reaction IDs classified as transport reactions in a GEM.

    Uses ``metatlas_gem_yaml_reactions`` to read the YAML manifest (already
    cached after ``metatlas_gem_network`` fetches it) and collects every
    reaction whose ``subsystem`` field equals ``'Transport reactions'``.

    Args:
        gem_name: GEM name (e.g. ``'Human-GEM'``).

    Returns:
        Frozenset of reaction ID strings (e.g. ``{'MAR01234', ...}``).
    """

    from pypath.inputs.metatlas._gem import metatlas_gem_yaml_reactions

    return frozenset(
        rxn.id
        for rxn in metatlas_gem_yaml_reactions(gem=gem_name)
        if rxn.subsystem == 'Transport reactions'
    )


def _strip_compartment(met_id: str, compartment: str) -> str:
    """
    Remove compartment suffix from a metabolite ID.

    Human-GEM encodes the compartment as a single-letter suffix on the
    metabolite ID (e.g. ``MAM01039c`` → ``MAM01039``).  This function
    strips that suffix when it matches ``compartment``.

    Args:
        met_id: Full metabolite ID including compartment suffix.
        compartment: Compartment code to strip (e.g. ``'c'``).

    Returns:
        Base metabolite ID without compartment suffix, or the original
        ID if the suffix does not match.
    """

    if compartment and met_id.endswith(compartment):
        return met_id[: -len(compartment)]

    return met_id


def gem_interactions(
    gem: str | list[str] = 'Human-GEM',
    organism: int = 9606,
    metab_max_degree: int = 400,
    include_reverse: bool = True,
    include_orphans: bool = True,
) -> Generator[Interaction, None, None]:
    """
    Yield GEM metabolite-enzyme interactions as uniform records.

    Reactions are decomposed into binary edges:

    - *metabolite → enzyme*: reactant consumed by the enzyme.
    - *enzyme → metabolite*: product released by the enzyme.

    Reversible reactions additionally produce a reversed pair with
    ``attrs['reverse'] = True``.

    High-degree metabolites (likely cofactors) are filtered out: any
    metabolite appearing in more than *metab_max_degree* edges across the
    full collected edge set is excluded.

    Orphan reactions (no gene rule) are retained when *include_orphans*
    is ``True``.  The reaction ID is used as a pseudo-enzyme node with
    ``id_type = 'reaction_id'`` and ``attrs['orphan'] = True``.  These
    edges pass through ID translation unchanged and can be identified and
    filtered downstream by ``attrs['orphan']``.

    Args:
        gem:
            GEM name or list of GEM names to load and merge.  Names must
            match entries in the Metabolic Atlas validation index (e.g.
            ``'Human-GEM'``, ``'Mouse-GEM'``).
        organism:
            NCBI taxonomy ID.  Accepted for API consistency with other
            resources; GEMs are inherently organism-specific so this
            parameter is not used for filtering.
        metab_max_degree:
            Maximum number of edges a metabolite may participate in.
            Metabolites exceeding this threshold are treated as cofactors
            and removed.  Default: 400 (matches OmnipathR reference).
        include_reverse:
            If ``True``, include reversed edges for reversible reactions
            (``attrs['reverse'] = True``).  Default: ``True``.
        include_orphans:
            If ``True``, include reactions with no gene rule, using the
            reaction ID as a pseudo-enzyme node.  Default: ``True``.

    Yields:
        :class:`Interaction` records.  For metabolite → enzyme edges,
        *source_type* is ``'small_molecule'`` and *target_type* is
        ``'protein'``.  For enzyme → metabolite edges the roles are
        swapped.  ``id_type_a`` / ``id_type_b`` are ``'metatlas'`` for
        metabolites, ``'ensembl'`` for enzymes, and ``'reaction_id'``
        for orphan reaction pseudo-enzyme nodes.
    """

    from pypath.inputs.metatlas._gem import (
        metatlas_gem_network,
        metatlas_gem_yaml_reactions,
        metatlas_gem_yaml_metabolites,
    )
    from pypath.inputs.metatlas._records import GemInteraction

    gems = [gem] if isinstance(gem, str) else list(gem)

    # Two-pass: collect everything first, then filter by metabolite degree.
    raw: list[tuple] = []  # (GemInteraction, gem_name)

    # Pre-build transport reaction ID sets (reuses already-cached YAML).
    transport_ids: dict[str, frozenset[str]] = {
        gem_name: _transport_reaction_ids(gem_name)
        for gem_name in gems
    }

    for gem_name in gems:

        for rec in metatlas_gem_network(gem=gem_name):

            if not include_reverse and rec.reverse:
                continue

            raw.append((rec, gem_name))

        if include_orphans:

            # Build metabolite → compartment lookup (YAML already cached).
            met_comp = {
                met.id: met.compartment
                for met in metatlas_gem_yaml_metabolites(gem=gem_name)
            }

            for rxn in metatlas_gem_yaml_reactions(gem=gem_name):

                if rxn.gene_reaction_rule and rxn.gene_reaction_rule.strip():
                    continue  # handled by metatlas_gem_network above

                rxn_id = rxn.id
                lb = rxn.lower_bound
                ub = rxn.upper_bound
                direction = 1 if lb + ub >= 0 else -1
                reversible = lb < 0 < ub

                mets = rxn.metabolites
                if isinstance(mets, list):
                    mets = dict(mets)

                reactants = [m for m, c in mets.items() if c * direction < 0]
                products = [m for m, c in mets.items() if c * direction > 0]

                for met_id in reactants:
                    raw.append((GemInteraction(
                        source=met_id, target=rxn_id,
                        source_type='metabolite', target_type='reaction',
                        source_compartment=met_comp.get(met_id, ''),
                        target_compartment='',
                        reaction_id=rxn_id, reverse=False,
                    ), gem_name))

                for met_id in products:
                    raw.append((GemInteraction(
                        source=rxn_id, target=met_id,
                        source_type='reaction', target_type='metabolite',
                        source_compartment='',
                        target_compartment=met_comp.get(met_id, ''),
                        reaction_id=rxn_id, reverse=False,
                    ), gem_name))

                if include_reverse and reversible:

                    for met_id in products:
                        raw.append((GemInteraction(
                            source=met_id, target=rxn_id,
                            source_type='metabolite', target_type='reaction',
                            source_compartment=met_comp.get(met_id, ''),
                            target_compartment='',
                            reaction_id=rxn_id, reverse=True,
                        ), gem_name))

                    for met_id in reactants:
                        raw.append((GemInteraction(
                            source=rxn_id, target=met_id,
                            source_type='reaction', target_type='metabolite',
                            source_compartment='',
                            target_compartment=met_comp.get(met_id, ''),
                            reaction_id=rxn_id, reverse=True,
                        ), gem_name))

    # Count each metabolite's total degree across all collected edges.
    metab_degree: Counter = Counter()

    for rec, _ in raw:

        if rec.source_type == 'metabolite':
            base_id = _strip_compartment(rec.source, rec.source_compartment)
            metab_degree[base_id] += 1

        if rec.target_type == 'metabolite':
            base_id = _strip_compartment(rec.target, rec.target_compartment)
            metab_degree[base_id] += 1

    # Yield interactions, skipping high-degree (cofactor) metabolites.
    for rec, gem_name in raw:

        is_orphan = rec.source_type == 'reaction' or rec.target_type == 'reaction'

        if rec.source_type == 'metabolite':
            met_id = _strip_compartment(rec.source, rec.source_compartment)
            enzyme_id = rec.target
            compartment = rec.source_compartment
            source = met_id
            target = enzyme_id
            source_type = 'small_molecule'
            target_type = 'protein'
            id_type_a = 'metatlas'
            id_type_b = 'reaction_id' if is_orphan else 'ensembl'

        else:
            enzyme_id = rec.source
            met_id = _strip_compartment(rec.target, rec.target_compartment)
            compartment = rec.target_compartment
            source = enzyme_id
            target = met_id
            source_type = 'protein'
            target_type = 'small_molecule'
            id_type_a = 'reaction_id' if is_orphan else 'ensembl'
            id_type_b = 'metatlas'

        if metab_degree[met_id] > metab_max_degree:
            continue

        locations = (compartment,) if compartment else ()

        is_complex = '_' in enzyme_id and not is_orphan

        resource_prefix = (
            'GEM_transporter'
            if rec.reaction_id in transport_ids[gem_name]
            else 'GEM'
        )

        attrs = {
            'reverse': rec.reverse,
            'reaction_id': rec.reaction_id,
            'enzyme_complex': is_complex,
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
            interaction_type='catalysis',
            resource=f'{resource_prefix}:{gem_name}',
            mor=0,
            locations=locations,
            attrs=attrs,
        )
