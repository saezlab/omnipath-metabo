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

import re
from collections import defaultdict
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


def _crossing_metabolites(edges: list) -> set[str]:
    """
    Return base metabolite IDs that cross a compartment boundary.

    Given all GemInteraction edges for a single reaction, finds metabolites
    whose base ID (compartment suffix stripped) appears on both the reactant
    side (``source_type == 'metabolite'``) and the product side
    (``target_type == 'metabolite'``) with *different* compartment codes.

    Cofactors consumed and regenerated within the same compartment (e.g.
    ``atp_c → adp_c`` within the cytoplasm) never satisfy this condition
    and are excluded.  ATP as a cofactor in an ABC transporter also fails
    because ATP and ADP have different base IDs — there is no base ID that
    appears on *both* sides with a compartment change.

    Args:
        edges: All GemInteraction records for a single reaction (may span
            multiple enzymes and/or forward + reverse edges).

    Returns:
        Set of base metabolite IDs (compartment suffix stripped) that
        physically cross a membrane in this reaction.
    """

    reactant_comps: dict[str, str] = {}
    product_comps: dict[str, str] = {}

    for rec in edges:

        if rec.source_type == 'metabolite':
            base = _strip_compartment(rec.source, rec.source_compartment)
            reactant_comps[base] = rec.source_compartment

        if rec.target_type == 'metabolite':
            base = _strip_compartment(rec.target, rec.target_compartment)
            product_comps[base] = rec.target_compartment

    return {
        base
        for base, comp in reactant_comps.items()
        if base in product_comps and product_comps[base] != comp
    }


_ENSG_RE = re.compile(r'^ENSG\d{11}$')


def _detect_gene_id_type(gem_name: str) -> str:
    """
    Detect the gene identifier type used in a GEM's gene-reaction rules.

    Human-GEM uses Ensembl gene IDs (``ENSG00000000000``).  Mouse-GEM and
    other species GEMs from MetAtlas use gene symbols (e.g. ``Adh1``).
    The type is determined by checking the first non-empty gene token
    against the ENSG pattern.

    Args:
        gem_name: GEM name (e.g. ``'Human-GEM'``, ``'Mouse-GEM'``).

    Returns:
        ``'ensembl'`` if genes look like ENSG IDs, ``'genesymbol'``
        otherwise.
    """

    from pypath.inputs.metatlas._gem import metatlas_gem_yaml_reactions

    for rxn in metatlas_gem_yaml_reactions(gem=gem_name):
        rule = rxn.gene_reaction_rule or ''
        # Extract the first gene token from the rule string.
        token = re.split(r'[\s()]+', rule.strip())[0]
        if token and token.lower() not in ('or', 'and'):
            return 'ensembl' if _ENSG_RE.match(token) else 'genesymbol'

    return 'ensembl'  # safe default


def gem_interactions(
    gem: str | list[str] = 'Human-GEM',
    organism: int = 9606,
    include_reverse: bool = True,
    include_orphans: bool = True,
    cell_surface_only: bool = False,
) -> Generator[Interaction, None, None]:
    """
    Yield GEM metabolite-enzyme interactions as uniform records.

    Reactions are decomposed into binary edges:

    - *metabolite → enzyme*: reactant consumed by the enzyme.
    - *enzyme → metabolite*: product released by the enzyme.

    Reversible reactions additionally produce a reversed pair with
    ``attrs['reverse'] = True``.

    Transport reactions (subsystem ``'Transport reactions'``) are filtered
    by compartment-crossing: only metabolites whose base ID appears on
    both the reactant and product side of the reaction with *different*
    compartment codes are kept.  This mechanistically excludes cofactors
    (e.g. ATP hydrolysed and regenerated within the same compartment) while
    preserving genuine transport substrates that physically cross a membrane.

    Metabolic reactions (all other subsystems) are not filtered: every
    metabolite in the reaction generates edges.

    Orphan reactions (no gene rule) are retained when *include_orphans*
    is ``True``.  The reaction ID is used as a pseudo-enzyme node with
    ``id_type = 'reaction_id'`` and ``attrs['orphan'] = True``.  Orphan
    transport reactions apply the same compartment-crossing filter as
    enzyme-annotated transport reactions.

    Args:
        gem:
            GEM name or list of GEM names to load and merge.  Names must
            match entries in the Metabolic Atlas validation index (e.g.
            ``'Human-GEM'``, ``'Mouse-GEM'``).
        organism:
            NCBI taxonomy ID.  Accepted for API consistency with other
            resources; GEMs are inherently organism-specific so this
            parameter is not used for filtering.
        include_reverse:
            If ``True``, include reversed edges for reversible reactions
            (``attrs['reverse'] = True``).  Default: ``True``.
        include_orphans:
            If ``True``, include reactions with no gene rule, using the
            reaction ID as a pseudo-enzyme node.  Default: ``True``.
        cell_surface_only:
            If ``True``, restrict transport reactions to those where at
            least one crossing metabolite touches the extracellular
            compartment (``'e'``).  Entire reaction groups are dropped
            before edge breakdown, so met→enzyme / enzyme→met pairs are
            always preserved.  Metabolic (non-transport) edges are
            unaffected.  Default: ``False``.

    Yields:
        :class:`Interaction` records.  For metabolite → enzyme edges,
        *source_type* is ``'small_molecule'`` and *target_type* is
        ``'protein'``.  For enzyme → metabolite edges the roles are
        swapped.  ``id_type_a`` / ``id_type_b`` are ``'metatlas'`` for
        metabolites, ``'ensembl'`` for enzymes, and ``'reaction_id'``
        for orphan reaction pseudo-enzyme nodes.

        Each yielded record carries ``attrs['gems']``: a sorted list of
        GEM names that contributed the edge.  When loading a single GEM
        this list has one element; when the same ``(source, target,
        reaction_id, reverse)`` key appears in multiple GEMs, the edge
        is deduplicated and ``attrs['gems']`` records all sources.

    Warning:
        Reverse edges are already included in the output when
        ``include_reverse=True`` (the default).  A downstream COSMOS
        formatter must NOT re-generate reverse edges.  Use
        ``attrs['reverse']`` to identify and label them, not to create
        new edges from them.
    """

    from pypath.inputs.metatlas._gem import (
        metatlas_gem_network,
        metatlas_gem_yaml_reactions,
        metatlas_gem_yaml_metabolites,
    )
    from pypath.inputs.metatlas._records import GemInteraction

    gems = [gem] if isinstance(gem, str) else list(gem)

    # Pre-build transport reaction ID sets (reuses already-cached YAML).
    transport_ids: dict[str, frozenset[str]] = {
        gem_name: _transport_reaction_ids(gem_name)
        for gem_name in gems
    }

    # Pre-detect gene ID type per GEM (samples first gene token from grRules).
    gene_id_types: dict[str, str] = {
        gem_name: _detect_gene_id_type(gem_name)
        for gem_name in gems
    }

    # Collect edges in two streams:
    #   - metabolic_raw: non-transport reactions, yielded without filtering.
    #   - transport_groups: transport reactions grouped by (reaction_id,
    #     gem_name) so compartment-crossing detection can be applied per
    #     reaction before yielding.
    metabolic_raw: list[tuple] = []
    transport_groups: dict[tuple, list] = defaultdict(list)

    for gem_name in gems:

        t_ids = transport_ids[gem_name]

        for rec in metatlas_gem_network(gem=gem_name):

            if not include_reverse and rec.reverse:
                continue

            if rec.reaction_id in t_ids:
                transport_groups[(rec.reaction_id, gem_name)].append(rec)
            else:
                metabolic_raw.append((rec, gem_name))

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

                reactant_ids = [m for m, c in mets.items() if c * direction < 0]
                product_ids = [m for m, c in mets.items() if c * direction > 0]

                is_transport = rxn_id in t_ids

                def _add_orphan(orphan_rec: GemInteraction) -> None:
                    if is_transport:
                        transport_groups[(rxn_id, gem_name)].append(orphan_rec)
                    else:
                        metabolic_raw.append((orphan_rec, gem_name))

                for met_id in reactant_ids:
                    _add_orphan(GemInteraction(
                        source=met_id, target=rxn_id,
                        source_type='metabolite', target_type='reaction',
                        source_compartment=met_comp.get(met_id, ''),
                        target_compartment='',
                        reaction_id=rxn_id, reverse=False,
                    ))

                for met_id in product_ids:
                    _add_orphan(GemInteraction(
                        source=rxn_id, target=met_id,
                        source_type='reaction', target_type='metabolite',
                        source_compartment='',
                        target_compartment=met_comp.get(met_id, ''),
                        reaction_id=rxn_id, reverse=False,
                    ))

                if include_reverse and reversible:

                    for met_id in product_ids:
                        _add_orphan(GemInteraction(
                            source=met_id, target=rxn_id,
                            source_type='metabolite', target_type='reaction',
                            source_compartment=met_comp.get(met_id, ''),
                            target_compartment='',
                            reaction_id=rxn_id, reverse=True,
                        ))

                    for met_id in reactant_ids:
                        _add_orphan(GemInteraction(
                            source=rxn_id, target=met_id,
                            source_type='reaction', target_type='metabolite',
                            source_compartment='',
                            target_compartment=met_comp.get(met_id, ''),
                            reaction_id=rxn_id, reverse=True,
                        ))

    # Apply compartment-crossing filter to all transport reactions
    # (both enzyme-annotated and orphan, since both land in transport_groups).
    transport_raw: list[tuple] = []

    for (rxn_id, gem_name), edges in transport_groups.items():
        crossing = _crossing_metabolites(edges)

        if cell_surface_only:
            crossing_comps: set[str] = set()
            for rec in edges:
                if rec.source_type == 'metabolite':
                    base = _strip_compartment(rec.source, rec.source_compartment)
                    if base in crossing:
                        crossing_comps.add(rec.source_compartment)
                if rec.target_type == 'metabolite':
                    base = _strip_compartment(rec.target, rec.target_compartment)
                    if base in crossing:
                        crossing_comps.add(rec.target_compartment)
            if 'e' not in crossing_comps:
                continue

        for rec in edges:
            if rec.source_type == 'metabolite':
                met_base = _strip_compartment(rec.source, rec.source_compartment)
            else:
                met_base = _strip_compartment(rec.target, rec.target_compartment)
            if met_base in crossing:
                transport_raw.append((rec, gem_name))

    # Build Interaction objects from both streams.
    # Collect first so deduplication across GEMs can attach attrs['gems'].
    all_interactions: list[tuple[Interaction, str]] = []

    for rec, gem_name in metabolic_raw + transport_raw:

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

    # Deduplicate across GEMs.
    #
    # Two edges from different GEMs are considered identical when their
    # (source, target, reaction_id, reverse) key matches.  For each unique
    # edge the first occurrence is kept; all contributing GEM names are
    # collected and stored in attrs['gems'] (always a sorted list, even for
    # single-GEM builds).
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
