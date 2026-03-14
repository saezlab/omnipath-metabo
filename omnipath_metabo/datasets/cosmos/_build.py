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
COSMOS PKN builder.

This module provides the main entry point for building the COSMOS
prior-knowledge network by combining data from multiple sources.

Resource categories:

1. Metabolite-protein interactions (current)
2. GEMs (genome-scale metabolic models)
3. PPIs (protein-protein interactions)
4. GRNs (gene regulatory networks)
5. Kinase-substrate networks

ID unification (when ``translate_ids`` is ``True``):
    - Metabolite source IDs → ChEBI (via UniChem or PubChem REST API)
    - Protein target IDs → UniProt accessions (via pypath BioMart)
"""

from __future__ import annotations

__all__ = ['build', 'build_transporters', 'build_receptors', 'build_allosteric', 'build_enzyme_metabolite']

import logging
from itertools import chain
from typing import TYPE_CHECKING

import pandas as pd

from ._bundle import CosmosBundle
from ._config import config
from ._record import CosmosMetabolite, CosmosProtein, CosmosReaction, Interaction

_log = logging.getLogger(__name__)

try:
    from tqdm import tqdm as _tqdm
except ImportError:
    _tqdm = None


def _progress(iterable, desc: str, **kwargs):
    """Wrap *iterable* with tqdm if available, otherwise return it unchanged."""
    if _tqdm is not None:
        return _tqdm(iterable, desc=desc, **kwargs)
    return iterable
from .resources import (
    brenda_regulations,
    gem_interactions,
    mrclinksdb_interactions,
    recon3d_transporter_interactions,
    slc_interactions,
    stitch_interactions,
    tcdb_interactions,
)

if TYPE_CHECKING:
    from pathlib import Path


PROCESSORS = {
    'stitch': stitch_interactions,
    'tcdb': tcdb_interactions,
    'slc': slc_interactions,
    'brenda': brenda_regulations,
    'mrclinksdb': mrclinksdb_interactions,
    'gem': gem_interactions,
    'recon3d': recon3d_transporter_interactions,
}


# ---------------------------------------------------------------------------
# Provenance helpers
# ---------------------------------------------------------------------------

def _collect_metabolites(prov: pd.DataFrame) -> list:
    """
    Build :class:`~._record.CosmosMetabolite` provenance records.

    Args:
        prov: Merged DataFrame with columns for translated IDs (``source``,
            ``target``) and original IDs (``_orig_source``, ``_orig_target``,
            ``_orig_id_type_a``, ``_orig_id_type_b``), plus ``source_type``,
            ``target_type``, ``resource``.

    Returns:
        Deduplicated list of :class:`CosmosMetabolite` records.
    """
    records = []
    seen: set = set()

    for entity_col, chebi_col, orig_col, id_type_col in [
        ('source_type', 'source', '_orig_source', '_orig_id_type_a'),
        ('target_type', 'target', '_orig_target', '_orig_id_type_b'),
    ]:
        mask = prov[entity_col] == 'small_molecule'
        for _, row in prov[mask].iterrows():
            chebi = row[chebi_col]
            orig = row[orig_col]
            id_type = row[id_type_col]
            resource = row['resource']
            key = (chebi, orig, id_type, resource)
            if key not in seen:
                seen.add(key)
                records.append(CosmosMetabolite(
                    chebi=chebi,
                    original_id=orig,
                    id_type=id_type,
                    resource=resource,
                ))

    return records


def _collect_proteins(prov: pd.DataFrame) -> list:
    """
    Build :class:`~._record.CosmosProtein` provenance records.

    Orphan reaction-ID pseudo-enzymes (``id_type == 'reaction_id'``) are
    excluded — they are not real proteins.

    Args:
        prov: Merged DataFrame (same structure as for
            :func:`_collect_metabolites`).

    Returns:
        Deduplicated list of :class:`CosmosProtein` records.
    """
    records = []
    seen: set = set()

    for entity_col, uniprot_col, orig_col, id_type_col in [
        ('source_type', 'source', '_orig_source', '_orig_id_type_a'),
        ('target_type', 'target', '_orig_target', '_orig_id_type_b'),
    ]:
        mask = (
            (prov[entity_col] == 'protein') &
            (prov[id_type_col] != 'reaction_id')
        )
        for _, row in prov[mask].iterrows():
            uniprot = row[uniprot_col]
            orig = row[orig_col]
            id_type = row[id_type_col]
            resource = row['resource']
            key = (uniprot, orig, id_type, resource)
            if key not in seen:
                seen.add(key)
                records.append(CosmosProtein(
                    uniprot=uniprot,
                    original_id=orig,
                    id_type=id_type,
                    resource=resource,
                ))

    return records


def _collect_reactions(df: pd.DataFrame) -> list:
    """
    Build :class:`~._record.CosmosReaction` metadata records from GEM rows.

    Groups translated GEM / Recon3D rows by ``(reaction_id, gem_name)`` and
    collects the set of UniProt ACs (genes) and ChEBI IDs (metabolites) that
    participate in each reaction.

    Args:
        df: Translated PKN DataFrame (may include ``_row_id`` column).

    Returns:
        List of :class:`CosmosReaction` records, one per unique
        ``(reaction_id, gem)`` pair.
    """
    gem_mask = (
        df['resource'].str.startswith('GEM') | (df['resource'] == 'Recon3D')
    )
    if not gem_mask.any():
        return []

    seen: dict = {}  # (reaction_id, gem) → {'genes': set, 'metabolites': set}

    for _, row in df[gem_mask].iterrows():
        attrs = row['attrs'] if isinstance(row['attrs'], dict) else {}
        reaction_id = attrs.get('reaction_id', '')
        if not reaction_id:
            continue

        resource = row['resource']
        gem = resource.split(':', 1)[1] if ':' in resource else resource

        key = (reaction_id, gem)
        if key not in seen:
            seen[key] = {'genes': set(), 'metabolites': set()}

        if row['source_type'] == 'protein' and row['id_type_a'] != 'reaction_id':
            seen[key]['genes'].add(row['source'])
        if row['target_type'] == 'protein' and row['id_type_b'] != 'reaction_id':
            seen[key]['genes'].add(row['target'])
        if row['source_type'] == 'small_molecule':
            seen[key]['metabolites'].add(row['source'])
        if row['target_type'] == 'small_molecule':
            seen[key]['metabolites'].add(row['target'])

    return [
        CosmosReaction(
            reaction_id=rxn_id,
            gem=gem,
            genes=tuple(sorted(gm['genes'])),
            metabolites=tuple(sorted(gm['metabolites'])),
        )
        for (rxn_id, gem), gm in seen.items()
    ]


def _filter_bundle(bundle: CosmosBundle, predicate) -> CosmosBundle:
    """
    Return a copy of *bundle* keeping only network rows that satisfy *predicate*.

    Provenance lists (``metabolites``, ``proteins``, ``reactions``) are
    filtered to include only entries referenced by the surviving network rows.

    Args:
        bundle: Source :class:`CosmosBundle`.
        predicate: Callable ``(Interaction) → bool``.

    Returns:
        Filtered :class:`CosmosBundle`.
    """
    filtered = [row for row in bundle.network if predicate(row)]

    met_ids: set = set()
    prot_ids: set = set()
    rxn_ids: set = set()

    for row in filtered:
        if row.source_type == 'small_molecule':
            met_ids.add(row.source)
        elif row.source_type == 'protein' and row.id_type_a != 'reaction_id':
            prot_ids.add(row.source)
        if row.target_type == 'small_molecule':
            met_ids.add(row.target)
        elif row.target_type == 'protein' and row.id_type_b != 'reaction_id':
            prot_ids.add(row.target)
        if isinstance(row.attrs, dict) and row.attrs.get('reaction_id'):
            rxn_ids.add(row.attrs['reaction_id'])

    return CosmosBundle(
        network=filtered,
        metabolites=[m for m in bundle.metabolites if m.chebi in met_ids],
        proteins=[p for p in bundle.proteins if p.uniprot in prot_ids],
        reactions=[r for r in bundle.reactions if r.reaction_id in rxn_ids],
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build(
    *args: dict | 'Path' | str,
    **kwargs,
) -> CosmosBundle:
    """
    Build the COSMOS prior-knowledge network from multiple sources.

    Calls each requested resource processor, collects the yielded
    :class:`Interaction` records, optionally translates IDs to canonical
    form (ChEBI / UniProt), and returns a :class:`CosmosBundle`.

    When ``translate_ids`` is ``True`` (default), the bundle's
    ``metabolites``, ``proteins``, and ``reactions`` lists are populated
    with provenance records linking original source IDs to their canonical
    counterparts.  The ``network`` list holds :class:`Interaction`
    namedtuples with translated IDs.

    The ``resources`` dict in the config controls both which resources are
    active and their parameters.  The top-level ``organism`` is injected
    into each resource unless overridden.

    Args:
        *args:
            Configuration overrides as dicts or YAML file paths.
        **kwargs:
            Top-level config keys.  Resource names are accepted as
            shorthand, e.g.::

                build(stitch={'score_threshold': 400})
                build(tcdb=False, slc=False)
                build(organism=10090)

    Returns:
        :class:`CosmosBundle` with all interactions in ``network``.
        When ``translate_ids`` is ``True``, provenance records are
        populated in ``metabolites``, ``proteins``, and ``reactions``.
    """

    cfg = config(*args, **kwargs)
    organism = cfg.get('organism', 9606)
    resources = cfg.get('resources', {})
    translate = cfg.get('translate_ids', True)

    generators = []
    active_names = []

    for name, params in resources.items():

        if params is False:
            continue

        if name not in PROCESSORS:
            raise ValueError(
                f'Unknown resource: {name!r}. '
                f'Available: {list(PROCESSORS)}'
            )

        resource_args = params if isinstance(params, dict) else {}
        resource_args.setdefault('organism', organism)
        active_names.append(name)
        generators.append(PROCESSORS[name](**resource_args))

    for name in _progress(active_names, desc='[COSMOS] collecting resources'):
        _log.info('[COSMOS] Building %s...', name)

    df = pd.DataFrame(
        chain.from_iterable(generators),
        columns=Interaction._fields,
    )

    _log.info(
        '[COSMOS] Collected %d edges from %d resources.',
        len(df),
        len(active_names),
    )

    metabolites: list = []
    proteins: list = []
    reactions: list = []

    if translate:
        from ._translate import translate_pkn

        n_before = len(df)
        df['_row_id'] = range(len(df))
        df_raw = df[['_row_id', 'source', 'target', 'id_type_a', 'id_type_b']].copy()

        _log.info('[COSMOS] Translating IDs...')
        df = translate_pkn(df, organism=organism)
        _log.info(
            '[COSMOS] Translation complete: %d/%d edges retained.',
            len(df),
            n_before,
        )

        if cfg.get('apply_blacklist', True):
            from ._blacklist import apply_blacklist
            df = apply_blacklist(df)

        # Merge to recover original IDs for provenance tracking.
        prov = df[
            ['_row_id', 'source', 'target', 'source_type', 'target_type', 'resource']
        ].merge(
            df_raw.rename(columns={
                'source': '_orig_source',
                'target': '_orig_target',
                'id_type_a': '_orig_id_type_a',
                'id_type_b': '_orig_id_type_b',
            }),
            on='_row_id',
        )

        metabolites = _collect_metabolites(prov)
        proteins = _collect_proteins(prov)
        reactions = _collect_reactions(df)
        df = df.drop(columns=['_row_id'])

    else:
        if cfg.get('apply_blacklist', True):
            from ._blacklist import apply_blacklist
            df = apply_blacklist(df)

    network = [Interaction(*row) for row in df.itertuples(index=False, name=None)]

    return CosmosBundle(
        network=network,
        metabolites=metabolites,
        proteins=proteins,
        reactions=reactions,
    )


def build_transporters(*args, **kwargs) -> CosmosBundle:
    """
    Build the transporter subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only
    transporter-relevant resources (TCDB, SLC, GEM, Recon3D) and
    post-filters to keep only transporter interactions.

    STITCH is excluded because it is a general chemical-protein interaction
    database (last updated 2015) that does not annotate transport specifically.
    Classifying STITCH interactions as transport by checking whether the
    protein appears in TCDB is methodologically unsound: a protein can be a
    transporter for some substrates while acting as a kinase, receptor, or
    enzyme for others. Dedicated transporter databases (TCDB, SLC, GEM,
    Recon3D) provide direct, mechanistically-annotated transport interactions
    and are actively maintained. See ADR 0002 in saezverse.

    Post-filter predicate:
        - ``interaction_type == 'transport'``
        - ``resource.startswith('GEM_transporter')``

    Args:
        *args: Passed through to :func:`build`.
        **kwargs: Passed through to :func:`build`.  ``brenda``,
            ``mrclinksdb``, and ``stitch`` are disabled unless explicitly
            re-enabled.

    Returns:
        :class:`CosmosBundle` containing only transporter interactions,
        with provenance filtered to the surviving edges.
    """
    kwargs.setdefault('brenda', False)
    kwargs.setdefault('mrclinksdb', False)
    kwargs.setdefault('stitch', False)
    bundle = build(*args, **kwargs)
    return _filter_bundle(bundle, lambda row: (
        row.interaction_type == 'transport' or
        row.resource.startswith('GEM_transporter')
    ))


def build_receptors(*args, **kwargs) -> CosmosBundle:
    """
    Build the receptor subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only
    receptor-relevant resources (MRCLinksDB, STITCH) and post-filters
    to keep only receptor/ligand interactions.

    Post-filter predicate:
        - ``interaction_type == 'ligand_receptor'``
        - STITCH rows where ``interaction_type == 'receptor'``

    Args:
        *args: Passed through to :func:`build`.
        **kwargs: Passed through to :func:`build`.  ``tcdb``, ``slc``,
            ``brenda``, ``gem``, and ``recon3d`` are disabled unless
            explicitly re-enabled.

    Returns:
        :class:`CosmosBundle` containing only receptor interactions,
        with provenance filtered to the surviving edges.
    """
    kwargs.setdefault('tcdb', False)
    kwargs.setdefault('slc', False)
    kwargs.setdefault('brenda', False)
    kwargs.setdefault('gem', False)
    kwargs.setdefault('recon3d', False)
    bundle = build(*args, **kwargs)
    return _filter_bundle(bundle, lambda row: (
        row.interaction_type == 'ligand_receptor' or
        (row.resource == 'STITCH' and row.interaction_type == 'receptor')
    ))


def build_allosteric(*args, **kwargs) -> CosmosBundle:
    """
    Build the allosteric-regulation subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only
    allosteric-relevant resources (BRENDA, STITCH) and post-filters
    to keep only allosteric interactions.

    Corresponds to the *Metabolite-protein interaction* category in the
    COSMOS PKN planning document: small molecules that activate or
    inhibit proteins through allosteric binding, distinct from
    stoichiometric enzymatic metabolism.

    Post-filter predicate:
        - ``interaction_type == 'allosteric_regulation'`` (BRENDA)
        - STITCH rows where ``interaction_type == 'other'``

    Args:
        *args: Passed through to :func:`build`.
        **kwargs: Passed through to :func:`build`.  ``tcdb``, ``slc``,
            ``mrclinksdb``, ``gem``, and ``recon3d`` are disabled unless
            explicitly re-enabled.

    Returns:
        :class:`CosmosBundle` containing only allosteric regulation
        interactions, with provenance filtered to the surviving edges.
    """
    kwargs.setdefault('tcdb', False)
    kwargs.setdefault('slc', False)
    kwargs.setdefault('mrclinksdb', False)
    kwargs.setdefault('gem', False)
    kwargs.setdefault('recon3d', False)
    bundle = build(*args, **kwargs)
    return _filter_bundle(bundle, lambda row: (
        row.interaction_type == 'allosteric_regulation' or
        (row.resource == 'STITCH' and row.interaction_type == 'other')
    ))


def build_enzyme_metabolite(*args, **kwargs) -> CosmosBundle:
    """
    Build the enzyme-metabolite (metabolic) subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only GEM
    resources and post-filters to keep only stoichiometric
    enzyme-metabolite interactions from genome-scale metabolic models.

    Corresponds to the *Enzyme-metabolite* category in the COSMOS PKN
    planning document: direct metabolic reactions where enzymes act on
    substrates and products, as opposed to allosteric regulation.

    Post-filter predicate:
        - ``resource.startswith('GEM:')`` (metabolic GEM edges,
          distinct from ``'GEM_transporter:'`` transport edges)

    Note:
        ``'GEM_transporter:...'`` resources are excluded because
        ``'GEM_transporter:...'.startswith('GEM:')`` is ``False``.
        For allosteric regulation (BRENDA, STITCH), use
        :func:`build_allosteric`.

    Args:
        *args: Passed through to :func:`build`.
        **kwargs: Passed through to :func:`build`.  ``tcdb``, ``slc``,
            ``brenda``, ``mrclinksdb``, ``recon3d``, and ``stitch`` are
            disabled unless explicitly re-enabled.

    Returns:
        :class:`CosmosBundle` containing only stoichiometric
        enzyme-metabolite interactions from GEMs, with provenance
        filtered to the surviving edges.
    """
    kwargs.setdefault('tcdb', False)
    kwargs.setdefault('slc', False)
    kwargs.setdefault('brenda', False)
    kwargs.setdefault('mrclinksdb', False)
    kwargs.setdefault('recon3d', False)
    kwargs.setdefault('stitch', False)
    bundle = build(*args, **kwargs)
    return _filter_bundle(bundle, lambda row: row.resource.startswith('GEM:'))
