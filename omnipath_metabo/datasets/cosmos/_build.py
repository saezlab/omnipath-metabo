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

def _collect_entities(
    prov: pd.DataFrame,
    entity_type: str,
    make_record,
    exclude_reaction_id: bool = False,
) -> list:
    """
    Build provenance records for one entity type from both source and target columns.

    Args:
        prov: Merged DataFrame with translated IDs (``source``, ``target``),
            original IDs (``_orig_source``, ``_orig_target``,
            ``_orig_id_type_a``, ``_orig_id_type_b``), plus ``source_type``,
            ``target_type``, ``resource``.
        entity_type: ``'small_molecule'`` or ``'protein'``.
        make_record: Callable ``(canonical_id, orig_id, id_type, resource)``
            that returns the appropriate record namedtuple.
        exclude_reaction_id: If ``True``, skip rows where ``id_type``
            is ``'reaction_id'`` (used for orphan pseudo-enzymes).

    Returns:
        Deduplicated list of records.
    """
    records = []
    seen: set = set()

    for entity_col, id_col, orig_col, id_type_col in [
        ('source_type', 'source', '_orig_source', '_orig_id_type_a'),
        ('target_type', 'target', '_orig_target', '_orig_id_type_b'),
    ]:
        mask = prov[entity_col] == entity_type
        if exclude_reaction_id:
            mask = mask & (prov[id_type_col] != 'reaction_id')
        for _, row in prov[mask].iterrows():
            canonical = row[id_col]
            orig = row[orig_col]
            id_type = row[id_type_col]
            resource = row['resource']
            key = (canonical, orig, id_type, resource)
            if key not in seen:
                seen.add(key)
                records.append(make_record(canonical, orig, id_type, resource))

    return records


def _collect_metabolites(prov: pd.DataFrame) -> list:
    """Build :class:`~._record.CosmosMetabolite` provenance records."""
    return _collect_entities(
        prov, 'small_molecule',
        lambda c, o, t, r: CosmosMetabolite(chebi=c, original_id=o, id_type=t, resource=r),
    )


def _collect_proteins(prov: pd.DataFrame) -> list:
    """
    Build :class:`~._record.CosmosProtein` provenance records.

    Orphan reaction-ID pseudo-enzymes (``id_type == 'reaction_id'``) are
    excluded — they are not real proteins.
    """
    return _collect_entities(
        prov, 'protein',
        lambda c, o, t, r: CosmosProtein(uniprot=c, original_id=o, id_type=t, resource=r),
        exclude_reaction_id=True,
    )


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


def _report_resource_overlaps(
    bundle: CosmosBundle,
    category: str,
    translate: bool = True,
) -> None:
    """
    Log unique edge counts and pairwise overlaps between resources.

    Diagnostic only — no edges are removed.  Called at the end of each
    ``build_*()`` function when ``translate_ids=True``.

    For the transporter and receptor categories the edge identity key is
    ``(source, target, compartment)`` with one entry per location
    (multi-compartment rows are exploded), using only
    ``source_type == 'small_molecule'`` rows.  For transporters this
    also avoids double-counting bidirectional GEM / Recon3D edges.

    For all other categories the key is ``(source, target)``.

    Args:
        bundle: Finalized :class:`CosmosBundle` after filtering.
        category: Human-readable name used in log messages
            (e.g. ``'transporter'``, ``'receptor'``).
        translate: If ``False``, skip reporting — IDs are not yet
            in canonical form so counts would be misleading.
    """
    if not translate or not bundle.network:
        return

    def _label(resource: str) -> str:
        if resource.startswith('GEM_transporter:'):
            return 'GEM_transporter'
        if resource.startswith('GEM:'):
            return 'GEM'
        return resource

    resource_edges: dict[str, set] = {}

    for row in bundle.network:
        label = _label(row.resource)
        if category in ('transporter', 'receptor'):
            if row.source_type != 'small_molecule':
                continue
            locs = row.locations if row.locations else (None,)
            keys = {(row.source, row.target, loc) for loc in locs}
        else:
            keys = {(row.source, row.target)}

        resource_edges.setdefault(label, set()).update(keys)

    if not resource_edges:
        return

    labels = sorted(resource_edges)
    lines = [f'[COSMOS] {category} edge counts:']

    for label in labels:
        lines.append(f'  {label:<20}: {len(resource_edges[label]):,} unique')

    for i, la in enumerate(labels):
        for lb in labels[i + 1:]:
            n = len(resource_edges[la] & resource_edges[lb])
            lines.append(f'  {la} \u2229 {lb}: {n:,}')

    _log.info('\n'.join(lines))


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
# Location enrichment
# ---------------------------------------------------------------------------

def _enrich_stitch_locations(df: pd.DataFrame, organism: int) -> pd.DataFrame:
    """
    Add subcellular location abbreviations to STITCH rows post-translation.

    Called after ``translate_pkn()`` so protein IDs are already UniProt ACs
    — no additional ID mapping is needed.  Location data is loaded once and
    applied via dict lookup, keeping this step fast.

    Rows where no location can be resolved retain an empty tuple (consistent
    with other resources that have no location information).

    Args:
        df: Translated PKN DataFrame containing a ``'resource'`` column.
        organism: NCBI taxonomy ID used to filter UniProt location data.

    Returns:
        DataFrame with ``'locations'`` column updated for STITCH rows.
    """
    from .location import resolve_protein_locations, tcdb_locations, uniprot_locations

    stitch_mask = df['resource'] == 'STITCH'

    if not stitch_mask.any():
        return df

    _log.info('[COSMOS] Enriching STITCH locations...')

    all_locations = uniprot_locations(organism=organism, reviewed=True)
    loc_map = tcdb_locations()

    def _loc(uniprot: str) -> tuple:
        abbr = resolve_protein_locations(uniprot, all_locations, loc_map)
        return tuple(sorted(abbr)) if abbr else ()

    df = df.copy()
    df.loc[stitch_mask, 'locations'] = (
        df.loc[stitch_mask, 'target'].apply(_loc)
    )

    return df


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build(
    *args: dict | 'Path' | str,
    row_filter=None,
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
        row_filter:
            Optional callable ``(Interaction) → bool``.  When provided,
            applied *before* ID translation so that rows known to be
            unneeded are dropped early and never translated.  This avoids
            paying the cost of translating metabolic GEM edges when only
            transporter edges are requested, for example.
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

    if row_filter is not None:
        n_before_filter = len(df)
        mask = df.apply(lambda r: row_filter(Interaction(*r)), axis=1)
        df = df[mask].reset_index(drop=True)
        _log.info(
            '[COSMOS] Pre-translation filter: %d/%d edges retained.',
            len(df),
            n_before_filter,
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
        df = _enrich_stitch_locations(df, organism)

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


def build_transporters(*args, cell_surface_only: bool = False, **kwargs) -> CosmosBundle:
    """
    Build the transporter subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only
    transporter-relevant resources (TCDB, SLC, GEM, Recon3D) and
    filters to keep only transporter interactions.

    STITCH is excluded because it is a general chemical-protein interaction
    database (last updated 2015) that does not annotate transport specifically.
    Classifying STITCH interactions as transport by checking whether the
    protein appears in TCDB is methodologically unsound: a protein can be a
    transporter for some substrates while acting as a kinase, receptor, or
    enzyme for others. Dedicated transporter databases (TCDB, SLC, GEM,
    Recon3D) provide direct, mechanistically-annotated transport interactions
    and are actively maintained. See ADR 0002 in saezverse.

    The filter is applied *before* ID translation so that metabolic GEM edges
    (``resource='GEM:<gem>'``) are discarded early and never translated,
    avoiding redundant computation.

    Filter predicate:
        - ``interaction_type in ('transport', 'transporter')``
        - ``resource.startswith('GEM_transporter')``

    Args:
        *args: Passed through to :func:`build`.
        **kwargs: Passed through to :func:`build`.  ``brenda`` and
            ``stitch`` are disabled unless explicitly re-enabled.
            ``mrclinksdb`` is enabled: its transport-classified records
            (``interaction_type='transport'``) are included.

    Returns:
        :class:`CosmosBundle` containing only transporter interactions,
        with provenance filtered to the surviving edges.
    """
    def _is_transport(row: Interaction) -> bool:
        is_transport_type = (
            row.interaction_type in ('transport', 'transporter') or
            row.resource.startswith('GEM_transporter')
        )
        if not is_transport_type:
            return False
        if cell_surface_only:
            return 'e' in row.locations
        return True

    kwargs.setdefault('brenda', False)
    kwargs.setdefault('stitch', False)
    bundle = build(*args, row_filter=_is_transport, **kwargs)
    _report_resource_overlaps(bundle, 'transporter', kwargs.get('translate_ids', True))
    return bundle


def build_receptors(*args, cell_surface_only: bool = False, **kwargs) -> CosmosBundle:
    """
    Build the receptor subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only
    receptor-relevant resources (MRCLinksDB, STITCH) and post-filters
    to keep only receptor/ligand interactions.

    Unlike :func:`build_transporters`, the ``cell_surface_only`` filter is
    applied *after* ID translation (via :func:`_filter_bundle`) rather than
    as a pre-translation ``row_filter``.  This is necessary because STITCH
    proteins have no location data at yield time — locations are assigned
    post-translation by :func:`_enrich_stitch_locations` inside
    :func:`build`.  MRCLinksDB locations are set at yield time and would
    support pre-filtering, but using a single post-translation pass keeps
    the implementation uniform.

    Post-filter predicate:
        - ``interaction_type == 'ligand_receptor'``
        - If ``cell_surface_only=True``: also ``'e' in row.locations``

    Args:
        *args: Passed through to :func:`build`.
        cell_surface_only:
            If ``True``, retain only interactions where the receptor
            protein is annotated to the plasma membrane / cell surface
            (``'e'`` in ``locations``).  Useful for cell-cell
            communication models (COSMOS intercellular layer, NicheNet,
            etc.) where only surface-exposed receptors are relevant.
            Locations are derived from UniProt subcellular location
            annotations mapped via the TCDB location table.
        **kwargs: Passed through to :func:`build`.  ``tcdb``, ``slc``,
            ``brenda``, ``gem``, and ``recon3d`` are disabled unless
            explicitly re-enabled.

    Returns:
        :class:`CosmosBundle` containing only receptor interactions,
        with provenance filtered to the surviving edges.
    """
    def _is_receptor(row: Interaction) -> bool:
        is_rec = (
            row.interaction_type == 'ligand_receptor' or
            (row.resource == 'STITCH' and row.interaction_type == 'receptor')
        )
        if not is_rec:
            return False
        if cell_surface_only:
            return 'e' in row.locations
        return True

    kwargs.setdefault('tcdb', False)
    kwargs.setdefault('slc', False)
    kwargs.setdefault('brenda', False)
    kwargs.setdefault('gem', False)
    kwargs.setdefault('recon3d', False)
    bundle = build(*args, **kwargs)
    bundle = _filter_bundle(bundle, _is_receptor)
    _report_resource_overlaps(bundle, 'receptor', kwargs.get('translate_ids', True))
    return bundle


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
    bundle = _filter_bundle(bundle, lambda row: (
        row.interaction_type == 'allosteric_regulation' or
        (row.resource == 'STITCH' and row.interaction_type == 'other')
    ))
    _report_resource_overlaps(bundle, 'allosteric', kwargs.get('translate_ids', True))
    return bundle


def build_enzyme_metabolite(*args, **kwargs) -> CosmosBundle:
    """
    Build the enzyme-metabolite (metabolic) subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only GEM
    resources and post-filters to keep only stoichiometric
    enzyme-metabolite interactions from genome-scale metabolic models.

    Corresponds to the *Enzyme-metabolite* category in the COSMOS PKN
    planning document: direct metabolic reactions where enzymes act on
    substrates and products, as opposed to allosteric regulation.

    Pre-translation filter predicate:
        - ``resource.startswith('GEM:')`` (metabolic GEM edges only,
          distinct from ``'GEM_transporter:'`` transport edges)

    The filter is applied *before* ID translation so that transport GEM
    edges (``resource='GEM_transporter:<gem>'``) are discarded early and
    never translated — avoiding spurious drop-rate warnings.

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
    bundle = build(*args, row_filter=lambda row: row.resource.startswith('GEM:'), **kwargs)
    _report_resource_overlaps(bundle, 'enzyme-metabolite', kwargs.get('translate_ids', True))
    return bundle
