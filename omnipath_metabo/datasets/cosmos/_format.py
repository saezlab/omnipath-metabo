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
COSMOS PKN formatter.

Applies COSMOS R package-compatible node-ID prefixes/suffixes to a
translated PKN DataFrame.  The formatter assumes the DataFrame has
already been passed through ``translate_pkn()`` so that all metabolite
IDs are ChEBI and all protein IDs are UniProt accessions.

Node-ID format
--------------
- Metabolite: ``Metab__CHEBI:xxxx_c``  (double underscore; single-
  underscore + compartment letter; no compartment suffix when unknown)
- Gene forward: ``Gene{N}__UniProtAC...``
- Gene reverse: ``Gene{N}__UniProtAC..._rev``
- Orphan reaction: ``Gene{N}__orphanReac<reaction_id>``

N is a sequential integer assigned per unique reaction within each
category (transporter / receptor / other).  The counter resets to 1
when moving to the next category.  For pre-expanded resources (GEM,
Recon3D) all rows sharing the same ``reaction_id`` receive the same N.
For every other resource each input row counts as one reaction.

Connector edges
---------------
After formatting, unique connector edges are appended::

    UniProtAC      →  Gene{N}__AC...              (one per unique formatted gene)
    UniProtAC      →  Gene{N}__AC..._rev          (transporter reverse genes)
    <reaction_id>  →  Gene{N}__orphanReac<id>     (orphan pseudo-enzyme)
    CHEBI:xxxx     →  Metab__CHEBI:xxxx_c         (one per unique formatted metabolite)

These edges allow downstream tools that key measurements by bare IDs
(e.g. transcriptomics → ENSG, metabolomics → CHEBI) to attach to the
formatted network nodes.
"""

from __future__ import annotations

__all__ = [
    'format_pkn',
    'format_transporters',
    'format_receptors',
    'format_allosteric',
    'format_enzyme_metabolite',
]

import logging

import pandas as pd

_log = logging.getLogger(__name__)

# Sources that already emit both forward and reverse rows (pre-expanded).
_PRE_EXPANDED_EXACT = frozenset({'Recon3D'})

_CONNECTOR_RESOURCE = 'COSMOS_formatter'
_CONNECTOR_ITYPE = 'connector'
_CATEGORIES = ('transporter', 'receptor', 'other')


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _is_pre_expanded(resource: str) -> bool:
    """True for resources that already emit both forward and reverse rows."""
    return resource.startswith('GEM') or resource in _PRE_EXPANDED_EXACT


def _row_category(itype: str, resource: str) -> str:
    """Classify an interaction row as 'transporter', 'receptor', or 'other'."""
    if (
        itype == 'transport'
        or resource.startswith('GEM_transporter')
        or (resource == 'STITCH' and itype == 'transporter')
    ):
        return 'transporter'
    if itype == 'ligand_receptor' or (resource == 'STITCH' and itype == 'receptor'):
        return 'receptor'
    return 'other'


def _fmt_met(chebi_id: str, comp: str) -> str:
    """Format a metabolite COSMOS node ID."""
    return f'Metab__{chebi_id}_{comp}' if comp else f'Metab__{chebi_id}'


def _fmt_gene(gene_id: str, n: int, reverse: bool = False) -> str:
    """Format a gene COSMOS node ID."""
    node = f'Gene{n}__{gene_id}'
    return node + '_rev' if reverse else node


def _fmt_rxn(rxn_id: str, n: int, reverse: bool = False) -> str:
    """Format an orphan reaction COSMOS node ID.

    Used for transport reactions with no gene rule (``attrs['orphan'] = True``).
    The ``Gene{N}__orphanReac`` prefix keeps these nodes in the same namespace
    as real gene product nodes while making their pseudo-enzyme origin clear.
    """
    node = f'Gene{n}__orphanReac{rxn_id}'
    return node + '_rev' if reverse else node


def _add_gene_connectors(
    bare_gene: str,
    fmt_gene: str,
    connectors: set[tuple[str, str]],
) -> None:
    """Add one connector from the bare gene ID to the formatted gene node."""
    connectors.add((bare_gene, fmt_gene))



def _assign_n(df: pd.DataFrame) -> pd.Series:
    """
    Assign sequential reaction index N to each row.

    N resets to 1 at the start of each category (transporter / receptor /
    other).  For pre-expanded resources (GEM, Recon3D), rows that share the
    same ``reaction_id`` within a category receive the same N.  All other
    rows each consume one N value.

    Implementation uses ``pd.factorize`` (which preserves first-occurrence
    order) to avoid Python-level row iteration.

    Args:
        df: DataFrame with ``_category`` column already set.

    Returns:
        Integer Series aligned to ``df.index``.
    """
    # Vectorised pre-expanded check (no apply needed)
    pre_exp_mask = (
        df['resource'].str.startswith('GEM', na=False) |
        df['resource'].isin(_PRE_EXPANDED_EXACT)
    )

    # Build a reaction key per row:
    #   - Pre-expanded with a reaction_id → use the reaction_id string so that
    #     all rows of the same reaction share a key (and thus the same N).
    #   - Everything else → use a unique per-row sentinel so each row gets
    #     its own N.  We prefix with '__idx__' to avoid collisions with
    #     real reaction IDs.
    rxn_key = pd.Series(
        '__idx__' + df.index.astype(str), index=df.index, dtype=object
    )

    if pre_exp_mask.any():
        rxn_ids = df.loc[pre_exp_mask, 'attrs'].apply(
            lambda a: a.get('reaction_id') if isinstance(a, dict) else None
        )
        valid = rxn_ids.dropna()
        if not valid.empty:
            rxn_key.loc[valid.index] = valid.astype(str)

    # factorize within each category: preserves first-occurrence order,
    # returns 0-based codes → add 1 for 1-based N.
    n_series = pd.Series(0, index=df.index, dtype=int)
    for cat in _CATEGORIES:
        cat_mask = df['_category'] == cat
        if cat_mask.any():
            codes, _ = pd.factorize(rxn_key[cat_mask])
            n_series[cat_mask] = codes + 1

    return n_series


# ---------------------------------------------------------------------------
# Per-row formatters
# ---------------------------------------------------------------------------

def _format_pre_expanded_row(
    row: dict,
    n: int,
    connectors: set[tuple[str, str]],
) -> list[dict]:
    """
    Format one row from a pre-expanded resource (GEM / Recon3D).

    The row already exists in the correct direction; this function only
    applies the COSMOS prefix/suffix to the node IDs.

    Each compartment in ``locations`` and each UniProt AC in a
    frozenset gene ID is expanded into a separate output row.

    For orphan reactions (``attrs['orphan'] == True``, ``id_type ==
    'reaction_id'``), the protein-side node receives a
    ``Gene{N}__orphanReac`` prefix, and a connector is emitted
    from the bare reaction ID to the formatted node.  No UniProt ID
    translation is attempted for orphan nodes.
    """
    attrs = dict(row['attrs']) if isinstance(row['attrs'], dict) else {}
    is_rev = attrs.get('reverse', False)
    is_orphan = attrs.get('orphan', False)
    locs = row['locations'] if isinstance(row['locations'], tuple) else ()
    comps = locs if locs else ('',)

    is_met_source = row['source_type'] == 'small_molecule'
    bare_met = row['source'] if is_met_source else row['target']
    bare_gene_raw = row['target'] if is_met_source else row['source']
    gene_ids = sorted(bare_gene_raw) if isinstance(bare_gene_raw, frozenset) else [bare_gene_raw]

    output = []
    for comp in comps:
        fmt_met = _fmt_met(bare_met, comp)
        for gid in gene_ids:
            if is_orphan:
                fmt_gene = _fmt_rxn(gid, n, is_rev)
                connectors.add((gid, fmt_gene))
            else:
                fmt_gene = _fmt_gene(gid, n, is_rev)
                _add_gene_connectors(gid, fmt_gene, connectors)

            out = dict(row)
            out['locations'] = (comp,) if comp else ()
            out['attrs'] = {**attrs, 'cosmos_formatted': True}
            if is_met_source:
                out['source'] = fmt_met
                out['target'] = fmt_gene
            else:
                out['source'] = fmt_gene
                out['target'] = fmt_met
            output.append(out)
    return output


def _format_transporter_row(
    row: dict,
    n: int,
    connectors: set[tuple[str, str]],
) -> list[dict]:
    """
    Expand one non-pre-expanded transporter row into directed edges.

    Input row is always ``small_molecule → protein`` (met → gene) for
    TCDB, SLC, and STITCH transporter rows.

    Each compartment in ``locations`` and each UniProt AC in a frozenset
    gene ID is expanded separately.  For each ``(src_comp, gene_id)``
    combination, four directed edges are emitted following the COSMOS
    convention::

        met[src_comp] → Gene{N}__AC          (forward: met enters)
        Gene{N}__AC   → met[c]               (forward: met exits to cytoplasm)
        met[c]        → Gene{N}__AC_rev      (reverse: met enters from cytoplasm)
        Gene{N}__AC_rev → met[src_comp]      (reverse: met exits)

    """
    attrs = dict(row['attrs']) if isinstance(row['attrs'], dict) else {}
    locs = row['locations'] if isinstance(row['locations'], tuple) else ()
    comps = locs if locs else ('',)

    bare_met = row['source']
    bare_gene_raw = row['target']
    id_type_met = row['id_type_a']
    id_type_gene = row['id_type_b']
    gene_ids = sorted(bare_gene_raw) if isinstance(bare_gene_raw, frozenset) else [bare_gene_raw]

    dest_comp = 'c'  # cytoplasm is always the intracellular destination

    def _met_gene(src_met: str, tgt_gene: str, met_comp: str, a: dict) -> dict:
        r = dict(row)
        r['source'] = src_met
        r['target'] = tgt_gene
        r['source_type'] = 'small_molecule'
        r['target_type'] = 'protein'
        r['id_type_a'] = id_type_met
        r['id_type_b'] = id_type_gene
        r['locations'] = (met_comp,) if met_comp else ()
        r['attrs'] = a
        return r

    def _gene_met(src_gene: str, tgt_met: str, met_comp: str, a: dict) -> dict:
        r = dict(row)
        r['source'] = src_gene
        r['target'] = tgt_met
        r['source_type'] = 'protein'
        r['target_type'] = 'small_molecule'
        r['id_type_a'] = id_type_gene
        r['id_type_b'] = id_type_met
        r['locations'] = (met_comp,) if met_comp else ()
        r['attrs'] = a
        return r

    output = []
    for src_comp in comps:
        met_source = _fmt_met(bare_met, src_comp)
        met_dest = _fmt_met(bare_met, dest_comp)
        for gid in gene_ids:
            gene_fwd = _fmt_gene(gid, n, reverse=False)
            gene_rev = _fmt_gene(gid, n, reverse=True)
            _add_gene_connectors(gid, gene_fwd, connectors)
            _add_gene_connectors(gid, gene_rev, connectors)

            fwd_attrs = {**attrs, 'cosmos_formatted': True, 'reverse': False}
            rev_attrs = {**attrs, 'cosmos_formatted': True, 'reverse': True}

            output.extend([
                _met_gene(met_source, gene_fwd, src_comp, fwd_attrs),
                _gene_met(gene_fwd, met_dest, dest_comp, fwd_attrs),
                _met_gene(met_dest, gene_rev, dest_comp, rev_attrs),
                _gene_met(gene_rev, met_source, src_comp, rev_attrs),
            ])
    return output


def _format_receptor_row(
    row: dict,
    n: int,
    connectors: set[tuple[str, str]],
) -> list[dict]:
    """
    Format one receptor or allosteric row.

    The protein node uses the bare UniProt AC — no ``Gene{N}__`` prefix.
    Each compartment in ``locations`` and each UniProt AC in a frozenset
    gene ID is expanded into a separate output row.

    Example: ``locations=('e', 'c'), target='P00533'`` →
    two rows::

        Metab__CHEBI:xxxx_e  →  P00533
        Metab__CHEBI:xxxx_c  →  P00533
    """
    attrs = dict(row['attrs']) if isinstance(row['attrs'], dict) else {}
    locs = row['locations'] if isinstance(row['locations'], tuple) else ()
    comps = locs if locs else ('',)

    is_met_source = row['source_type'] == 'small_molecule'
    bare_met = row['source'] if is_met_source else row['target']
    bare_gene_raw = row['target'] if is_met_source else row['source']
    gene_ids = sorted(bare_gene_raw) if isinstance(bare_gene_raw, frozenset) else [bare_gene_raw]

    output = []
    for comp in comps:
        fmt_met = _fmt_met(bare_met, comp)
        for gid in gene_ids:
            out = dict(row)
            out['locations'] = (comp,) if comp else ()
            out['attrs'] = {**attrs, 'cosmos_formatted': True}
            if is_met_source:
                out['source'] = fmt_met
                out['target'] = gid
            else:
                out['source'] = gid
                out['target'] = fmt_met
            output.append(out)
    return output




def _make_connector_rows(
    connectors: set[tuple[str, str]],
    columns: list[str],
) -> pd.DataFrame:
    """Build a DataFrame of connector edges from (bare_id, formatted_id) pairs."""
    if not connectors:
        return pd.DataFrame(columns=columns)

    rows = []
    for src, tgt in sorted(connectors):
        r = {col: None for col in columns}
        r['source'] = src
        r['target'] = tgt
        r['mor'] = 1
        r['interaction_type'] = _CONNECTOR_ITYPE
        r['resource'] = _CONNECTOR_RESOURCE
        r['locations'] = ()
        r['attrs'] = {'cosmos_formatted': True}
        rows.append(r)

    return pd.DataFrame(rows, columns=columns)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def format_pkn(
    source,
    include_orphans: bool = True,
) -> 'CosmosBundle':
    """
    Apply COSMOS node-ID formatting to a translated PKN.

    Converts ChEBI / UniProt IDs to COSMOS R package-compatible node IDs
    and appends connector edges linking bare IDs to their formatted
    counterparts.

    The input must have been produced by
    :func:`~omnipath_metabo.datasets.cosmos._translate.translate_pkn`
    (i.e. all metabolite IDs are ChEBI, all protein IDs are UniProt
    accessions).

    Args:
        source:
            Either a :class:`~._bundle.CosmosBundle` (as returned by
            :func:`~omnipath_metabo.datasets.cosmos.build`) or a
            translated PKN DataFrame.  When a bundle is provided, the
            ``metabolites``, ``proteins``, and ``reactions`` components
            are carried through unchanged into the returned bundle.
        include_orphans:
            If ``True`` (default), keep rows where the enzyme node is an
            orphan pseudo-enzyme (``attrs['orphan'] == True``, i.e. a
            ``reaction_id`` string used in place of a UniProt AC).
            If ``False``, drop those rows before formatting.

    Returns:
        :class:`~._bundle.CosmosBundle` with COSMOS-formatted
        ``source`` / ``target`` node IDs in the ``network`` list (as
        :class:`~._record.CosmosEdge` namedtuples) and connector edges
        appended.  ``attrs['cosmos_formatted']`` is ``True`` on every
        edge.  Provenance lists (``metabolites``, ``proteins``,
        ``reactions``) are carried through from the input bundle, or
        are empty when a plain DataFrame was provided.

    Warning:
        Do **not** call this function more than once on the same input.
        The presence of ``attrs['cosmos_formatted'] == True`` can be
        used to detect accidental double-application.
    """
    from ._bundle import CosmosBundle
    from ._record import CosmosEdge, Interaction

    if isinstance(source, CosmosBundle):
        in_bundle = source
        if not in_bundle.network:
            return CosmosBundle(
                metabolites=in_bundle.metabolites,
                proteins=in_bundle.proteins,
                reactions=in_bundle.reactions,
            )
        df = pd.DataFrame(in_bundle.network)
    else:
        in_bundle = None
        df = source

    if df.empty:
        return CosmosBundle()

    if not include_orphans:
        mask_orphan = df['attrs'].apply(
            lambda a: isinstance(a, dict) and a.get('orphan', False)
        )
        df = df[~mask_orphan]

    if df.empty:
        return CosmosBundle()

    df = df.copy()
    df['_category'] = df.apply(
        lambda r: _row_category(r['interaction_type'], r['resource']),
        axis=1,
    )
    df['_n'] = _assign_n(df)

    output_rows: list[dict] = []
    connectors: set[tuple[str, str]] = set()
    # Columns to keep in output (drop internal helpers)
    cols = [c for c in df.columns if not c.startswith('_')]

    n_pre_expanded_in = 0
    n_transporter_in = 0
    n_receptor_in = 0
    n_simple_in = 0

    for _, row in df.iterrows():
        row_dict = row.to_dict()
        cat = row['_category']
        n = int(row['_n'])

        if _is_pre_expanded(row['resource']):
            output_rows.extend(
                _format_pre_expanded_row(row_dict, n, connectors)
            )
            n_pre_expanded_in += 1
        elif cat == 'transporter':
            output_rows.extend(
                _format_transporter_row(row_dict, n, connectors)
            )
            n_transporter_in += 1
        elif cat == 'receptor':
            output_rows.extend(
                _format_receptor_row(row_dict, n, connectors)
            )
            n_receptor_in += 1
        else:
            output_rows.extend(
                _format_receptor_row(row_dict, n, connectors)
            )
            n_simple_in += 1

    _log.info(
        '[COSMOS format] %d pre-expanded input rows → %d output rows; '
        '%d transporter input rows expanded; '
        '%d receptor + %d other input rows expanded; '
        '%d total main rows, %d connector edges added.',
        n_pre_expanded_in,
        n_pre_expanded_in,  # pre-expanded rows expand per comp/gene
        n_transporter_in,
        n_receptor_in,
        n_simple_in,
        len(output_rows),
        len(connectors),
    )

    if not output_rows:
        return CosmosBundle(
            metabolites=in_bundle.metabolites if in_bundle else [],
            proteins=in_bundle.proteins if in_bundle else [],
            reactions=in_bundle.reactions if in_bundle else [],
        )

    result = pd.DataFrame(output_rows)[cols]
    conn_df = _make_connector_rows(connectors, cols)
    final_df = pd.concat([result, conn_df], ignore_index=True)

    _edge_fields = CosmosEdge._fields
    edges = [
        CosmosEdge(**{f: row[f] for f in _edge_fields})
        for _, row in final_df.iterrows()
    ]

    return CosmosBundle(
        network=edges,
        metabolites=in_bundle.metabolites if in_bundle else [],
        proteins=in_bundle.proteins if in_bundle else [],
        reactions=in_bundle.reactions if in_bundle else [],
    )


# ---------------------------------------------------------------------------
# Category-specific format wrappers
# ---------------------------------------------------------------------------

def _filter_bundle_network(bundle, predicate) -> 'CosmosBundle':
    """
    Return a copy of *bundle* with ``network`` filtered by *predicate*.

    Provenance lists (``metabolites``, ``proteins``, ``reactions``) are
    carried through unchanged — they reflect the build-step scope, which
    is already category-specific when the bundle comes from a
    ``build_*()`` call.

    Args:
        bundle: :class:`~._bundle.CosmosBundle` to filter.
        predicate: Callable ``(Interaction) → bool``.

    Returns:
        New :class:`~._bundle.CosmosBundle` with filtered ``network``.
    """
    from ._bundle import CosmosBundle

    return CosmosBundle(
        network=[row for row in bundle.network if predicate(row)],
        metabolites=bundle.metabolites,
        proteins=bundle.proteins,
        reactions=bundle.reactions,
    )


def format_transporters(source) -> 'CosmosBundle':
    """
    Format the transporter category of a COSMOS PKN bundle.

    Convenience wrapper around :func:`format_pkn` that pre-filters
    *source* to transporter rows before formatting.  When *source* already
    comes from :func:`~._build.build_transporters`, the filter is a no-op.

    Orphan transport reactions included in *source* (from a build step
    with ``include_orphans=True``) are formatted as ``Rxn{N}__<reaction_id>``
    nodes.  To exclude them, pass ``include_orphans=False`` to the
    upstream build call before formatting.

    Args:
        source:
            :class:`~._bundle.CosmosBundle` or translated PKN DataFrame.

    Returns:
        :class:`~._bundle.CosmosBundle` with COSMOS-formatted transporter edges.
    """
    from ._bundle import CosmosBundle

    if isinstance(source, CosmosBundle):
        source = _filter_bundle_network(
            source,
            lambda row: _row_category(row.interaction_type, row.resource) == 'transporter',
        )
    return format_pkn(source)


def format_receptors(source) -> 'CosmosBundle':
    """
    Format the receptor category of a COSMOS PKN bundle.

    Convenience wrapper around :func:`format_pkn` that pre-filters
    *source* to receptor rows before formatting.  When *source* already
    comes from :func:`~._build.build_receptors`, the filter is a no-op.

    Args:
        source:
            :class:`~._bundle.CosmosBundle` or translated PKN DataFrame.

    Returns:
        :class:`~._bundle.CosmosBundle` with COSMOS-formatted receptor edges.
    """
    from ._bundle import CosmosBundle

    if isinstance(source, CosmosBundle):
        source = _filter_bundle_network(
            source,
            lambda row: _row_category(row.interaction_type, row.resource) == 'receptor',
        )
    return format_pkn(source)


def format_allosteric(source) -> 'CosmosBundle':
    """
    Format the allosteric-regulation category of a COSMOS PKN bundle.

    Convenience wrapper around :func:`format_pkn` that pre-filters
    *source* to allosteric rows before formatting.  When *source* already
    comes from :func:`~._build.build_allosteric`, the filter is a no-op.

    Filter predicate matches :func:`~._build.build_allosteric`:

    - ``interaction_type == 'allosteric_regulation'`` (BRENDA)
    - ``resource == 'STITCH'`` and ``interaction_type == 'other'`` (STITCH other)

    Args:
        source:
            :class:`~._bundle.CosmosBundle` or translated PKN DataFrame.

    Returns:
        :class:`~._bundle.CosmosBundle` with COSMOS-formatted allosteric edges.
    """
    from ._bundle import CosmosBundle

    if isinstance(source, CosmosBundle):
        source = _filter_bundle_network(
            source,
            lambda row: (
                row.interaction_type == 'allosteric_regulation' or
                (row.resource == 'STITCH' and row.interaction_type == 'other')
            ),
        )
    return format_pkn(source)


def format_enzyme_metabolite(source) -> 'CosmosBundle':
    """
    Format the enzyme-metabolite (metabolic GEM) category of a COSMOS PKN bundle.

    Convenience wrapper around :func:`format_pkn` that pre-filters
    *source* to metabolic GEM rows before formatting.  When *source* already
    comes from :func:`~._build.build_enzyme_metabolite`, the filter is a no-op.

    Filter predicate matches :func:`~._build.build_enzyme_metabolite`:

    - ``resource.startswith('GEM:')`` — metabolic GEM edges only.
      ``'GEM_transporter:...'`` resources are excluded because
      ``'GEM_transporter:...'.startswith('GEM:')`` is ``False``.

    Args:
        source:
            :class:`~._bundle.CosmosBundle` or translated PKN DataFrame.

    Returns:
        :class:`~._bundle.CosmosBundle` with COSMOS-formatted enzyme-metabolite edges.
    """
    from ._bundle import CosmosBundle

    if isinstance(source, CosmosBundle):
        source = _filter_bundle_network(
            source,
            lambda row: row.resource.startswith('GEM:'),
        )
    return format_pkn(source)
