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

N is a sequential integer assigned per unique reaction within each
category (transporter / receptor / other).  The counter resets to 1
when moving to the next category.  For pre-expanded resources (GEM,
Recon3D) all rows sharing the same ``reaction_id`` receive the same N.
For every other resource each input row counts as one reaction.

Four-row transporter pattern
----------------------------
Resources that do not pre-generate reverse edges (TCDB, SLC, STITCH)
contribute one input row per metabolite-transporter pair.  The
formatter expands each such row into four directed edges following the
COSMOS transporter convention::

    Metab__X_other  →  Gene{N}__ENSG         (met enters transporter)
    Gene{N}__AC   →  Metab__X_c            (met exits on cytoplasm side)
    Metab__X_c      →  Gene{N}__AC_rev     (reverse: met enters)
    Gene{N}__AC_rev → Metab__X_other       (reverse: met exits)

The non-cytoplasm compartment (``other``) is the first non-``'c'``
entry in the ``locations`` tuple.  When ``locations`` is empty (STITCH)
both sides use the bare ``Metab__CHEBI:xxxx`` node without a compartment
suffix.

For GEM / Recon3D resources all rows are already present in the
DataFrame; the formatter only renames the node IDs.

Connector edges
---------------
After formatting, unique connector edges are appended::

    UniProtAC  →  Gene{N}__AC...       (one per unique formatted gene)
    UniProtAC  →  Gene{N}__AC..._rev   (transporter reverse genes)
    CHEBI:xxxx →  Metab__CHEBI:xxxx_c   (one per unique formatted metabolite)

These edges allow downstream tools that key measurements by bare IDs
(e.g. transcriptomics → ENSG, metabolomics → CHEBI) to attach to the
formatted network nodes.
"""

from __future__ import annotations

__all__ = ['format_pkn']

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


def _fmt_gene(gene_id, n: int, reverse: bool = False) -> str:
    """Format a gene COSMOS node ID.

    When *gene_id* is a ``frozenset`` (ambiguous 1-to-many mapping), the
    alphabetically first accession is used for the node label so that the
    formatted ID is deterministic and unique per reaction.
    """
    if isinstance(gene_id, frozenset):
        gene_id = min(gene_id)
    node = f'Gene{n}__{gene_id}'
    return node + '_rev' if reverse else node


def _fmt_rxn(rxn_id: str, n: int, reverse: bool = False) -> str:
    """Format an orphan reaction COSMOS node ID.

    Used for transport reactions with no gene rule (``attrs['orphan'] = True``).
    The reaction ID is used as the node identifier with an ``Rxn{N}__`` prefix
    to distinguish these pseudo-nodes from real gene products.
    """
    node = f'Rxn{n}__{rxn_id}'
    return node + '_rev' if reverse else node


def _add_gene_connectors(
    bare_gene,
    fmt_gene: str,
    connectors: set[tuple[str, str]],
) -> None:
    """Add connector(s) from bare gene ID(s) to the formatted gene node.

    When *bare_gene* is a ``frozenset``, one connector is added per AC so
    that every UniProt accession in the set maps to the formatted node.
    """
    if isinstance(bare_gene, frozenset):
        for ac in bare_gene:
            connectors.add((ac, fmt_gene))
    else:
        connectors.add((bare_gene, fmt_gene))


def _other_comp(locations: tuple) -> str:
    """Return the first non-cytoplasm compartment, or '' if none."""
    for comp in locations:
        if comp != 'c':
            return comp
    return ''


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
) -> dict:
    """
    Format one row from a pre-expanded resource (GEM / Recon3D).

    The row already exists in the correct direction; this function only
    applies the COSMOS prefix/suffix to the node IDs.

    For orphan reactions (``attrs['orphan'] == True``, ``id_type ==
    'reaction_id'``), the protein-side node receives an ``Rxn{N}__``
    prefix instead of ``Gene{N}__``, and a single connector is emitted
    from the bare reaction ID to the formatted node.  No UniProt ID
    translation is attempted for orphan nodes.
    """
    attrs = dict(row['attrs']) if isinstance(row['attrs'], dict) else {}
    is_rev = attrs.get('reverse', False)
    is_orphan = attrs.get('orphan', False)
    locs = row['locations'] if isinstance(row['locations'], tuple) else ()
    comp = locs[0] if locs else ''

    is_met_source = row['source_type'] == 'small_molecule'
    bare_met = row['source'] if is_met_source else row['target']
    bare_gene = row['target'] if is_met_source else row['source']

    fmt_met = _fmt_met(bare_met, comp)

    if is_orphan:
        fmt_gene = _fmt_rxn(bare_gene, n, is_rev)
        connectors.add((bare_gene, fmt_gene))
    else:
        fmt_gene = _fmt_gene(bare_gene, n, is_rev)
        _add_gene_connectors(bare_gene, fmt_gene, connectors)

    out = dict(row)
    if is_met_source:
        out['source'] = fmt_met
        out['target'] = fmt_gene
    else:
        out['source'] = fmt_gene
        out['target'] = fmt_met

    connectors.add((bare_met, fmt_met))

    attrs['cosmos_formatted'] = True
    out['attrs'] = attrs
    return out


def _format_transporter_row(
    row: dict,
    n: int,
    connectors: set[tuple[str, str]],
) -> list[dict]:
    """
    Expand one non-pre-expanded transporter row into four directed edges.

    Input row is always ``small_molecule → protein`` (met → gene) for
    TCDB, SLC, and STITCH transporter rows.

    The four output rows follow the COSMOS convention::

        met[other] → Gene{N}__ENSG         (forward: met enters)
        Gene{N}__ENSG → met[c]             (forward: met exits)
        met[c] → Gene{N}__ENSG_rev         (reverse: met enters)
        Gene{N}__ENSG_rev → met[other]     (reverse: met exits)

    When ``locations`` is empty (STITCH), both metabolite nodes use the
    same bare ``Metab__CHEBI:xxxx`` ID without a compartment suffix.
    """
    attrs = dict(row['attrs']) if isinstance(row['attrs'], dict) else {}
    locs = row['locations'] if isinstance(row['locations'], tuple) else ()

    bare_met = row['source']   # CHEBI:xxxx
    bare_gene = row['target']  # ENSG...
    id_type_met = row['id_type_a']
    id_type_gene = row['id_type_b']

    other_comp = _other_comp(locs)
    c_comp = 'c' if 'c' in locs else ''

    met_other = _fmt_met(bare_met, other_comp)
    met_c = _fmt_met(bare_met, c_comp) if locs else met_other  # STITCH: same node
    gene_fwd = _fmt_gene(bare_gene, n, reverse=False)
    gene_rev = _fmt_gene(bare_gene, n, reverse=True)

    connectors.add((bare_met, met_other))
    if met_c != met_other:
        connectors.add((bare_met, met_c))
    _add_gene_connectors(bare_gene, gene_fwd, connectors)
    _add_gene_connectors(bare_gene, gene_rev, connectors)

    fwd_attrs = {**attrs, 'cosmos_formatted': True, 'reverse': False}
    rev_attrs = {**attrs, 'cosmos_formatted': True, 'reverse': True}

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

    return [
        _met_gene(met_other, gene_fwd, other_comp, fwd_attrs),
        _gene_met(gene_fwd, met_c, c_comp, fwd_attrs),
        _met_gene(met_c, gene_rev, c_comp, rev_attrs),
        _gene_met(gene_rev, met_other, other_comp, rev_attrs),
    ]


def _format_simple_row(
    row: dict,
    n: int,
    connectors: set[tuple[str, str]],
) -> dict:
    """
    Format one receptor or other (enzyme-metabolite) row.

    Applies prefix to the gene node (no ``_rev`` suffix — these categories
    do not have a reverse representation in the COSMOS PKN).
    """
    attrs = dict(row['attrs']) if isinstance(row['attrs'], dict) else {}
    locs = row['locations'] if isinstance(row['locations'], tuple) else ()
    comp = locs[0] if locs else ''

    is_met_source = row['source_type'] == 'small_molecule'
    bare_met = row['source'] if is_met_source else row['target']
    bare_gene = row['target'] if is_met_source else row['source']

    fmt_met = _fmt_met(bare_met, comp)
    fmt_gene = _fmt_gene(bare_gene, n, reverse=False)

    out = dict(row)
    if is_met_source:
        out['source'] = fmt_met
        out['target'] = fmt_gene
    else:
        out['source'] = fmt_gene
        out['target'] = fmt_met

    connectors.add((bare_met, fmt_met))
    _add_gene_connectors(bare_gene, fmt_gene, connectors)

    attrs['cosmos_formatted'] = True
    out['attrs'] = attrs
    return out


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

    n_pre_expanded = 0
    n_transporter = 0
    n_simple = 0

    for _, row in df.iterrows():
        row_dict = row.to_dict()
        cat = row['_category']
        n = int(row['_n'])

        if _is_pre_expanded(row['resource']):
            output_rows.append(
                _format_pre_expanded_row(row_dict, n, connectors)
            )
            n_pre_expanded += 1
        elif cat == 'transporter':
            output_rows.extend(
                _format_transporter_row(row_dict, n, connectors)
            )
            n_transporter += 1
        else:
            output_rows.append(
                _format_simple_row(row_dict, n, connectors)
            )
            n_simple += 1

    _log.info(
        '[COSMOS format] %d pre-expanded rows renamed, '
        '%d transporter rows expanded (→ %d rows), '
        '%d simple rows formatted, '
        '%d connector edges added.',
        n_pre_expanded,
        n_transporter,
        n_transporter * 4,
        n_simple,
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
