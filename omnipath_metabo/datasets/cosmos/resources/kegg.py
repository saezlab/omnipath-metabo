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
KEGG enzyme-metabolite interactions for COSMOS PKN.

First inputs_v2-native COSMOS resource: data acquisition is fully
delegated to :func:`pypath.inputs_v2.kegg_metabolic.make_kegg_resource`,
which handles download and caching via the inputs_v2 infrastructure.
This module is responsible only for converting the flat reaction dicts
produced by the KEGG parser into COSMOS :class:`~.._record.Interaction`
records.

Each KEGG reaction record carries:

- ``'Reaction'``:  bare KEGG reaction ID (e.g. ``'R00001'``)
- ``'UniProt'``:   semicolon-joined UniProt ACs for all catalysing genes
- ``'Substrate'``: list of ``{name, kegg_id, chebi}`` dicts
- ``'Product'``:   list of ``{name, kegg_id, chebi}`` dicts

Two directed edges are emitted per (compound, enzyme) pair:

- substrate → enzyme  (``source_type='small_molecule'``)
- enzyme → product    (``source_type='protein'``)

Metabolite ID priority (per compound dict):

1. ``chebi`` present  → ``id_type='chebi'`` (pass-through in translation)
2. ``kegg_id`` present → ``id_type='kegg'`` (resolved via omnipath-utils)
3. name only          → ``id_type='synonym'`` (name lookup fallback)
"""

from __future__ import annotations

__all__ = ['kegg_interactions']

import logging
from collections.abc import Generator

from .._record import Interaction

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Metabolite ID helpers
# ---------------------------------------------------------------------------

def _compound_id(compound: dict) -> tuple[str, str]:
    """
    Return ``(id_value, id_type)`` for a compound dict.

    Priority: direct ChEBI annotation > KEGG compound ID > compound name.

    Args:
        compound: Dict with keys ``'chebi'``, ``'kegg_id'``, ``'name'``.

    Returns:
        Tuple of ``(identifier_string, id_type_string)``.
    """
    chebi = compound.get('chebi') or ''
    kegg_id = compound.get('kegg_id') or ''
    name = compound.get('name') or ''

    if chebi:
        # Normalise to 'CHEBI:xxxxx'; parser emits 'chebi:xxxxx' (lowercase).
        numeric = chebi.split(':', 1)[-1]
        return f'CHEBI:{numeric}', 'chebi'

    if kegg_id:
        # Strip 'cpd:' prefix to bare compound ID (e.g. 'C00002').
        bare = kegg_id.split(':', 1)[-1] if ':' in kegg_id else kegg_id
        return bare, 'kegg'

    return name, 'synonym'


# ---------------------------------------------------------------------------
# Record converter
# ---------------------------------------------------------------------------

def _record_to_interactions(record: dict) -> Generator[Interaction, None, None]:
    """
    Yield COSMOS Interaction records from one KEGG reaction dict.

    Args:
        record: Flat dict from the KEGG parser with keys ``'Reaction'``,
            ``'UniProt'``, ``'Substrate'``, ``'Product'``.

    Yields:
        :class:`~.._record.Interaction` records with
        ``interaction_type='catalysis'``, ``resource='KEGG'``, ``mor=1``.
    """
    rxn_id = record.get('Reaction', '')
    uniprot_raw = record.get('UniProt', '')
    substrates = record.get('Substrate') or []
    products = record.get('Product') or []

    if not uniprot_raw:
        return

    uniprots = [u.strip() for u in uniprot_raw.split(';') if u.strip()]

    if not uniprots:
        return

    attrs = {'reaction_id': rxn_id}

    for compound in substrates:

        cid, cid_type = _compound_id(compound)

        if not cid:
            continue

        for uniprot in uniprots:

            yield Interaction(
                source=cid,
                target=uniprot,
                source_type='small_molecule',
                target_type='protein',
                id_type_a=cid_type,
                id_type_b='uniprot',
                interaction_type='catalysis',
                resource='KEGG',
                mor=1,
                locations=(),
                attrs=attrs,
            )

    for compound in products:

        cid, cid_type = _compound_id(compound)

        if not cid:
            continue

        for uniprot in uniprots:

            yield Interaction(
                source=uniprot,
                target=cid,
                source_type='protein',
                target_type='small_molecule',
                id_type_a='uniprot',
                id_type_b=cid_type,
                interaction_type='catalysis',
                resource='KEGG',
                mor=1,
                locations=(),
                attrs=attrs,
            )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def kegg_interactions(
    organism: int = 9606,
    **_kwargs,
) -> Generator[Interaction, None, None]:
    """
    Yield KEGG enzyme-metabolite interactions as COSMOS Interaction records.

    Data acquisition is delegated to
    :func:`pypath.inputs_v2.kegg_metabolic.make_kegg_resource`, which
    downloads and caches the KEGG REST endpoints via the inputs_v2
    infrastructure.  The raw reaction dicts are then converted to
    :class:`~.._record.Interaction` records by :func:`_record_to_interactions`.

    Supported organisms: human (9606), mouse (10090), rat (10116), and
    others defined in :data:`pypath.inputs_v2.kegg_metabolic._KEGG_ORGANISM_CODES`.

    Two directed edges are emitted per (compound, enzyme) pair:

    - ``substrate → enzyme`` (``source_type='small_molecule'``)
    - ``enzyme → product``   (``source_type='protein'``)

    Protein IDs come out of the KEGG parser as UniProt ACs (already resolved
    via the gene→EC→reaction→UniProt chain) so no protein translation is
    needed downstream.  Metabolite IDs use ChEBI when available in the KEGG
    annotation, KEGG compound ID otherwise, compound name as last resort.

    Args:
        organism:
            NCBI taxonomy ID.  Mapped to a KEGG organism code via
            :func:`pypath.inputs_v2.kegg_metabolic.kegg_organism_code`.
            Yields nothing and logs a warning for unmapped organisms.

    Yields:
        :class:`~.._record.Interaction` records with
        ``interaction_type='catalysis'``, ``resource='KEGG'``, ``mor=1``.
    """
    try:
        from pypath.inputs_v2.kegg_metabolic import kegg_organism_code, make_kegg_resource
    except ImportError:
        _log.warning(
            "[COSMOS] KEGG: pypath.inputs_v2.kegg_metabolic is not available in "
            "the installed pypath version — skipping KEGG catalysis interactions."
        )
        return

    org_code = kegg_organism_code(organism)

    if org_code is None:
        _log.warning(
            '[COSMOS] KEGG: no organism code for taxon %d — skipping. '
            'Supported organisms: %s',
            organism,
            ', '.join(
                f'{name} ({taxon})'
                for taxon, name in {
                    9606: 'hsa', 10090: 'mmu', 10116: 'rno',
                    7955: 'dre', 6239: 'cel', 7227: 'dme', 4932: 'sce',
                }.items()
            ),
        )
        return

    _log.info('[COSMOS] KEGG: loading reactions for %s (taxon %d)...', org_code, organism)
    kegg_resource = make_kegg_resource(org_code)
    n_records = 0

    for record in kegg_resource.reactions.raw():

        n_records += 1
        yield from _record_to_interactions(record)

    _log.info('[COSMOS] KEGG: processed %d reaction records for %s.', n_records, org_code)
