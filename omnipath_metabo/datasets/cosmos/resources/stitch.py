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
STITCH chemical-protein interactions for COSMOS PKN.

STITCH integrates chemical-protein interaction data from multiple
evidence channels. This module filters by confidence score and
interaction mode, normalising orientation so the small molecule
is always the source.

Each STITCH interaction is classified into one of three protein-role
categories using Guide to Pharmacology (``guidetopharma``):

- ``'receptor'``: the protein is a receptor in Guide to Pharmacology
  (GPCR, ion channel, NHR, catalytic receptor).
- ``'transporter'``: the protein is annotated as a transporter.
- ``'other'``: allosteric regulators, enzymes, and proteins not present
  in Guide to Pharmacology.

The classification is stored in the ``interaction_type`` field.  The
original STITCH mode (``'activation'``, ``'inhibition'``, ``'binding'``,
etc.) is preserved in ``attrs['stitch_mode']``.
"""

from __future__ import annotations

__all__ = ['stitch_interactions']

from collections.abc import Generator
from functools import cache
from typing import TYPE_CHECKING

from .._record import Interaction

if TYPE_CHECKING:
    from collections.abc import Sequence


# Guide to Pharmacology target_type values that map to each COSMOS category.
# Values are lowercase with underscores as returned by guidetopharma._targets.
_G2P_RECEPTOR_TYPES = frozenset({
    'gpcr',            # G protein-coupled receptors
    'catalytic_receptor',  # receptor tyrosine kinases etc.
    'nhr',             # nuclear hormone receptors
    'vgic',            # voltage-gated ion channels
    'lgic',            # ligand-gated ion channels
    'other_ic',        # other ion channels
})
_G2P_TRANSPORTER_TYPES = frozenset({'transporter'})


@cache
def _g2p_uniprot_types(organism: int) -> dict[str, str]:
    """
    Build a UniProt → interaction_type dict from Guide to Pharmacology.

    Uses the ``guidetopharma`` subpackage (the newer, actively maintained
    implementation). ``target_type`` values from the
    ``targets_and_families`` table are collapsed to three categories:
    ``'receptor'``, ``'transporter'``, or ``'other'``.

    Downloaded once per session per organism and cached in memory.

    Args:
        organism: NCBI taxonomy ID.

    Returns:
        Dict mapping UniProt accession strings to category strings.
    """

    from pypath.inputs.guidetopharma import protein_targets

    result: dict[str, str] = {}

    for targets_list in protein_targets().values():

        for t in targets_list:

            if t.organism != organism or not t.uniprot:
                continue

            if t.target_type in _G2P_RECEPTOR_TYPES:
                result[t.uniprot] = 'receptor'
            elif t.target_type in _G2P_TRANSPORTER_TYPES:
                result[t.uniprot] = 'transporter'
            else:
                result[t.uniprot] = 'other'

    return result


def _classify_protein(ensp: str, organism: int) -> str:
    """
    Return the COSMOS interaction_type for a STITCH protein (ENSP).

    Translates the ENSP to UniProt via pypath's BioMart-backed mapping,
    then looks up the UniProt in the Guide to Pharmacology type table.
    Returns ``'other'`` if translation fails or the protein is absent
    from Guide to Pharmacology.

    Args:
        ensp: Ensembl protein identifier (ENSP...).
        organism: NCBI taxonomy ID.

    Returns:
        One of ``'receptor'``, ``'transporter'``, or ``'other'``.
    """

    import pypath.utils.mapping as mapping_mod

    g2p_types = _g2p_uniprot_types(organism)
    uniprots = mapping_mod.map_name(ensp, 'ensp', 'uniprot', ncbi_tax_id=organism)

    for uniprot in uniprots:
        itype = g2p_types.get(uniprot)

        if itype:
            return itype

    return 'other'


def stitch_interactions(
    organism: int = 9606,
    score_threshold: int = 700,
    mode: str | Sequence[str] | None = ('activation', 'inhibition', 'binding'),
    a_is_acting: bool = True,
) -> Generator[Interaction, None, None]:
    """
    Yield STITCH chemical-protein interactions as uniform records.

    Each protein target is classified by its biological role using
    Guide to Pharmacology: ``'receptor'``, ``'transporter'``, or
    ``'other'`` (allosteric / unknown binding).  The classification
    is stored in ``interaction_type``; the original STITCH mode is
    preserved in ``attrs['stitch_mode']``.

    Args:
        organism:
            NCBI taxonomy ID (default: 9606 for human).
        score_threshold:
            Minimum ``final_score`` for interactions.
        mode:
            Interaction mode(s) to keep.  Pass a single string, a
            sequence of strings, or ``None`` to keep all modes.
            Available modes: ``'binding'``, ``'pred_bind'``,
            ``'expression'``, ``'activation'``, ``'inhibition'``,
            ``'reaction'``, ``'catalysis'``.
        a_is_acting:
            If ``True`` (default), keep only interactions where the
            chemical is the directed causal agent (``a_is_acting`` flag
            in the STITCH actions file).  Undirected interactions are
            dropped.  Set to ``False`` to include all interactions
            regardless of directionality.  Following the legacy
            Generation-PKN-COSMOS behaviour.

    Yields:
        :class:`Interaction` records with *source_type*
        ``'small_molecule'`` and *target_type* ``'protein'``.
        ``interaction_type`` is one of ``'receptor'``,
        ``'transporter'``, or ``'other'``.
    """

    from pypath.inputs.new_stitch import interactions

    allowed = None

    if mode is not None:
        allowed = (mode,) if isinstance(mode, str) else tuple(mode)

    for rec in interactions(ncbi_tax_id=organism):

        if rec.final_score < score_threshold:
            continue

        if allowed and rec.mode not in allowed:
            continue

        if a_is_acting and not rec.directed:
            continue

        # Normalise orientation: small molecule as source
        if rec.source.type == 'small_molecule':
            chemical_id = rec.source.id
            protein_id = rec.target.id
        elif rec.target.type == 'small_molecule':
            chemical_id = rec.target.id
            protein_id = rec.source.id
        else:
            continue

        if rec.activation:
            mor = 1
        elif rec.inhibition:
            mor = -1
        else:
            mor = 0

        interaction_type = _classify_protein(protein_id, organism)

        yield Interaction(
            source=chemical_id,
            target=protein_id,
            source_type='small_molecule',
            target_type='protein',
            id_type_a='pubchem',
            id_type_b='ensp',
            interaction_type=interaction_type,
            resource='STITCH',
            mor=mor,
            attrs={'stitch_mode': rec.mode},
        )
