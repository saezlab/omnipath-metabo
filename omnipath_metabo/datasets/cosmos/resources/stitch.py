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
categories using a two-source annotation strategy:

1. **Guide to Pharmacology** (``guidetopharma``) classifies GPCRs, ion
   channels, nuclear hormone receptors, and catalytic receptors as
   ``'receptor'``; its ``'transporter'`` type maps to ``'transporter'``.

2. **TCDB** (Transporter Classification Database) provides exhaustive
   transporter coverage beyond what Guide to Pharmacology annotates.
   Any UniProt present in TCDB that is not already classified as a
   receptor by Guide to Pharmacology is assigned ``'transporter'``.

Priority order: receptor (G2P or Intercell) > transporter (G2P, TCDB,
or Intercell) > other.

- ``'receptor'``: protein is a GPCR, ion channel, NHR, or catalytic
  receptor in Guide to Pharmacology, or annotated as a receptor in
  OmniPath Intercell (aggregates CellPhoneDB, HPMR, CellChatDB, etc.).
- ``'transporter'``: protein is a transporter in Guide to Pharmacology,
  present in TCDB, or annotated as a transporter in OmniPath Intercell.
- ``'other'``: allosteric regulators, enzymes, and proteins absent from
  all three sources.

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


# Guide to Pharmacology target_type values that map to 'receptor'.
# Values are lowercase with underscores as returned by guidetopharma.
_G2P_RECEPTOR_TYPES = frozenset({
    'gpcr',                # G protein-coupled receptors
    'catalytic_receptor',  # receptor tyrosine kinases etc.
    'nhr',                 # nuclear hormone receptors
    'vgic',                # voltage-gated ion channels
    'lgic',                # ligand-gated ion channels
    'other_ic',            # other ion channels
})
_G2P_TRANSPORTER_TYPES = frozenset({'transporter'})


@cache
def _multidb_uniprot_types(organism: int) -> dict[str, str]:
    """
    Build a UniProt → interaction_type dict from OmniPath Intercell,
    TCDB, and Guide to Pharmacology (union).

    The result is cached to disk (pickle) in the pypath cache directory
    so that the expensive intercell database build is skipped on
    subsequent sessions.  Delete the file to force a rebuild.

    Classification priority (lowest → highest):

    1. OmniPath Intercell transporter set → ``'transporter'``.
    2. TCDB annotations → ``'transporter'``.
    3. OmniPath Intercell receptor set → ``'receptor'``.
    4. Guide to Pharmacology transporter type → ``'transporter'``
       (unless already ``'receptor'``).
    5. Guide to Pharmacology receptor types → ``'receptor'``
       (highest priority; overrides everything).

    OmniPath Intercell aggregates CellPhoneDB, HPMR, CellChatDB, and
    many other sources.  TCDB provides exhaustive transporter coverage.
    Guide to Pharmacology has the strongest receptor curation and
    highest classification priority.

    Intercell db build is expensive on first call (~minutes); subsequent
    calls within the same Python session are instant.  If intercell
    fails for any reason the function falls back to TCDB + G2P.

    Downloaded once per session per organism and cached in memory.

    Args:
        organism: NCBI taxonomy ID.

    Returns:
        Dict mapping UniProt accession strings to ``'receptor'`` or
        ``'transporter'``.  Proteins absent from all sources are not
        present in the dict (caller defaults to ``'other'``).
    """

    import pickle
    from pathlib import Path

    import pypath.share.settings as _settings

    cache_path = Path(_settings.get('cachedir')) / f'cosmos_protein_types_{organism}.pkl'

    if cache_path.exists():
        with cache_path.open('rb') as _f:
            return pickle.load(_f)

    from pypath.inputs.guidetopharma import protein_targets
    from pypath.inputs.tcdb import tcdb_classes

    result: dict[str, str] = {}

    # --- Intercell: broad multi-source baseline (lowest priority) ---
    try:
        from pypath.core import intercell as intercell_mod

        db = intercell_mod.get_db()

        for entity in db.select('transporter'):
            if isinstance(entity, str):
                result[entity] = 'transporter'

        for entity in db.select('receptor'):
            if isinstance(entity, str):
                result[entity] = 'receptor'

    except Exception:
        pass  # intercell unavailable; TCDB + G2P still applied below

    # --- TCDB: exhaustive transporter coverage ---
    # Every entry in TCDB is a transporter by definition; no family-name
    # lookup needed.
    for uniprot in tcdb_classes():
        result[uniprot] = 'transporter'

    # --- G2P: highest-priority receptor curation ---
    for targets_list in protein_targets().values():

        for t in targets_list:

            if t.organism != organism or not t.uniprot:
                continue

            if t.target_type in _G2P_RECEPTOR_TYPES:
                result[t.uniprot] = 'receptor'
            elif t.target_type in _G2P_TRANSPORTER_TYPES:
                if result.get(t.uniprot) != 'receptor':
                    result[t.uniprot] = 'transporter'

    with cache_path.open('wb') as _f:
        pickle.dump(result, _f)

    return result


def _classify_protein(ensp: str, organism: int) -> str:
    """
    Return the COSMOS interaction_type for a STITCH protein (ENSP).

    Translates the ENSP to UniProt via pypath's BioMart-backed mapping,
    then looks up the UniProt in the multi-source annotation table
    (Guide to Pharmacology + TCDB).  Returns ``'other'`` if translation
    fails or the protein is absent from both sources.

    Args:
        ensp: Ensembl protein identifier (ENSP...).
        organism: NCBI taxonomy ID.

    Returns:
        One of ``'receptor'``, ``'transporter'``, or ``'other'``.
    """

    import pypath.utils.mapping as mapping_mod

    types = _multidb_uniprot_types(organism)
    uniprots = mapping_mod.map_name(ensp, 'ensp', 'uniprot', ncbi_tax_id=organism)

    for uniprot in uniprots:
        itype = types.get(uniprot)

        if itype:
            return itype

    return 'other'


def stitch_interactions(
    organism: int = 9606,
    score_threshold: int = 700,
    mode: str | Sequence[str] | None = ('activation', 'inhibition'),
    a_is_acting: bool = True,
) -> Generator[Interaction, None, None]:
    """
    Yield STITCH chemical-protein interactions as uniform records.

    Each protein target is classified by its biological role using
    Guide to Pharmacology and TCDB: ``'receptor'``, ``'transporter'``,
    or ``'other'`` (allosteric / unknown binding).  The classification
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

        ptype = _classify_protein(protein_id, organism)
        interaction_type = 'ligand_receptor' if ptype == 'receptor' else ptype

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
