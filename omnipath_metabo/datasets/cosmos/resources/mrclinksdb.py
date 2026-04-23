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
MRCLinksDB receptor-metabolite processing for COSMOS PKN.

MRCLinksDB provides receptor-metabolite interaction data with
subcellular localization information.

Each protein is classified as ``'receptor'``, ``'transporter'``, or
``'other'`` using the same multi-source strategy as STITCH
(Guide to Pharmacology > OmniPath Intercell > TCDB).  The classification
is stored in ``interaction_type``:

- ``'ligand_receptor'``: receptor or unclassified protein
- ``'transport'``: transporter protein
"""

from __future__ import annotations

import logging

_log = logging.getLogger(__name__)

__all__ = ['mrclinksdb_interactions', 'mrclinksdb_transporter_protein_interactions']

from collections.abc import Generator

from .._record import Interaction
from ..location import ORGANISM_NAMES


def mrclinksdb_interactions(
    organism: int = 9606,
) -> Generator[Interaction, None, None]:
    """
    Yield MRCLinksDB receptor-metabolite interactions as uniform records.

    Args:
        organism:
            NCBI taxonomy ID (default: 9606 for human).

    Yields:
        :class:`Interaction` records with *source_type*
        ``'small_molecule'`` and *target_type* ``'protein'``.
    """

    from pypath.inputs.mrclinksdb import _interactions

    from ..location import (
        resolve_protein_locations,
        tcdb_locations,
        uniprot_locations,
    )

    from .stitch import _multidb_uniprot_types

    organism_name = ORGANISM_NAMES.get(organism, str(organism))
    location_mapping = tcdb_locations()
    # For non-human organisms, include unreviewed (TrEMBL) entries — many
    # non-human proteins are not in SwissProt and would otherwise all be dropped.
    reviewed = organism == 9606
    all_locations = uniprot_locations(organism=organism, reviewed=reviewed)
    protein_types = _multidb_uniprot_types(organism)

    try:
        mrclinks_data = list(
            _interactions.mrclinksdb_interaction(organism=organism_name)
        )
    except (TypeError, Exception):
        # MRCLinksDB doesn't support this organism (download returns None)
        _log.info(
            '[COSMOS] MRCLinksDB: no data for organism %s, skipping.',
            organism_name,
        )
        return

    for rec in mrclinks_data:

        receptor = str(rec.receptor_uniprot)
        pubchem_raw = str(rec.pubchem)

        # PubChem IDs arrive with a 'CID:' prefix (e.g. 'CID:13712').
        # Records that do not match this pattern (including the CSV header
        # row that pypath passes through) are skipped.
        if not pubchem_raw.startswith('CID:'):
            continue

        pubchem = pubchem_raw[4:]  # strip 'CID:' → plain numeric string

        if not receptor or not pubchem:
            continue

        # MRCLinksDB catalogues only plasma-membrane interactions, so 'e'
        # is the correct fallback for proteins without UniProt location data
        # (common for unreviewed non-human entries).
        abbreviations = resolve_protein_locations(receptor, all_locations, location_mapping) or {'e'}

        ptype = protein_types.get(receptor, 'other')
        interaction_type = 'transport' if ptype == 'transporter' else 'ligand_receptor'

        yield Interaction(
            source=pubchem,
            target=receptor,
            source_type='small_molecule',
            target_type='protein',
            id_type_a='pubchem',
            id_type_b='uniprot',
            interaction_type=interaction_type,
            resource='MRCLinksDB',
            mor=1,
            locations=tuple(sorted(abbreviations)),
        )


def mrclinksdb_transporter_protein_interactions(
    organism: int = 9606,
) -> Generator[Interaction, None, None]:
    """
    Yield MRCLinksDB transporter-protein interactions from the dedicated
    transporter protein file.

    Unlike :func:`mrclinksdb_interactions`, which parses the ligand-receptor
    file and infers transport vs. receptor class from protein annotations, this
    function reads the separate organism-specific transporter protein file
    (``<Organism> transporter protein.txt``).  Every protein in that file is
    a confirmed transporter, so no protein classification is needed.

    Metabolite IDs are HMDB accessions (``id_type_a='hmdb'``); protein IDs
    are UniProt accessions already present in the file (``id_type_b='uniprot'``
    is a pass-through in the translation step).

    Args:
        organism:
            NCBI taxonomy ID (default: 9606 for human).

    Yields:
        :class:`Interaction` records with *source_type* ``'small_molecule'``
        and *target_type* ``'protein'``, ``interaction_type='transport'``,
        ``resource='MRCLinksDB_transporter'``.
    """

    from pypath.inputs.mrclinksdb import _interactions

    from ..location import (
        resolve_protein_locations,
        tcdb_locations,
        uniprot_locations,
    )

    organism_name = ORGANISM_NAMES.get(organism, str(organism))
    location_mapping = tcdb_locations()
    reviewed = organism == 9606
    all_locations = uniprot_locations(organism=organism, reviewed=reviewed)

    try:
        data = list(
            _interactions.mrclinksdb_transporter_interaction(organism=organism_name)
        )
    except Exception:
        _log.info(
            '[COSMOS] MRCLinksDB transporter: no data for organism %s, skipping.',
            organism_name,
        )
        return

    for rec in data:

        uniprot_id = str(rec.transporter_uniprot)
        hmdb_id = str(rec.hmdb)

        if not hmdb_id or not uniprot_id:
            continue

        # All proteins in this file are plasma-membrane transporters; use 'e'
        # as fallback for entries without UniProt location data.
        abbreviations = (
            resolve_protein_locations(uniprot_id, all_locations, location_mapping)
            or {'e'}
        )

        yield Interaction(
            source=hmdb_id,
            target=uniprot_id,
            source_type='small_molecule',
            target_type='protein',
            id_type_a='hmdb',
            id_type_b='uniprot',
            interaction_type='transport',
            resource='MRCLinksDB_transporter',
            mor=1,
            locations=tuple(sorted(abbreviations)),
        )
