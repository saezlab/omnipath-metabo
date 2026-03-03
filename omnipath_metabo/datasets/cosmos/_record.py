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
Interaction record and typed output records for the COSMOS PKN.

``Interaction`` is the internal record yielded by all resource processors.

``CosmosEdge``, ``CosmosMetabolite``, ``CosmosProtein``, and
``CosmosReaction`` are the typed output records that make up the final
:class:`~omnipath_metabo.datasets.cosmos._bundle.CosmosBundle` returned
by the build and format pipeline.  All are plain :class:`~typing.NamedTuple`
subclasses so they compose naturally with pandas (``pd.DataFrame(records)``)
and downstream tools that iterate over typed records.
"""

from __future__ import annotations

__all__ = [
    'Interaction',
    'CosmosEdge',
    'CosmosMetabolite',
    'CosmosProtein',
    'CosmosReaction',
]

from typing import NamedTuple


class Interaction(NamedTuple):
    """
    A single interaction record from a COSMOS PKN resource.

    All resource processors yield instances of this record, providing
    a uniform format regardless of the upstream data source.
    """

    source: str
    """Source entity identifier (e.g. CID, ChEBI, UniProt)."""

    target: str
    """Target entity identifier."""

    source_type: str
    """Type of source entity: ``'small_molecule'`` or ``'protein'``."""

    target_type: str
    """Type of target entity: ``'small_molecule'`` or ``'protein'``."""

    id_type_a: str
    """Identifier type of source (e.g. ``'pubchem'``, ``'chebi'``, ``'uniprot'``)."""

    id_type_b: str
    """Identifier type of target."""

    interaction_type: str
    """Type of interaction (e.g. ``'transport'``, ``'signaling'``, ``'regulation'``)."""

    resource: str
    """Name of the upstream database."""

    mor: int
    """Mode of regulation: 1 (stimulation), -1 (inhibition), 0 (unknown)."""

    locations: tuple[str, ...] = ()
    """Subcellular compartment abbreviations (e.g. ``('e', 'r')``)."""

    attrs: dict = {}
    """Arbitrary extra attributes (e.g. ``{'reverse': True, 'reaction_id': 'MAR00001'}``)."""


# ---------------------------------------------------------------------------
# Output records — components of CosmosBundle
# ---------------------------------------------------------------------------

class CosmosEdge(NamedTuple):
    """
    One edge in the final formatted COSMOS PKN network.

    Yielded by :func:`~omnipath_metabo.datasets.cosmos._format.format_pkn`.
    Node IDs are COSMOS-formatted strings (e.g. ``'Metab__CHEBI:15422_c'``,
    ``'Gene1__P12345'``).  Convert a collection to a DataFrame with
    ``pd.DataFrame(edges)``.
    """

    source: str
    """COSMOS-formatted source node ID."""

    target: str
    """COSMOS-formatted target node ID."""

    mor: int
    """Mode of regulation: 1 (activation), -1 (inhibition), 0 (unknown)."""

    interaction_type: str
    """Interaction type (e.g. ``'transport'``, ``'ligand_receptor'``, ``'connector'``)."""

    resource: str
    """Originating resource name (e.g. ``'TCDB'``, ``'GEM:Human-GEM'``)."""

    source_type: str
    """Entity type of source: ``'small_molecule'``, ``'protein'``, or ``None`` for connectors."""

    target_type: str
    """Entity type of target."""

    locations: tuple
    """Subcellular compartments associated with the edge (may be empty)."""

    attrs: dict
    """Arbitrary metadata inherited from the resource processor."""


class CosmosMetabolite(NamedTuple):
    """
    Metabolite provenance record in the COSMOS PKN.

    One record per unique metabolite encountered during the build, providing
    the mapping from the canonical ChEBI ID back to the original source
    identifier and the database it came from.
    """

    chebi: str
    """Canonical ChEBI ID (e.g. ``'CHEBI:15422'``)."""

    original_id: str
    """Identifier as it appeared in the source database."""

    id_type: str
    """Type of *original_id* (e.g. ``'pubchem'``, ``'bigg'``, ``'metatlas'``)."""

    resource: str
    """Source database that provided this metabolite."""

    name: str = ''
    """Common name if available (e.g. from ChEBI or the source database)."""


class CosmosProtein(NamedTuple):
    """
    Protein provenance record in the COSMOS PKN.

    One record per unique protein encountered during the build, providing
    the mapping from the canonical UniProt AC back to the original source
    identifier.
    """

    uniprot: str
    """Canonical UniProt accession (e.g. ``'P00533'``)."""

    original_id: str
    """Identifier as it appeared in the source database."""

    id_type: str
    """Type of *original_id* (e.g. ``'ensp'``, ``'uniprot'``, ``'genesymbol'``)."""

    resource: str
    """Source database that provided this protein."""

    gene_symbol: str = ''
    """HGNC gene symbol if available."""


class CosmosReaction(NamedTuple):
    """
    Reaction metadata record from genome-scale metabolic models (GEMs).

    One record per unique GEM reaction that contributed edges to the PKN.
    Provides a link back to the GEM annotation for each edge in the network.
    """

    reaction_id: str
    """GEM reaction identifier (e.g. ``'MAR00001'``, ``'R_PFK'``)."""

    gem: str
    """GEM name (e.g. ``'Human-GEM'``, ``'Recon3D'``)."""

    subsystem: str = ''
    """Metabolic subsystem / pathway the reaction belongs to."""

    genes: tuple = ()
    """UniProt ACs of the catalysing gene products."""

    metabolites: tuple = ()
    """ChEBI IDs of the reaction's substrates and products."""
