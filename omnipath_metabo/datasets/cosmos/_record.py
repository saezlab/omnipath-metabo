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

"""Uniform interaction record for COSMOS PKN resources."""

from __future__ import annotations

__all__ = ['Interaction']

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
