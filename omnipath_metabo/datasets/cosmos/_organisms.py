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
Organism-specific configuration for COSMOS PKN builds.

Maps NCBI taxonomy IDs to genome-scale metabolic models (GEMs) and
tracks which resources support which organisms (directly or via
orthology translation from human).
"""

from __future__ import annotations

__all__ = [
    'default_gem',
    'available_gems',
    'organism_resources',
]

import logging

_log = logging.getLogger(__name__)


# Organism → default GEM name mapping.
# Each organism has a primary GEM from Metabolic Atlas / SysBioChalmers.
# When multiple GEMs exist, the first is the default.
_ORGANISM_GEMS: dict[int, list[str]] = {
    9606: ['Human-GEM'],
    10090: ['Mouse-GEM'],
    10116: ['Rat-GEM'],
    7955: ['Zebrafish-GEM'],
    7227: ['Fruitfly-GEM'],
    6239: ['Worm-GEM'],
    4932: ['yeast-GEM'],
    562: ['Ecoli-GEM'],
}


# Which resources support each organism directly.
# Resources not listed here are assumed to be human-only and require
# orthology translation.
_RESOURCE_ORGANISMS: dict[str, set[int] | str] = {
    'tcdb': 'all',         # UniProt-based, any organism
    'brenda': 'all',       # Queries by organism name
    'stitch': 'all',       # Queries by NCBI tax ID
    'mrclinksdb': {9606, 10090},  # Human + mouse
    'slc': {9606},         # Human only (curated)
    'gem': 'gem',          # Supported if organism has a GEM
    'recon3d': {9606},     # Human only
}


def default_gem(organism: int) -> str | None:
    """Return the default GEM name for an organism.

    Args:
        organism: NCBI taxonomy ID.

    Returns:
        GEM name string, or ``None`` if no GEM is available.
    """

    gems = _ORGANISM_GEMS.get(organism)
    return gems[0] if gems else None


def available_gems(organism: int) -> list[str]:
    """Return all available GEM names for an organism.

    Args:
        organism: NCBI taxonomy ID.

    Returns:
        List of GEM names (may be empty).
    """

    return list(_ORGANISM_GEMS.get(organism, []))


def organism_resources(organism: int) -> dict[str, bool]:
    """Check which resources support an organism directly.

    Args:
        organism: NCBI taxonomy ID.

    Returns:
        Dict mapping resource name to whether it supports the organism
        directly (True) or needs orthology translation (False).
    """

    result: dict[str, bool] = {}

    for resource, support in _RESOURCE_ORGANISMS.items():
        if support == 'all':
            result[resource] = True
        elif support == 'gem':
            result[resource] = default_gem(organism) is not None
        elif isinstance(support, set):
            result[resource] = organism in support
        else:
            result[resource] = False

    return result


def needs_orthology(organism: int) -> list[str]:
    """Return resources that need orthology translation for an organism.

    These are human-only resources whose data can be translated to the
    target organism via orthologous gene pairs.

    Args:
        organism: NCBI taxonomy ID.

    Returns:
        List of resource names needing orthology translation.
    """

    if organism == 9606:
        return []

    support = organism_resources(organism)
    return [r for r, direct in support.items() if not direct]
