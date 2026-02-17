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
Subcellular location lookup and mapping for COSMOS PKN.

This module provides utilities for querying UniProt subcellular locations
and mapping them to standardized compartment abbreviations used in
genome-scale metabolic models (e.g., BIGG database conventions).

Compartment abbreviations:
    - c: Cytoplasm/Cytosol
    - n: Nucleus
    - e: Extracellular/Cell membrane
    - r: Endoplasmic reticulum
    - g: Golgi apparatus
    - m: Mitochondrion
    - l: Lysosome
    - x: Peroxisome
    - v: Vesicle
    - eg: ER-Golgi intermediate compartment
"""

from __future__ import annotations

__all__ = [
    'uniprot_locations',
    'locations_to_abbreviations',
    'tcdb_locations',
    'tcdb_routes',
    'slc_locations',
    'slc_routes',
]

from functools import cache
from typing import TYPE_CHECKING

import pandas as pd

from .data import data_path

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pypath.inputs.uniprot import UniprotLocation


def uniprot_locations(
    organism: int = 9606,
    reviewed: bool = True,
) -> dict[str, set[UniprotLocation]]:
    """
    Query UniProt for subcellular location annotations.

    Uses the pypath UniProt API client to fetch location data for all
    proteins of the specified organism.

    Args:
        organism: NCBI taxonomy ID (default: 9606 for human).
        reviewed: If True, only query reviewed (SwissProt) proteins.

    Returns:
        Dictionary mapping UniProt IDs to sets of UniprotLocation
        namedtuples. Each namedtuple has `location` (compartment name)
        and `features` (additional annotations like membrane topology).

    Examples:
        Get locations for human proteins:

        >>> locations = uniprot_locations(organism=9606)
        >>> 'P00533' in locations  # EGFR
        True
    """

    from pypath.inputs.uniprot import uniprot_locations as _uniprot_locations

    return _uniprot_locations(organism=organism, reviewed=reviewed)


def locations_to_abbreviations(
    locations: Iterable[UniprotLocation],
    location_mapping: dict[str, str],
) -> set[str]:
    """
    Map UniProt location namedtuples to compartment abbreviations.

    Args:
        locations: Iterable of UniprotLocation namedtuples from
            UniProt location data.
        location_mapping: Dict mapping location names to abbreviations.

    Returns:
        Set of compartment abbreviations that match the input locations.

    Examples:
        >>> mapping = tcdb_locations()
        >>> from collections import namedtuple
        >>> Loc = namedtuple('UniprotLocation', ['location', 'features'])
        >>> locs = [Loc('Cytoplasm', None), Loc('Cell membrane', ('Multi-pass',))]
        >>> abbrevs = locations_to_abbreviations(locs, mapping)
        >>> 'c' in abbrevs
        True
        >>> 'e' in abbrevs
        True
    """

    return {
        location_mapping[loc.location]
        for loc in locations
        if loc.location in location_mapping
    }


@cache
def tcdb_locations() -> dict[str, str]:
    """
    TCDB location name to abbreviation mapping.

    Returns:
        Dict mapping location names (e.g., 'Cytoplasm') to
        compartment abbreviations (e.g., 'c').

    Examples:
        >>> mapping = tcdb_locations()
        >>> mapping['Cytoplasm']
        'c'
        >>> mapping['Mitochondrion']
        'm'
    """

    df = pd.read_csv(data_path('tcdb_locations.csv'))
    return {loc: abbrev for loc, abbrev in zip(df['location'], df['abbreviation'])}


@cache
def tcdb_routes() -> list[tuple[str, str]]:
    """
    TCDB transport routes between compartments.

    Transport routes define the possible directions of metabolite
    transport between subcellular compartments.

    Returns:
        List of (from, to) tuples representing transport routes.

    Examples:
        >>> routes = tcdb_routes()
        >>> ('e', 'c') in routes  # Extracellular to cytosol
        True
    """

    df = pd.read_csv(data_path('tcdb_routes.csv'))
    return [(row['from'], row['to']) for _, row in df.iterrows()]


@cache
def slc_locations() -> dict[str, str]:
    """
    SLC location name to abbreviation mapping.

    The SLC (Solute Carrier) table uses different location naming
    conventions than TCDB, often with compound locations separated
    by semicolons.

    Returns:
        Dict mapping location names (e.g., 'Plasma membrane') to
        compartment abbreviations (e.g., 'e').

    Examples:
        >>> mapping = slc_locations()
        >>> mapping['Plasma membrane']
        'e'
        >>> mapping['Lysosome']
        'l'
    """

    df = pd.read_csv(data_path('slc_locations.csv'))
    return {loc: abbrev for loc, abbrev in zip(df['location'], df['abbreviation'])}


@cache
def slc_routes() -> list[tuple[str, str]]:
    """
    SLC transport routes between compartments.

    Returns:
        List of (from, to) tuples representing transport routes.

    Examples:
        >>> routes = slc_routes()
        >>> any(to == 'c' for _, to in routes)  # Cytosol is a destination
        True
    """

    df = pd.read_csv(data_path('slc_routes.csv'))
    return [(row['from'], row['to']) for _, row in df.iterrows()]
