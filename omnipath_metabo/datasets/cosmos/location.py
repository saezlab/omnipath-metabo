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
Subcellular location parsing and mapping for COSMOS PKN.

This module provides utilities for parsing UniProt location strings
and mapping them to standardized compartment abbreviations.
"""

from __future__ import annotations

__all__ = [
    'parse_uniprot_locations',
    'load_location_mapping',
]

import re
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from typing import Sequence


def parse_uniprot_locations(location_string: str | None) -> pd.DataFrame:
    """
    Parse UniProt location strings and extract location and features.

    Args:
        location_string: UniProt location string in the format
            `{UniprotLocation(location='...', features=...)}`.

    Returns:
        DataFrame with 'location' and 'features' columns.
    """

    if pd.isna(location_string) or not isinstance(location_string, str):
        return pd.DataFrame(columns=['location', 'features'])

    pattern = r"UniprotLocation\(location='([^']+)', features=([^)]+)\)"
    cleaned = re.sub(r'^\{|\}$', '', str(location_string))
    matches = re.findall(pattern, cleaned)

    if not matches:
        return pd.DataFrame(columns=['location', 'features'])

    return pd.DataFrame([
        {
            'location': m[0],
            'features': re.sub(r"\('|'\)|'|,\s*$|^\(|\)$", '', m[1]),
        }
        for m in matches
    ])


def load_location_mapping(
    filepath: str | Path,
    sep: str = ';',
    columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Load and clean location abbreviation mapping from a file.

    Args:
        filepath: Path to the CSV/TSV file containing location mappings.
        sep: Column separator.
        columns: Column names to assign if file has no header.

    Returns:
        DataFrame with location mapping (typically 'location' and
        'abbreviation' columns).
    """

    df = pd.read_csv(filepath, sep=sep, header=None if columns else 0)
    df = df.iloc[:, :2]
    df.columns = columns or df.columns[:2]

    if 'location' in df.columns:
        df['location'] = df['location'].str.strip()

    return df
