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
SLC transporter processing for COSMOS PKN.

The SLC (Solute Carrier) table provides human transporter-metabolite
interactions with subcellular localization information.
"""

from __future__ import annotations

__all__ = [
    'slc_interactions',
]

import pandas as pd
from pypath.inputs import slc

from ..data import get_data_path
from ..location import load_location_mapping
from ..network import add_gene_prefix, add_metab_prefix


def slc_interactions() -> pd.DataFrame:
    """
    Fetch and process SLC interaction data.

    Note: SLC data is only available for human (NCBI taxonomy ID 9606).

    Returns:
        DataFrame with columns: Source (metabolite with location),
        Target (protein), database.
    """

    # Load location abbreviation mapping
    abb_data = load_location_mapping(
        get_data_path('location_abb_slc.txt'),
        sep='\t',
        columns=['localization', 'abbreviation'],
    )

    # Collect SLC interaction data
    data = pd.DataFrame([
        {
            'transporter': r.transporter.uniprot,
            'substrate': r.substrate.chebi,
            'localization': r.localization,
        }
        for r in slc.slc_interactions()
        if (
            r.substrate.chebi
            and r.localization
            and r.localization not in ('', 'Unknown', 'None')
        )
    ])

    if data.empty:
        return pd.DataFrame(columns=['Source', 'Target', 'database'])

    # Merge with location abbreviations and filter
    result = (
        data
        .merge(abb_data, on='localization', how='left')
        .dropna(subset=['abbreviation'])
    )

    # Add prefixes and format columns
    result['Source'] = result.apply(
        lambda row: add_metab_prefix(row['substrate'], row['abbreviation']),
        axis=1,
    )
    result['Target'] = result['transporter'].apply(add_gene_prefix)
    result['database'] = 'SLC'

    return result[['Source', 'Target', 'database']]
