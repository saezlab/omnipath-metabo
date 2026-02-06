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
"""

from __future__ import annotations

__all__ = [
    'mrclinksdb_interactions',
]

import pandas as pd
from pypath.inputs.mrclinksdb import _interactions

from ..data import get_data_path
from ..location import load_location_mapping, parse_uniprot_locations
from ..network import add_gene_prefix, add_metab_prefix


def mrclinksdb_interactions(organism: str = 'human') -> pd.DataFrame:
    """
    Process MRCLinksDB interactions for a specific organism.

    Args:
        organism: Organism name (e.g., 'human', 'mouse').

    Returns:
        DataFrame with columns: Source (metabolite with location),
        Target (protein), database.
    """

    # Load location abbreviation mapping (reuse TCDB mapping)
    abb_data = load_location_mapping(
        get_data_path('location_abb_tcdb.csv'),
        sep=';',
        columns=['location', 'abbreviation'],
    )

    # Fetch interactions from pypath
    interactions = list(
        _interactions.mrclinksdb_interaction(organism=organism)
    )

    # Convert to DataFrame
    data = pd.DataFrame([
        {
            'Target': rec.pubchem,
            'Source': str(rec.receptor_uniprot),
            'location': rec.receptor_location,
        }
        for rec in interactions
    ])

    if data.empty:
        return pd.DataFrame(columns=['Source', 'Target', 'database'])

    # Parse locations for each row
    all_parsed_locations = []

    for idx, row in data.iterrows():

        parsed = parse_uniprot_locations(row['location'])

        if parsed.empty:
            continue

        for _, loc_row in parsed.iterrows():
            all_parsed_locations.append({
                'location': loc_row['location'],
                'features': loc_row['features'],
                'original_row': str(idx),
            })

    if not all_parsed_locations:
        return pd.DataFrame(columns=['Source', 'Target', 'database'])

    final_result = pd.DataFrame(all_parsed_locations)

    # Merge with abbreviation data and filter
    final_result = (
        final_result
        .merge(abb_data, on='location', how='left')
        .dropna(subset=['abbreviation'])
    )

    # Group by original row to get unique locations
    location_summary = (
        final_result
        .groupby('original_row')['abbreviation']
        .apply(lambda x: list(x.unique()))
        .reset_index(name='locations')
    )

    # Merge back with original data
    data['row_id'] = data.index.astype(str)
    result_table = (
        data
        .merge(
            location_summary,
            left_on='row_id',
            right_on='original_row',
            how='left',
        )
        .dropna(subset=['locations'])
        .explode('locations')
        .reset_index(drop=True)
    )

    # Create Gene/Metab identifiers
    result_table['Source_new'] = result_table.apply(
        lambda row: add_metab_prefix(row['Target'], row['locations']),
        axis=1,
    )
    result_table['Target_new'] = result_table.apply(
        lambda row: add_gene_prefix(row['Source']),
        axis=1,
    )
    result_table['database'] = 'MRCLinksDB'

    return result_table[['Source_new', 'Target_new', 'database']].rename(
        columns={'Source_new': 'Source', 'Target_new': 'Target'}
    )
