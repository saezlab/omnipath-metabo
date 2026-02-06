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
TCDB transporter processing for COSMOS PKN.

TCDB (Transporter Classification Database) is a multi-species database
with transporter substrate information. This module extracts species-specific
transporter-metabolite interactions.
"""

from __future__ import annotations

__all__ = [
    'tcdb_interactions',
]

import pandas as pd
from pypath.inputs.tcdb import _substrates as tcdb_substrates
from pypath.utils import reflists

from ..data import get_data_path
from ..location import load_location_mapping, parse_uniprot_locations
from ..network import add_gene_prefix, add_metab_prefix


def tcdb_interactions(ncbi_tax_id: int = 9606) -> pd.DataFrame:
    """
    Fetch and process TCDB substrate data for a specific species.

    TCDB is a multi-species database where entries from different species
    are mixed together. This function extracts species-specific interactions
    based on protein UniProt IDs.

    Args:
        ncbi_tax_id: NCBI taxonomy ID (default: 9606 for human).

    Returns:
        DataFrame with columns: Source (metabolite with location),
        Target (protein), database.
    """

    # Load location abbreviation mapping
    abb_data = load_location_mapping(
        get_data_path('location_abb_tcdb.csv'),
        sep=';',
        columns=['location', 'abbreviation'],
    )

    # Get species-specific protein set
    species_proteins = set(
        reflists.get_reflist('uniprot', ncbi_tax_id=ncbi_tax_id)
    )

    # Collect species-specific protein substrate data
    data = pd.DataFrame([
        {
            'Source': r.transporter_uniprot,
            'Target': r.substrate_id,
            'substrate_name': r.substrate_name,
            'location': r.location,
        }
        for r in tcdb_substrates.tcdb_substrate()
        if r.transporter_uniprot in species_proteins
    ])

    if data.empty:
        return pd.DataFrame(columns=['Source', 'Target', 'database'])

    # Parse and map locations to abbreviations
    location_results = [
        parse_uniprot_locations(loc).assign(original_row=str(idx))
        for idx, loc in enumerate(data['location'])
    ]
    location_results = [df for df in location_results if not df.empty]

    if location_results:

        final_result = (
            pd.concat(location_results, ignore_index=True)
            .merge(abb_data, on='location', how='left')
            .dropna(subset=['abbreviation'])
        )
        location_summary = (
            final_result
            .groupby('original_row')['abbreviation']
            .apply(lambda x: list(x.unique()))
            .reset_index(name='locations')
        )

    else:
        location_summary = pd.DataFrame(columns=['original_row', 'locations'])

    # Create final TCDB dataframe
    result = (
        data
        .assign(row_id=data.index.astype(str))
        .merge(
            location_summary,
            left_on='row_id',
            right_on='original_row',
            how='left',
        )
        .dropna(subset=['locations'])
        .explode('locations')
    )

    # Add prefixes and format columns
    result['Source'] = result.apply(
        lambda row: add_metab_prefix(row['Target'], row['locations']),
        axis=1,
    )
    result['Target'] = result['Source_x'].apply(add_gene_prefix)
    result['database'] = 'TCDB'

    return result[['Source', 'Target', 'database']]
