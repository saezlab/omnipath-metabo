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
COSMOS PKN builder.

This module provides the main entry point for building the COSMOS
prior-knowledge network by combining data from multiple sources.
"""

from __future__ import annotations

__all__ = [
    'build',
]

from typing import TYPE_CHECKING

import pandas as pd

from .network import node_mappings, reverse_reactions
from .resources import (
    brenda_regulations,
    mrclinksdb_interactions,
    rhea_reactions,
    slc_interactions,
    tcdb_interactions,
)

if TYPE_CHECKING:
    from typing import Sequence


DEFAULT_SOURCES = ('tcdb', 'slc', 'brenda', 'mrclinksdb')


def build(
    sources: Sequence[str] | None = None,
    ncbi_tax_id: int = 9606,
    include_reverse: bool = True,
) -> dict[str, pd.DataFrame]:
    """
    Build the COSMOS prior-knowledge network from multiple sources.

    Args:
        sources: List of source names to include. Options are:
            'tcdb', 'slc', 'brenda', 'mrclinksdb', 'rhea'.
            Defaults to all sources except 'rhea' (requires lipinet).
        ncbi_tax_id: NCBI taxonomy ID (default: 9606 for human).
        include_reverse: Whether to add reverse reaction edges.

    Returns:
        Dictionary containing:
        - 'edges': Combined edge table from all sources
        - 'gene_mapping': Gene node mappings
        - 'metabolite_mapping': Metabolite node mappings
        - 'sources': Dictionary of individual source DataFrames
    """

    if sources is None:
        sources = DEFAULT_SOURCES

    source_dfs = {}
    source_processors = {
        'tcdb': lambda: tcdb_interactions(ncbi_tax_id=ncbi_tax_id),
        'slc': slc_interactions,
        'brenda': brenda_regulations,
        'mrclinksdb': lambda: mrclinksdb_interactions(
            organism='human' if ncbi_tax_id == 9606 else str(ncbi_tax_id)
        ),
        'rhea': lambda: rhea_reactions(ncbi_tax_id=ncbi_tax_id),
    }

    # Process each requested source
    for source_name in sources:

        if source_name not in source_processors:
            raise ValueError(
                f"Unknown source: {source_name}. "
                f"Available sources: {list(source_processors.keys())}"
            )

        # Skip SLC for non-human species
        if source_name == 'slc' and ncbi_tax_id != 9606:
            continue

        source_dfs[source_name] = source_processors[source_name]()

    # Combine transporter-like sources (TCDB, SLC, MRCLinksDB)
    transporter_sources = ['tcdb', 'slc', 'mrclinksdb']
    transporter_dfs = [
        source_dfs[s]
        for s in transporter_sources
        if s in source_dfs and not source_dfs[s].empty
    ]

    if transporter_dfs:

        combined = pd.concat(transporter_dfs, ignore_index=True)
        combined = combined.drop_duplicates()
        combined['mor'] = 1

        if include_reverse:
            edges = reverse_reactions(combined)
        else:
            edges = combined

    else:
        edges = pd.DataFrame(columns=['Source', 'Target', 'database'])

    # Create node mappings
    mappings = node_mappings(edges)

    return {
        'edges': edges,
        'gene_mapping': mappings['gene_mapping'],
        'metabolite_mapping': mappings['metabolite_mapping'],
        'sources': source_dfs,
    }
