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
Rhea reaction processing for COSMOS PKN.

Rhea is a comprehensive expert-curated database of biochemical reactions.
This module parses Rhea equations and creates edges between substrates,
enzymes, and products.
"""

from __future__ import annotations

__all__ = [
    'rhea_reactions',
]

import re
from io import StringIO
from itertools import zip_longest
from typing import TYPE_CHECKING

import pandas as pd
import pypath.share.curl as curl
import pypath.utils.reflists as reflists

if TYPE_CHECKING:
    from typing import Literal


def _smart_split(equation_side: str) -> list[str]:
    """Split equation side by + and =, handling whitespace."""

    tokens = re.split(r'\s+\+\s+|\s+=\s+', equation_side)
    return [t.strip() for t in tokens if t.strip()]


def _parse_equation_to_edges(
    row: pd.Series,
    counter: int,
) -> tuple[list[tuple[str, str, int]], int]:
    """
    Parse chemical equation into edges based on reaction direction.

    Direction logic:
    - BI/UN: Bidirectional (substrates <-> protein <-> products)
    - LR: Left to Right (substrates -> protein -> products)
    - RL: Right to Left (products -> protein -> substrates)

    Args:
        row: DataFrame row with Equation, ID, DIRECTION, chebi_name,
            ChEBI identifier columns.
        counter: Current counter for unique protein node IDs.

    Returns:
        Tuple of (edges list, updated counter).
    """

    equation = row['Equation']
    protein_id = row['ID']
    direction = row['DIRECTION']

    # Split equation into substrates (left) and products (right)
    parts = equation.split('=')

    if len(parts) != 2:
        return [], counter

    substrates = _smart_split(parts[0])
    products = _smart_split(parts[1])

    edges = []

    # Restrict to matching pairs
    min_len = min(len(substrates), len(products))
    substrates = substrates[:min_len]
    products = products[:min_len]

    # Create ChEBI ID mapping from names
    chebi_mapping = {}
    names = str(row['chebi_name']).split(';')
    ids = str(row['ChEBI identifier']).split(';')

    for name, chebi_id in zip(names, ids):

        name = name.strip()
        chebi_id = chebi_id.strip()

        if name and chebi_id:
            chebi_mapping[name] = chebi_id

    if direction in ('BI', 'UN'):
        # Bidirectional: substrates <-> protein <-> products
        forward_counters = []

        for sub, prod in zip_longest(substrates, products, fillvalue=None):

            protein_node = f'{protein_id}_{counter}'
            forward_counters.append(counter)

            if sub is not None:
                edges.append((sub, protein_node, 1))

            if prod is not None:
                edges.append((protein_node, prod, 1))

            counter += 1

        for i, (sub, prod) in enumerate(
            zip_longest(substrates, products, fillvalue=None)
        ):

            protein_node = f'{protein_id}_{forward_counters[i]}_rev'

            if prod is not None:
                edges.append((prod, protein_node, 1))

            if sub is not None:
                edges.append((protein_node, sub, 1))

    elif direction == 'LR':
        # Left to Right: only substrates -> protein -> products
        for sub, prod in zip_longest(substrates, products, fillvalue=None):

            protein_node = f'{protein_id}_{counter}'

            if sub is not None:
                edges.append((sub, protein_node, 1))

            if prod is not None:
                edges.append((protein_node, prod, 1))

            counter += 1

    elif direction == 'RL':
        # Right to Left: only products -> protein -> substrates
        for sub, prod in zip_longest(substrates, products, fillvalue=None):

            protein_node = f'{protein_id}_{counter}'

            if prod is not None:
                edges.append((prod, protein_node, 1))

            if sub is not None:
                edges.append((protein_node, sub, 1))

            counter += 1

    # Apply ChEBI ID mapping
    mapped_edges = []

    for src, dst, w in edges:

        src_mapped = chebi_mapping.get(src, src)
        dst_mapped = chebi_mapping.get(dst, dst)
        mapped_edges.append((src_mapped, dst_mapped, w))

    return mapped_edges, counter


def rhea_reactions(
    ncbi_tax_id: int = 9606,
    rhea_nodes: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Fetch and process Rhea reaction data for a specific species.

    Args:
        ncbi_tax_id: NCBI taxonomy ID (default: 9606 for human).
        rhea_nodes: Pre-loaded Rhea nodes DataFrame. If None, will be
            fetched (requires lipinet package).

    Returns:
        DataFrame with columns: source, target, weight.
    """

    if rhea_nodes is None:
        # Try to import lipinet for Rhea parsing
        try:
            from lipinet.parse_rhea import parse_rhea_data

            rhea_results = parse_rhea_data(verbose=False, use_cache=True)
            df_rhea_nodes = rhea_results['df_nodes']

        except ImportError:
            raise ImportError(
                'lipinet package required for Rhea parsing. '
                'Provide rhea_nodes DataFrame directly or install lipinet.'
            )

    else:
        df_rhea_nodes = rhea_nodes

    # Fetch Rhea to UniProt mapping
    c = curl.Curl(
        'https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv',
        large=False,
        silent=True,
    )
    df_rhea2uniprot = pd.read_csv(StringIO(c.result), sep='\t')

    # Filter to RHEA IDs only
    df_rhea_only = df_rhea_nodes[
        df_rhea_nodes['node_id'].str.startswith('RHEA:', na=False)
    ].copy()

    df_rhea_only['RHEA_ID'] = (
        df_rhea_only['node_id']
        .str.replace('RHEA:', '', regex=False)
        .astype(int)
    )

    # Merge with UniProt mapping
    df_merged = pd.merge(
        df_rhea2uniprot,
        df_rhea_only,
        on='RHEA_ID',
        how='inner',
    )

    # Filter to species-specific proteins
    unique_ids = set(df_merged['ID'].tolist())
    valid_ids = [
        protein_id
        for protein_id in unique_ids
        if reflists.check(protein_id, 'uniprot', ncbi_tax_id=ncbi_tax_id)
    ]
    df_merged = df_merged[df_merged['ID'].isin(valid_ids)]

    # Parse equations to edges
    all_edges = []
    counter = 1

    for _, row in df_merged.iterrows():

        edges, counter = _parse_equation_to_edges(row, counter)
        all_edges.extend(edges)

    return pd.DataFrame(all_edges, columns=['source', 'target', 'weight'])
