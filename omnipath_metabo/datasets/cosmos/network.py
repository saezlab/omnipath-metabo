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
Network edge/node creation and direction handling for COSMOS PKN.

This module provides utilities for node prefixing, edge creation,
and handling reaction directionality.
"""

from __future__ import annotations

__all__ = [
    'gene_prefix',
    'metab_prefix',
    'strip_prefix',
    'reverse_reactions',
    'node_mappings',
]

import re
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from typing import Literal


GENE_PREFIX = 'Gene__'
METAB_PREFIX = 'Metab__'


def gene_prefix(identifier: str) -> str:
    """
    Add the Gene__ prefix to a protein identifier.

    Args:
        identifier: Protein/gene identifier (e.g., UniProt ID).

    Returns:
        Prefixed identifier (e.g., 'Gene__P12345').
    """

    return f'{GENE_PREFIX}{identifier}'


def metab_prefix(
    identifier: str,
    compartment: str | None = None,
) -> str:
    """
    Add the Metab__ prefix to a metabolite identifier.

    Args:
        identifier: Metabolite identifier (e.g., ChEBI ID).
        compartment: Optional compartment abbreviation to append.

    Returns:
        Prefixed identifier (e.g., 'Metab__CHEBI:12345_c').
    """

    result = f'{METAB_PREFIX}{identifier}'

    if compartment:
        result = f'{result}_{compartment}'

    return result


def strip_prefix(identifier: str) -> str:
    """
    Remove Gene__ or Metab__ prefix and any compartment suffix.

    Args:
        identifier: Prefixed identifier.

    Returns:
        Original identifier without prefix or suffix.
    """

    result = re.sub(r'^(Gene__|Metab__)', '', identifier)
    result = re.sub(r'_[a-z]+$', '', result)
    result = re.sub(r'_reverse$', '', result)

    return result


def reverse_reactions(
    table: pd.DataFrame,
    source_col: str = 'Source',
    target_col: str = 'Target',
) -> pd.DataFrame:
    """
    Add forward and reverse reaction pairs with _c suffix for metabolites.

    For transporter interactions, this creates the full reaction cycle:
    - metabolite_compartment -> protein
    - protein -> metabolite_c (cytoplasm)
    - metabolite_c -> protein_reverse
    - protein_reverse -> metabolite_compartment

    Args:
        table: DataFrame with source and target columns.
        source_col: Name of source column.
        target_col: Name of target column.

    Returns:
        Extended DataFrame with reverse reactions added.
    """

    table_list = []

    for _, row in table.iterrows():

        source = row[source_col]
        target = row[target_col]
        database = row.get('database', '')

        # Replace the last compartment with _c in the source (metabolite)
        target_modified = re.sub(r'_[a-z]+$', '_c', source)

        # Add the two rows (forward direction pair)
        table_list.append({
            source_col: source,
            target_col: target,
            'database': database,
        })
        table_list.append({
            source_col: target,
            target_col: target_modified,
            'database': database,
        })

    table_pairs = pd.DataFrame(table_list)

    # Create extended table with reverse reactions
    table_extended_list = []

    for i in range(0, len(table_pairs), 2):

        if i + 1 >= len(table_pairs):
            break

        source1 = table_pairs.iloc[i][source_col]
        target1 = table_pairs.iloc[i][target_col]
        db1 = table_pairs.iloc[i]['database']

        source2 = table_pairs.iloc[i + 1][source_col]
        target2 = table_pairs.iloc[i + 1][target_col]
        db2 = table_pairs.iloc[i + 1]['database']

        # Add four rows per batch (forward + reverse)
        table_extended_list.append({
            source_col: source1,
            target_col: target1,
            'database': db1,
        })
        table_extended_list.append({
            source_col: source2,
            target_col: target2,
            'database': db2,
        })
        table_extended_list.append({
            source_col: target2,
            target_col: f'{target1}_reverse',
            'database': db2,
        })
        table_extended_list.append({
            source_col: f'{target1}_reverse',
            target_col: source1,
            'database': db1,
        })

    return pd.DataFrame(table_extended_list)


def node_mappings(
    edges: pd.DataFrame,
    source_col: str = 'Source',
    target_col: str = 'Target',
) -> dict[str, pd.DataFrame]:
    """
    Create node mappings for genes and metabolites from edge table.

    Args:
        edges: DataFrame with source and target columns.
        source_col: Name of source column.
        target_col: Name of target column.

    Returns:
        Dictionary with 'gene_mapping' and 'metabolite_mapping' DataFrames.
    """

    all_ids = pd.concat([edges[source_col], edges[target_col]]).unique()

    # Gene mappings
    gene_ids = [gid for gid in all_ids if str(gid).startswith(GENE_PREFIX)]
    gene_mapping_list = []

    for gene_id in gene_ids:
        source = strip_prefix(gene_id)
        gene_mapping_list.append({'Source': source, 'Target': gene_id})

    gene_mapping = pd.DataFrame(gene_mapping_list)

    # Metabolite mappings
    metab_ids = [mid for mid in all_ids if str(mid).startswith(METAB_PREFIX)]
    metabolite_mapping_list = []

    for metab_id in metab_ids:
        source = strip_prefix(metab_id)
        metabolite_mapping_list.append({'Source': source, 'Target': metab_id})

    metabolite_mapping = pd.DataFrame(metabolite_mapping_list)

    return {
        'gene_mapping': gene_mapping,
        'metabolite_mapping': metabolite_mapping,
    }
