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
BRENDA allosteric regulation processing for COSMOS PKN.

BRENDA provides enzyme allosteric regulation data, including
activators and inhibitors.
"""

from __future__ import annotations

__all__ = [
    'brenda_regulations',
]

from typing import TYPE_CHECKING

import pandas as pd
from pypath.inputs.brenda._main import allosteric_regulation

from ..network import add_gene_prefix, add_metab_prefix

if TYPE_CHECKING:
    from typing import Sequence


def brenda_regulations(
    organisms: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Generate BRENDA allosteric regulation data for specified organisms.

    Args:
        organisms: List of organism names (e.g., ['human', 'mouse']).
            Defaults to ['human'].

    Returns:
        DataFrame with columns: compound, protein_id, action, id_type, pubmeds.
    """

    if organisms is None:
        organisms = ['human']

    records = []

    for record in allosteric_regulation(organisms=organisms, limit=None):

        if not record.protein:
            continue

        for protein_id in record.protein:
            records.append({
                'compound': add_metab_prefix(record.compound),
                'protein_id': add_gene_prefix(protein_id),
                'action': record.action,
                'id_type': record.id_type,
                'pubmeds': record.pubmeds,
            })

    return pd.DataFrame(
        records,
        columns=['compound', 'protein_id', 'action', 'id_type', 'pubmeds'],
    )
