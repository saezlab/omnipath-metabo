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

Resource categories:

1. Metabolite-protein interactions (current)
2. GEMs (genome-scale metabolic models)
3. PPIs (protein-protein interactions)
4. GRNs (gene regulatory networks)
5. Kinase-substrate networks
"""

from __future__ import annotations

__all__ = ['build']

from itertools import chain
from typing import TYPE_CHECKING

import pandas as pd

from ._record import Interaction
from .resources import (
    brenda_regulations,
    mrclinksdb_interactions,
    slc_interactions,
    stitch_interactions,
    tcdb_interactions,
)

if TYPE_CHECKING:
    from collections.abc import Sequence


PROCESSORS = {
    'stitch': stitch_interactions,
    'tcdb': tcdb_interactions,
    'slc': slc_interactions,
    'brenda': brenda_regulations,
    'mrclinksdb': mrclinksdb_interactions,
}

DEFAULT_ARGS: dict[str, dict] = {
    'stitch': {'ncbi_tax_id': 9606, 'score_threshold': 700},
    'tcdb': {'ncbi_tax_id': 9606},
    'slc': {},
    'brenda': {'organisms': ['human']},
    'mrclinksdb': {'organism': 'human'},
}

DEFAULT_SOURCES = tuple(PROCESSORS)


def build(
    sources: Sequence[str] | None = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Build the COSMOS prior-knowledge network from multiple sources.

    Calls each requested resource processor, collects the yielded
    :class:`Interaction` records, and returns them as a single
    DataFrame.

    Args:
        sources:
            Resource names to include. Defaults to all available:
            ``('stitch', 'tcdb', 'slc', 'brenda', 'mrclinksdb')``.
        **kwargs:
            Override default arguments for individual processors.
            Pass a dict keyed by argument name, e.g.
            ``build(stitch={'score_threshold': 400})``.

    Returns:
        DataFrame with one row per interaction, columns matching
        :class:`Interaction` fields.
    """

    if sources is None:
        sources = DEFAULT_SOURCES

    generators = []

    for name in sources:

        if name not in PROCESSORS:
            raise ValueError(
                f'Unknown source: {name!r}. '
                f'Available: {list(PROCESSORS)}'
            )

        args = {**DEFAULT_ARGS.get(name, {}), **kwargs.get(name, {})}
        generators.append(PROCESSORS[name](**args))

    return pd.DataFrame(
        chain.from_iterable(generators),
        columns=Interaction._fields,
    )
