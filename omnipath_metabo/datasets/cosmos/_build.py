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

ID unification (when ``translate_ids`` is ``True``):
    - Metabolite source IDs → ChEBI (via UniChem or PubChem REST API)
    - Protein target IDs → Ensembl gene IDs, ENSG (via pypath BioMart)
"""

from __future__ import annotations

__all__ = ['build']

from itertools import chain
from typing import TYPE_CHECKING

import pandas as pd

from ._config import config
from ._record import Interaction
from .resources import (
    brenda_regulations,
    gem_interactions,
    mrclinksdb_interactions,
    slc_interactions,
    stitch_interactions,
    tcdb_interactions,
)

if TYPE_CHECKING:
    from pathlib import Path


PROCESSORS = {
    'stitch': stitch_interactions,
    'tcdb': tcdb_interactions,
    'slc': slc_interactions,
    'brenda': brenda_regulations,
    'mrclinksdb': mrclinksdb_interactions,
    'gem': gem_interactions,
}


def build(
    *args: dict | Path | str,
    **kwargs,
) -> pd.DataFrame:
    """
    Build the COSMOS prior-knowledge network from multiple sources.

    Calls each requested resource processor, collects the yielded
    :class:`Interaction` records, and returns them as a single
    DataFrame.  The ``resources`` dict in the config controls both
    which resources are active and their parameters.  The top-level
    ``organism`` is injected into each resource unless overridden.

    Args:
        *args:
            Configuration overrides as dicts or YAML file paths.
        **kwargs:
            Top-level config keys.  Resource names are accepted as
            shorthand, e.g.::

                build(stitch={'score_threshold': 400})
                build(tcdb=False, slc=False)
                build(organism=10090)

    Returns:
        DataFrame with one row per interaction, columns matching
        :class:`Interaction` fields.
    """

    cfg = config(*args, **kwargs)
    organism = cfg.get('organism', 9606)
    resources = cfg.get('resources', {})
    translate = cfg.get('translate_ids', True)

    generators = []

    for name, params in resources.items():

        if params is False:
            continue

        if name not in PROCESSORS:
            raise ValueError(
                f'Unknown resource: {name!r}. '
                f'Available: {list(PROCESSORS)}'
            )

        resource_args = params if isinstance(params, dict) else {}
        resource_args.setdefault('organism', organism)
        generators.append(PROCESSORS[name](**resource_args))

    df = pd.DataFrame(
        chain.from_iterable(generators),
        columns=Interaction._fields,
    )

    if translate:
        from ._translate import translate_pkn
        df = translate_pkn(df, organism=organism)

    return df
