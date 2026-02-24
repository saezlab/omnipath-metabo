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

__all__ = ['build', 'build_transporters', 'build_receptors', 'build_enzyme_metabolite']

import logging
from itertools import chain
from typing import TYPE_CHECKING

import pandas as pd

from ._config import config
from ._record import Interaction

_log = logging.getLogger(__name__)

try:
    from tqdm import tqdm as _tqdm
except ImportError:
    _tqdm = None


def _progress(iterable, desc: str, **kwargs):
    """Wrap *iterable* with tqdm if available, otherwise return it unchanged."""
    if _tqdm is not None:
        return _tqdm(iterable, desc=desc, **kwargs)
    return iterable
from .resources import (
    brenda_regulations,
    gem_interactions,
    mrclinksdb_interactions,
    recon3d_transporter_interactions,
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
    'recon3d': recon3d_transporter_interactions,
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
    active_names = []

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
        active_names.append(name)
        generators.append(PROCESSORS[name](**resource_args))

    for name in _progress(active_names, desc='[COSMOS] collecting resources'):
        _log.info('[COSMOS] Building %s...', name)

    df = pd.DataFrame(
        chain.from_iterable(generators),
        columns=Interaction._fields,
    )

    _log.info(
        '[COSMOS] Collected %d edges from %d resources.',
        len(df),
        len(active_names),
    )

    if translate:
        from ._translate import translate_pkn
        n_before = len(df)
        _log.info('[COSMOS] Translating IDs...')
        df = translate_pkn(df, organism=organism)
        _log.info(
            '[COSMOS] Translation complete: %d/%d edges retained.',
            len(df),
            n_before,
        )

    if cfg.get('apply_blacklist', True):
        from ._blacklist import apply_blacklist
        df = apply_blacklist(df)

    return df


def build_transporters(*args, **kwargs) -> pd.DataFrame:
    """
    Build the transporter subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only
    transporter-relevant resources (TCDB, SLC, GEM, Recon3D, STITCH)
    and post-filters to keep only transporter interactions.

    Post-filter predicate:
        - ``interaction_type == 'transport'``
        - ``resource.startswith('GEM_transporter')``
        - STITCH rows where ``interaction_type == 'transporter'``

    Args:
        *args: Passed through to :func:`build`.
        **kwargs: Passed through to :func:`build`.  ``brenda`` and
            ``mrclinksdb`` are disabled unless explicitly re-enabled.

    Returns:
        DataFrame containing only transporter interactions.
    """
    kwargs.setdefault('brenda', False)
    kwargs.setdefault('mrclinksdb', False)
    df = build(*args, **kwargs)
    mask = (
        df['interaction_type'].eq('transport') |
        df['resource'].str.startswith('GEM_transporter') |
        (df['resource'].eq('STITCH') & df['interaction_type'].eq('transporter'))
    )
    return df[mask].reset_index(drop=True)


def build_receptors(*args, **kwargs) -> pd.DataFrame:
    """
    Build the receptor subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only
    receptor-relevant resources (MRCLinksDB, STITCH) and post-filters
    to keep only receptor/ligand interactions.

    Post-filter predicate:
        - ``interaction_type == 'ligand_receptor'``
        - STITCH rows where ``interaction_type == 'receptor'``

    Args:
        *args: Passed through to :func:`build`.
        **kwargs: Passed through to :func:`build`.  ``tcdb``, ``slc``,
            ``brenda``, ``gem``, and ``recon3d`` are disabled unless
            explicitly re-enabled.

    Returns:
        DataFrame containing only receptor interactions.
    """
    kwargs.setdefault('tcdb', False)
    kwargs.setdefault('slc', False)
    kwargs.setdefault('brenda', False)
    kwargs.setdefault('gem', False)
    kwargs.setdefault('recon3d', False)
    df = build(*args, **kwargs)
    mask = (
        df['interaction_type'].eq('ligand_receptor') |
        (df['resource'].eq('STITCH') & df['interaction_type'].eq('receptor'))
    )
    return df[mask].reset_index(drop=True)


def build_enzyme_metabolite(*args, **kwargs) -> pd.DataFrame:
    """
    Build the enzyme-metabolite (metabolic) subset of the COSMOS PKN.

    Convenience wrapper around :func:`build` that enables only
    metabolic resources (BRENDA, GEM, STITCH) and post-filters to
    keep only enzyme-metabolite interactions.

    Post-filter predicate:
        - ``interaction_type == 'allosteric_regulation'``
        - ``resource.startswith('GEM:')`` (metabolic GEM edges,
          distinct from ``'GEM_transporter:'`` transport edges)
        - STITCH rows where ``interaction_type == 'other'``

    Note:
        ``'GEM_transporter:...'`` resources are excluded because
        ``'GEM_transporter:...'.startswith('GEM:')`` is ``False``.

    Args:
        *args: Passed through to :func:`build`.
        **kwargs: Passed through to :func:`build`.  ``tcdb``, ``slc``,
            ``mrclinksdb``, and ``recon3d`` are disabled unless
            explicitly re-enabled.

    Returns:
        DataFrame containing only enzyme-metabolite interactions.
    """
    kwargs.setdefault('tcdb', False)
    kwargs.setdefault('slc', False)
    kwargs.setdefault('mrclinksdb', False)
    kwargs.setdefault('recon3d', False)
    df = build(*args, **kwargs)
    mask = (
        df['interaction_type'].eq('allosteric_regulation') |
        df['resource'].str.startswith('GEM:') |
        (df['resource'].eq('STITCH') & df['interaction_type'].eq('other'))
    )
    return df[mask].reset_index(drop=True)
