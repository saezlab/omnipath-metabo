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
Expert-curation blacklist for the COSMOS PKN.

The blacklist is a list of filter entries loaded from
``data/blacklist.yaml``.  Each entry is a dict of column → value pairs;
an interaction is removed if **all** fields in an entry match (AND logic
within an entry).  Multiple entries are combined with OR logic — i.e. a
row is dropped if it matches *any* entry.

Matching is exact string equality on the translated DataFrame columns
(``source``, ``target``, ``resource``, ``interaction_type``,
``source_type``, ``target_type``).

Example entries in ``blacklist.yaml``::

    blacklist:
      # Remove one specific STITCH edge known to be a false positive:
      - source: CHEBI:15422
        target: ENSG00000001234
        resource: STITCH

      # Remove all Recon3D transport edges for a known promiscuous metabolite:
      - source: CHEBI:57618
        resource: Recon3D
"""

from __future__ import annotations

__all__ = ['apply_blacklist']

import logging
from typing import TYPE_CHECKING

import pandas as pd
import yaml

from .data import data_path

if TYPE_CHECKING:
    pass

_log = logging.getLogger(__name__)

_BLACKLIST_FILE = 'blacklist.yaml'


def _load_default_blacklist() -> list[dict]:
    """Load the built-in expert-curation blacklist entries."""

    path = data_path(_BLACKLIST_FILE)

    with path.open() as f:
        data = yaml.safe_load(f) or {}

    return data.get('blacklist', [])


def apply_blacklist(
    df: pd.DataFrame,
    entries: list[dict] | None = None,
) -> pd.DataFrame:
    """
    Remove blacklisted interactions from the PKN DataFrame.

    Each entry in *entries* is a dict of column → value pairs.  An
    interaction (row) is removed if **all** field conditions in an entry
    match.  Multiple entries are combined with OR logic.

    Applies best after ID translation so that IDs are in the canonical
    form (ChEBI, ENSG) that the blacklist entries reference.

    Args:
        df:
            PKN DataFrame as returned by
            :func:`~omnipath_metabo.datasets.cosmos.build` (possibly
            after translation).
        entries:
            Blacklist entries.  Each entry is a dict mapping a column
            name to an exact match value.  If ``None``, entries are
            loaded from the built-in ``data/blacklist.yaml``.

    Returns:
        DataFrame with blacklisted interactions removed.  Index is
        reset.  If no entries match, the original DataFrame is returned
        unchanged (copy).
    """

    if entries is None:
        entries = _load_default_blacklist()

    entries = [e for e in entries if isinstance(e, dict) and e]

    if not entries:
        return df

    blacklist_mask = pd.Series(False, index=df.index)

    for entry in entries:
        entry_mask = pd.Series(True, index=df.index)

        for col, val in entry.items():
            if col not in df.columns:
                _log.warning(
                    'Blacklist entry references unknown column %r — skipping field.',
                    col,
                )
                continue

            entry_mask &= df[col].astype(str) == str(val)

        blacklist_mask |= entry_mask

    n_removed = int(blacklist_mask.sum())

    if n_removed:
        _log.info('Blacklist removed %d interaction(s).', n_removed)

    return df[~blacklist_mask].reset_index(drop=True)
