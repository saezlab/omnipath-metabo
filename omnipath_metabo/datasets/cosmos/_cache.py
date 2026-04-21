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
Pre-build cache for COSMOS PKN components.

Builds each COSMOS PKN category for each organism as a Parquet file,
so the HTTP server can assemble them quickly without re-running the
full build pipeline.

Cache layout::

    <cache_dir>/
        transporters_9606.parquet
        receptors_9606.parquet
        allosteric_9606.parquet
        enzyme_metabolite_9606.parquet
"""

from __future__ import annotations

__all__ = ['build_cache', 'load_cached', 'list_cached', 'DEFAULT_CACHE_DIR']

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

_log = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path.home() / '.cache' / 'omnipath-metabo' / 'cosmos'

# Build functions per category — values are attribute names on ``_build``.
_CATEGORIES = {
    'transporters': 'build_transporters',
    'receptors': 'build_receptors',
    'allosteric': 'build_allosteric',
    'enzyme_metabolite': 'build_enzyme_metabolite',
    'ppi': 'build_ppi',
    'grn': 'build_grn',
}


def build_cache(
    organisms: list[int] | None = None,
    categories: list[str] | None = None,
    cache_dir: Path | str | None = None,
) -> None:
    """Build and cache PKN components as Parquet files.

    Args:
        organisms: NCBI taxonomy IDs to build for.  Defaults to ``[9606]``.
        categories: Subset of category names.  Defaults to all categories.
        cache_dir: Directory to write Parquet files to.
    """
    cache_dir = Path(cache_dir or DEFAULT_CACHE_DIR)
    cache_dir.mkdir(parents=True, exist_ok=True)

    organisms = organisms or [9606]
    categories = categories or list(_CATEGORIES)

    from omnipath_metabo.datasets.cosmos import _build

    for organism in organisms:
        for category in categories:
            build_fn_name = _CATEGORIES.get(category)
            if not build_fn_name:
                _log.warning('Unknown category: %s', category)
                continue

            build_fn = getattr(_build, build_fn_name)
            _log.info('Building %s for organism %d...', category, organism)

            try:
                bundle = build_fn(organism=organism)
                _save_bundle(bundle, category, organism, cache_dir)
                _log.info(
                    'Cached %s/%d: %d edges',
                    category, organism, len(bundle.network),
                )
            except Exception:
                _log.exception('Failed %s/%d', category, organism)


def _save_bundle(
    bundle: object,
    category: str,
    organism: int,
    cache_dir: Path,
) -> None:
    """Save the network edges of a :class:`CosmosBundle` as Parquet."""
    import pandas as pd

    records = [row._asdict() for row in bundle.network]
    df = pd.DataFrame(records)

    # Convert frozenset values to semicolon-joined strings for Parquet
    for col in ('source', 'target'):
        if col in df.columns:
            df[col] = df[col].apply(
                lambda x: ';'.join(sorted(x)) if isinstance(x, frozenset) else str(x)
            )

    path = cache_dir / f'{category}_{organism}.parquet'
    df.to_parquet(path, index=False)
    _log.info('Saved %s (%d edges)', path, len(df))


def load_cached(
    category: str,
    organism: int,
    cache_dir: Path | str | None = None,
) -> pd.DataFrame | None:
    """Load a cached PKN component, or ``None`` if not cached.

    Args:
        category: One of the category keys (e.g. ``'transporters'``).
        organism: NCBI taxonomy ID.
        cache_dir: Cache root directory.

    Returns:
        DataFrame of :class:`~._record.Interaction` rows, or ``None``.
    """
    import pandas as pd

    cache_dir = Path(cache_dir or DEFAULT_CACHE_DIR)
    path = cache_dir / f'{category}_{organism}.parquet'

    if not path.exists():
        return None

    return pd.read_parquet(path)


def list_cached(cache_dir: Path | str | None = None) -> list[dict]:
    """List available cached components.

    Args:
        cache_dir: Cache root directory.

    Returns:
        List of dicts with keys ``category``, ``organism``, ``size_mb``,
        ``modified``.
    """
    cache_dir = Path(cache_dir or DEFAULT_CACHE_DIR)

    if not cache_dir.exists():
        return []

    result = []
    for f in sorted(cache_dir.glob('*.parquet')):
        parts = f.stem.rsplit('_', 1)
        if len(parts) == 2:
            try:
                organism = int(parts[1])
            except ValueError:
                continue
            result.append({
                'category': parts[0],
                'organism': organism,
                'size_mb': round(f.stat().st_size / 1024 / 1024, 3),
                'modified': f.stat().st_mtime,
            })
    return result
