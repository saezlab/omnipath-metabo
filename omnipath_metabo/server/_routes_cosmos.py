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
COSMOS PKN route handlers.

Serves pre-built PKN data from cached Parquet files.
"""

from __future__ import annotations

__all__ = ['CosmosController']

from pathlib import Path
from typing import Any

import pandas as pd
from litestar import Controller, Request, get
from litestar.params import Parameter
from litestar.response import Response

from omnipath_metabo.datasets.cosmos._cache import (
    DEFAULT_CACHE_DIR,
    list_cached,
    load_cached,
)

_ALL_CATEGORIES = [
    'transporters',
    'receptors',
    'allosteric',
    'enzyme_metabolite',
]


class CosmosController(Controller):
    """Controller for COSMOS PKN endpoints."""

    path = '/cosmos'

    def _cache_dir(self, request: Request) -> Path:
        """Resolve cache directory from app state."""
        return request.app.state.get('cache_dir', DEFAULT_CACHE_DIR)

    @get('/')
    async def cosmos_root(self, request: Request) -> dict[str, Any]:
        """Root endpoint — returns cache status."""
        cache_dir = self._cache_dir(request)
        cached = list_cached(cache_dir)
        total_mb = sum(e['size_mb'] for e in cached)
        return {
            'cache_dir': str(cache_dir),
            'components': cached,
            'total_size_mb': round(total_mb, 3),
            'n_components': len(cached),
        }

    @get('/pkn')
    async def get_pkn(
        self,
        request: Request,
        organism: int = Parameter(default=9606, description='NCBI taxonomy ID'),
        categories: str = Parameter(
            default='all',
            description='Comma-separated category names, or "all"',
        ),
        resources: str | None = Parameter(
            default=None,
            required=False,
            description='Comma-separated resource filter',
        ),
        format: str = Parameter(
            default='json',
            description='Response format: "json" or "parquet"',
        ),
    ) -> Any:
        """Get assembled PKN from cached components.

        Loads pre-built Parquet files for the requested categories and
        organism, optionally filters by resource, and returns the result.
        """
        cache_dir = self._cache_dir(request)

        if categories == 'all':
            cats = _ALL_CATEGORIES
        else:
            cats = [c.strip() for c in categories.split(',')]

        frames = []
        for cat in cats:
            df = load_cached(cat, organism, cache_dir)
            if df is not None:
                frames.append(df)

        if not frames:
            return Response(
                content={'error': 'No cached data found', 'organism': organism},
                status_code=404,
            )

        result = pd.concat(frames, ignore_index=True)

        if resources:
            resource_list = [r.strip() for r in resources.split(',')]
            result = result[
                result['resource'].isin(resource_list)
            ].reset_index(drop=True)

        if format == 'parquet':
            import io

            buf = io.BytesIO()
            result.to_parquet(buf, index=False)
            return Response(
                content=buf.getvalue(),
                media_type='application/octet-stream',
                headers={
                    'Content-Disposition': (
                        f'attachment; filename=cosmos_pkn_{organism}.parquet'
                    ),
                },
            )

        import numpy as np

        # Recursively convert non-serializable values for JSON output
        def _to_json_safe(x):
            if isinstance(x, np.ndarray):
                return x.tolist()
            if isinstance(x, np.integer):
                return int(x)
            if isinstance(x, np.floating):
                return float(x)
            if isinstance(x, dict):
                return {k: _to_json_safe(v) for k, v in x.items()}
            if isinstance(x, (list, tuple)):
                return [_to_json_safe(i) for i in x]
            if hasattr(x, '__iter__') and not isinstance(x, str):
                return [_to_json_safe(i) for i in x]
            return x

        for col in result.columns:
            if result[col].dtype == object:
                result[col] = result[col].apply(_to_json_safe)

        import json as _json

        records = result.to_dict(orient='records')
        body = _json.dumps({
            'network': records,
            'meta': {
                'organism': organism,
                'categories': cats,
                'resources': resources.split(',') if resources else None,
                'total_edges': len(records),
            },
        })
        return Response(content=body, media_type='application/json')

    @get('/categories')
    async def categories(self) -> list[str]:
        """List available PKN categories."""
        return _ALL_CATEGORIES

    @get('/organisms')
    async def organisms(self, request: Request) -> list[int]:
        """List organisms with cached data."""
        cache_dir = self._cache_dir(request)
        cached = list_cached(cache_dir)
        return sorted({entry['organism'] for entry in cached})

    @get('/resources')
    async def resources(self, request: Request) -> dict[str, list[str]]:
        """List resources available per category in the cache."""
        cache_dir = self._cache_dir(request)
        cached = list_cached(cache_dir)
        result: dict[str, list[str]] = {}

        for entry in cached:
            cat = entry['category']
            df = load_cached(cat, entry['organism'], cache_dir)
            if df is not None and 'resource' in df.columns:
                resources_in_cat = sorted(df['resource'].unique().tolist())
                result[cat] = resources_in_cat

        return result

    @get('/status')
    async def status(self, request: Request) -> dict[str, Any]:
        """Cache status: available components and sizes."""
        cache_dir = self._cache_dir(request)
        cached = list_cached(cache_dir)
        total_mb = sum(e['size_mb'] for e in cached)
        return {
            'cache_dir': str(cache_dir),
            'components': cached,
            'total_size_mb': round(total_mb, 3),
            'n_components': len(cached),
        }
