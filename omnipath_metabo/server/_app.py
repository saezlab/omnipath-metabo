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
Litestar application factory for omnipath-metabo.

Serves pre-built COSMOS PKN data from Parquet cache files.
No database required; reads Parquet, filters, returns JSON.
"""

from __future__ import annotations

__all__ = ['create_app']

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from litestar import Litestar

_SERVER_DIR = Path(__file__).parent
_STATIC_DIR = _SERVER_DIR / 'static'
_LANDING_HTML = (_SERVER_DIR / 'templates' / 'landing.html').read_text(
    encoding='utf-8',
)


def create_app(cache_dir: Path | str | None = None) -> 'Litestar':
    """Create the Litestar application.

    Args:
        cache_dir: Path to the pre-built Parquet cache directory.
            Defaults to ``~/.cache/omnipath-metabo/cosmos``.

    Returns:
        Configured :class:`~litestar.Litestar` instance.
    """
    import os

    from litestar import Litestar, get
    from litestar.config.cors import CORSConfig
    from litestar.openapi import OpenAPIConfig
    from litestar.openapi.plugins import SwaggerRenderPlugin
    from litestar.static_files import create_static_files_router

    from omnipath_metabo.datasets.cosmos._cache import DEFAULT_CACHE_DIR

    from ._routes_cosmos import CosmosController
    from ._routes_networks import NetworksController

    resolved_cache_dir = Path(
        cache_dir
        or os.environ.get('OMNIPATH_METABO_CACHE_DIR')
        or DEFAULT_CACHE_DIR
    )
    omnipath_db_url = (
        os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')
    )

    @get('/', media_type='text/html', include_in_schema=False)
    async def landing() -> str:
        """Landing page."""
        return _LANDING_HTML

    @get('/health')
    async def health() -> dict[str, str]:
        """Health check endpoint."""
        return {'status': 'ok'}

    cors_config = CORSConfig(allow_origins=['*'])

    openapi_config = OpenAPIConfig(
        title='omnipath-metabo API',
        version='0.0.3',
        path='/schema',
        render_plugins=[SwaggerRenderPlugin()],
    )

    static_router = create_static_files_router(
        path='/static',
        directories=[_STATIC_DIR],
        name='static',
    )

    app = Litestar(
        route_handlers=[
            landing,
            health,
            static_router,
            CosmosController,
            NetworksController,
        ],
        cors_config=cors_config,
        openapi_config=openapi_config,
        state={
            'cache_dir': resolved_cache_dir,
            'omnipath_db_url': omnipath_db_url,
        },
    )

    return app


def _create_app_from_env() -> 'Litestar':
    """Create the app using environment-variable configuration.

    Used as the uvicorn target via module-level attribute.
    """
    return create_app()


# Module-level instance for uvicorn ``app`` reference.
# This import-time call is safe because this module is only loaded
# when the ``server`` extras are installed (litestar, uvicorn, pyarrow).
_app_instance = _create_app_from_env()

