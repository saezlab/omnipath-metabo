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

    from omnipath_metabo.datasets.cosmos._cache import DEFAULT_CACHE_DIR

    from ._routes_cosmos import CosmosController

    resolved_cache_dir = Path(
        cache_dir
        or os.environ.get('OMNIPATH_METABO_CACHE_DIR')
        or DEFAULT_CACHE_DIR
    )

    from litestar.response import Response as LitestarResponse

    @get('/', media_type='text/html', include_in_schema=False)
    async def landing() -> LitestarResponse:
        """Landing page."""
        return LitestarResponse(content=_LANDING_HTML, media_type='text/html')

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

    app = Litestar(
        route_handlers=[landing, health, CosmosController],
        cors_config=cors_config,
        openapi_config=openapi_config,
        state={'cache_dir': resolved_cache_dir},
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


_LANDING_HTML = """\
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>OmniPath Metabolomics</title>
    <link href="https://fonts.googleapis.com/css?family=Raleway:400,300,600" rel="stylesheet">
    <link rel="stylesheet" href="https://omnipathdb.org/css/normalize.css">
    <link rel="stylesheet" href="https://omnipathdb.org/css/barebones.css">
    <link rel="stylesheet" href="https://omnipathdb.org/css/omnipath.css">
    <link rel="icon" type="image/vnd.microsoft.icon" href="https://omnipathdb.org/favicon.ico">
    <style>
        h5 { font-weight: 600; }
        .box code, .box a.code { font-size: 85%; }
        pre code { font-size: 85%; }
    </style>
</head>
<body>
    <!-- Header: logo + title -->
    <div class="grid-container u-align-left thirds">
        <div>
            <img src="https://omnipathdb.org/img/omnipath_logo.png"
                 title="OmniPath" class="full-width" />
        </div>
        <div class="span2 u-align-right">
            <h2>Metabolomics prior knowledge for systems biology</h2>
        </div>
    </div>

    <!-- Navigation -->
    <div class="grid-container u-align-left full">
        <nav>
            <a class="topmenu" href="#cosmos"><span class="nav">COSMOS PKN</span></a>
            <a class="topmenu" href="/schema/swagger"><span class="nav">API docs</span></a>
            <a class="topmenu" href="/cosmos/status"><span class="nav">status</span></a>
            <a class="topmenu" href="https://saezlab.github.io/omnipath-metabo"><span class="nav">documentation</span></a>
            <a class="topmenu" href="https://github.com/saezlab/omnipath-metabo"><span class="nav">GitHub</span></a>
            <a class="topmenu" href="https://omnipathdb.org"><span class="nav">OmniPath</span></a>
            <a class="topmenu" href="https://utils.omnipathdb.org"><span class="nav">Utils</span></a>
        </nav>
    </div>

    <!-- About -->
    <div class="grid-container u-align-left full">
        <div>
            <p>
                OmniPath Metabolomics provides pre-built prior-knowledge networks
                for multi-omics causal reasoning (COSMOS).  The COSMOS PKN combines
                metabolite&ndash;protein interactions from TCDB, SLC, BRENDA, MRCLinksDB,
                STITCH, genome-scale metabolic models (Human-GEM, Recon3D), and more.
                Available as a Python library and as this HTTP API with
                <a href="/schema/swagger">interactive documentation</a>.
            </p>
            <ul>
                <li><strong>COSMOS PKN</strong> &mdash; pre-built prior-knowledge network combining transporters, receptors, allosteric regulators, and enzyme&ndash;metabolite interactions</li>
                <li><strong>ID translation</strong> &mdash; all metabolites unified to ChEBI, all proteins to primary SwissProt ACs</li>
                <li><strong>Multiple resources</strong> &mdash; TCDB, SLC, BRENDA, MRCLinksDB, STITCH, Human-GEM, Recon3D</li>
                <li><strong>Multi-organism</strong> &mdash; human, mouse, and additional organisms (coming soon)</li>
            </ul>
        </div>
    </div>

    <!-- COSMOS PKN -->
    <div class="grid-container u-align-left full" id="cosmos">
        <div>
            <h5>COSMOS PKN API</h5>
            <p>Get the full human transporter PKN:</p>
            <div class="box codebox">
                <a href="/cosmos/pkn?organism=9606&categories=transporters" class="no-uline code">/cosmos/pkn?organism=9606&amp;categories=transporters</a>
            </div>
            <p>Get transporters + receptors:</p>
            <div class="box codebox">
                <a href="/cosmos/pkn?organism=9606&categories=transporters,receptors" class="no-uline code">/cosmos/pkn?organism=9606&amp;categories=transporters,receptors</a>
            </div>
            <p>Filter by resource:</p>
            <div class="box codebox">
                <a href="/cosmos/pkn?organism=9606&categories=transporters&resources=TCDB,SLC" class="no-uline code">/cosmos/pkn?organism=9606&amp;categories=transporters&amp;resources=TCDB,SLC</a>
            </div>
            <p>More endpoints:</p>
            <ul>
                <li><a href="/cosmos/categories">/cosmos/categories</a> &mdash; available PKN categories</li>
                <li><a href="/cosmos/organisms">/cosmos/organisms</a> &mdash; organisms with pre-built PKNs</li>
                <li><a href="/cosmos/resources">/cosmos/resources</a> &mdash; resources per category</li>
                <li><a href="/cosmos/status">/cosmos/status</a> &mdash; cache status and build info</li>
            </ul>
        </div>
    </div>

    <!-- Python client -->
    <div class="grid-container u-align-left full" id="python">
        <div>
            <h5>Python usage</h5>
            <div class="box codebox">
                <pre><code>pip install omnipath-metabo

from omnipath_metabo.datasets.cosmos import build_transporters
bundle = build_transporters()
# 51,067 edges, 3,045 metabolites, 2,317 proteins</code></pre>
            </div>
        </div>
    </div>

    <!-- Footer -->
    <div class="grid-container u-align-left full">
        <p class="small">
            OmniPath Metabolomics is developed by the
            <a href="https://saezlab.org/">Saez Lab</a> at
            Heidelberg University.
            Source code:
            <a href="https://github.com/saezlab/omnipath-metabo">GitHub</a>.
        </p>
    </div>
</body>
</html>
"""
