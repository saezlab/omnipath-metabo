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
Server launcher for the omnipath-metabo web service.

Wraps uvicorn to serve the Litestar app from
:mod:`omnipath_metabo.server._app`.
"""

from __future__ import annotations

__all__ = ['run_server']

import argparse


def run_server(args: argparse.Namespace) -> None:
    """Start the uvicorn server with the Litestar app.

    Args:
        args: Parsed CLI arguments with ``cache_dir``, ``host``,
            ``port``, and ``reload`` attributes.
    """
    import os

    import uvicorn

    # Pass cache_dir via environment so the app factory can pick it up.
    if args.cache_dir:
        os.environ['OMNIPATH_METABO_CACHE_DIR'] = args.cache_dir

    uvicorn.run(
        'omnipath_metabo.server._app:_app_instance',
        factory=False,
        host=args.host,
        port=args.port,
        reload=args.reload,
    )
