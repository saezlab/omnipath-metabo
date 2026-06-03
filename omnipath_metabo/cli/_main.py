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
Main CLI entry point for omnipath-metabo.

Usage::

    omnipath-metabo serve --port 8084
    omnipath-metabo serve --cache-dir /path/to/cache
"""

from __future__ import annotations

__all__ = ['main']

import argparse
import sys


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog='omnipath-metabo',
        description='omnipath-metabo command-line interface.',
    )
    sub = p.add_subparsers(dest='command')

    # ── serve ─────────────────────────────────────────────────────────
    serve_p = sub.add_parser(
        'serve',
        help='Start the Litestar web server.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    serve_p.add_argument(
        '--cache-dir', default=None,
        help='Path to the Parquet cache directory.',
    )
    serve_p.add_argument(
        '--host', default='127.0.0.1',
        help='Bind address.',
    )
    serve_p.add_argument(
        '--port', type=int, default=8084,
        help='Bind port.',
    )
    serve_p.add_argument(
        '--reload', action='store_true',
        help='Enable auto-reload (development).',
    )

    # ── post-build-metabo ─────────────────────────────────────────────
    pb = sub.add_parser(
        'post-build-metabo',
        help='Build the in-DB chemistry layer (RDKit cartridge) on the main DB.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    pb.add_argument(
        '--db-url', default=None,
        help='Main OmniPath Postgres URL (else OMNIPATH_DB_URL / DATABASE_URL).',
    )
    pb.add_argument(
        '--schema', default='public',
        help='Schema holding the built OmniPath tables.',
    )
    pb.add_argument(
        '--force', action='store_true',
        help='Recompute even if the main build id changed since the last run.',
    )
    pb.add_argument(
        '--no-conflicts', action='store_true',
        help='Skip the RaMP multi-InChIKey conflict classifier (Milestone F).',
    )
    pb.add_argument(
        '--conflict-max-records', type=int, default=None,
        help='Cap RaMP rows fetched for the conflict classifier (dev iteration).',
    )

    args = p.parse_args()
    if args.command is None:
        p.print_help()
        sys.exit(1)

    return args


def main() -> None:
    """Entry point for the ``omnipath-metabo`` command."""
    args = _parse_args()

    if args.command == 'serve':
        from ._serve import run_server
        run_server(args)
    elif args.command == 'post-build-metabo':
        from omnipath_metabo.postbuild import post_build_metabo
        post_build_metabo(
            db_url=args.db_url,
            schema=args.schema,
            force=args.force,
            conflicts=not args.no_conflicts,
            conflict_max_records=args.conflict_max_records,
        )
