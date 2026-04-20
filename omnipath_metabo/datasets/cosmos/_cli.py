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
Command-line interface for building and exporting the COSMOS PKN.

Entry point: ``cosmos-pkn`` (installed via ``pip install omnipath-metabo``).

Pipeline::

    build() → format_pkn() → write file

Examples::

    # Full human PKN, default output cosmos_pkn.csv
    cosmos-pkn

    # Transporters only, TSV output
    cosmos-pkn --subset transporters --output cosmos_transporters.tsv

    # Mouse, no STITCH, drop connector edges
    cosmos-pkn --organism 10090 --no-stitch --no-connector-edges

    # Full PKN with every column (for inspection)
    cosmos-pkn --all-columns

Output columns (default)
------------------------
``source``  COSMOS-formatted node ID (e.g. ``Metab__CHEBI:15422_c``,
            ``Gene1__P00533``)
``target``  COSMOS-formatted node ID
``sign``    Mode of regulation (1 = activation, -1 = inhibition, 0 = unknown)

Use ``--all-columns`` to also include ``interaction_type``, ``resource``,
``source_type``, ``target_type``, and ``locations``.
"""

from __future__ import annotations

__all__ = ['main']

import argparse
import sys
from pathlib import Path


def _parse_args() -> argparse.Namespace:

    p = argparse.ArgumentParser(
        description='Build and export the COSMOS prior-knowledge network.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest='command')

    # ── export (default, backward-compatible) ─────────────────────────
    export_p = sub.add_parser(
        'export',
        help='Build PKN and export to CSV/TSV.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_export_args(export_p)

    # ── build-cache ───────────────────────────────────────────────────
    cache_p = sub.add_parser(
        'build-cache',
        help='Pre-build PKN categories as Parquet files for the server.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    cache_p.add_argument(
        '--organism', type=int, nargs='+', default=[9606],
        help='NCBI taxonomy ID(s) to build for.',
    )
    cache_p.add_argument(
        '--category', nargs='+', default=None,
        choices=['transporters', 'receptors', 'allosteric', 'enzyme_metabolite'],
        help='Categories to build (default: all four).',
    )
    cache_p.add_argument(
        '--cache-dir', default=None,
        help='Directory for Parquet cache files.',
    )

    args = p.parse_args()

    # When invoked without a subcommand, fall back to ``export``
    # for backward-compatibility (``cosmos-pkn --organism 9606`` still works).
    if args.command is None:
        args = export_p.parse_args(sys.argv[1:])
        args.command = 'export'

    return args


def _add_export_args(p: argparse.ArgumentParser) -> None:
    """Add the export-specific CLI arguments to *p*."""

    # --- organism / resources ------------------------------------------------
    p.add_argument(
        '--organism', type=int, default=9606,
        help='NCBI taxonomy ID.',
    )
    p.add_argument(
        '--subset',
        choices=['all', 'transporters', 'receptors', 'allosteric', 'enzyme_metabolite'],
        default='all',
        help=(
            'Functional subset to build. '
            '"transporters": membrane transport (TCDB, SLC, GEM transporters, '
            'Recon3D, STITCH transporters). '
            '"receptors": ligand-receptor signalling (MRCLinksDB, STITCH receptors). '
            '"allosteric": small-molecule allosteric regulation (BRENDA, STITCH other). '
            '"enzyme_metabolite": stoichiometric enzyme-metabolite reactions from GEMs. '
            '"all": union of all four subsets.'
        ),
    )
    p.add_argument(
        '--score-threshold', type=int, default=700,
        help='STITCH combined confidence score threshold (0–1000).',
    )
    p.add_argument('--no-stitch',     action='store_true', help='Exclude STITCH.')
    p.add_argument('--no-tcdb',       action='store_true', help='Exclude TCDB.')
    p.add_argument('--no-slc',        action='store_true', help='Exclude SLC tables.')
    p.add_argument('--no-brenda',     action='store_true', help='Exclude BRENDA.')
    p.add_argument('--no-mrclinksdb', action='store_true', help='Exclude MRCLinksDB.')
    p.add_argument('--no-gem',        action='store_true', help='Exclude Human-GEM.')
    p.add_argument('--no-recon3d',    action='store_true', help='Exclude Recon3D.')

    # --- formatter options ---------------------------------------------------
    p.add_argument(
        '--no-orphans', action='store_true',
        help=(
            'Drop orphan pseudo-enzyme nodes (reaction_id entries with no '
            'mapped gene product).'
        ),
    )
    p.add_argument(
        '--no-connector-edges', action='store_true',
        help=(
            'Exclude connector edges (bare ChEBI/UniProt ID → COSMOS node). '
            'Connector edges are required by the COSMOS R package but can be '
            'omitted for manual inspection.'
        ),
    )

    # --- output --------------------------------------------------------------
    p.add_argument(
        '--output', default='cosmos_pkn.csv',
        help=(
            'Output file path. Extension determines the separator: '
            '.tsv/.txt → tab-separated; everything else → comma-separated.'
        ),
    )
    p.add_argument(
        '--sep', default=None,
        help='Override the column separator (e.g. "\\t" for TSV regardless of extension).',
    )
    p.add_argument(
        '--all-columns', action='store_true',
        help=(
            'Write all PKN columns (source, target, sign, interaction_type, '
            'resource, source_type, target_type, locations) instead of the '
            'minimal three.'
        ),
    )


def _build(args: argparse.Namespace):
    """Call the appropriate build function and return a CosmosBundle."""

    from omnipath_metabo.datasets.cosmos import (
        build,
        build_allosteric,
        build_enzyme_metabolite,
        build_receptors,
        build_transporters,
    )

    kwargs: dict = dict(
        organism=args.organism,
        stitch={'score_threshold': args.score_threshold},
    )

    if args.no_stitch:
        kwargs['stitch'] = False
    if args.no_tcdb:
        kwargs['tcdb'] = False
    if args.no_slc:
        kwargs['slc'] = False
    if args.no_brenda:
        kwargs['brenda'] = False
    if args.no_mrclinksdb:
        kwargs['mrclinksdb'] = False
    if args.no_gem:
        kwargs['gem'] = False
    if args.no_recon3d:
        kwargs['recon3d'] = False

    fn = {
        'all':              build,
        'transporters':     build_transporters,
        'receptors':        build_receptors,
        'allosteric':       build_allosteric,
        'enzyme_metabolite': build_enzyme_metabolite,
    }[args.subset]

    return fn(**kwargs)


def _separator(args: argparse.Namespace) -> str:

    if args.sep is not None:
        return args.sep.encode('raw_unicode_escape').decode('unicode_escape')

    suffix = Path(args.output).suffix.lower()

    if suffix in {'.tsv', '.txt'}:
        return '\t'

    return ','


def main() -> None:
    """Entry point for the ``cosmos-pkn`` command."""

    args = _parse_args()

    if args.command == 'build-cache':
        _run_build_cache(args)
        return

    _run_export(args)


def _run_build_cache(args: argparse.Namespace) -> None:
    """Execute the ``build-cache`` subcommand."""
    import logging

    logging.basicConfig(level=logging.INFO, format='%(message)s')

    from omnipath_metabo.datasets.cosmos._cache import build_cache

    build_cache(
        organisms=args.organism,
        categories=args.category,
        cache_dir=args.cache_dir,
    )


def _run_export(args: argparse.Namespace) -> None:
    """Execute the ``export`` subcommand (default)."""

    import pandas as pd

    from omnipath_metabo.datasets.cosmos._format import format_pkn
    from omnipath_metabo.datasets.cosmos._record import CosmosEdge

    # ------------------------------------------------------------------
    # 1. Build translated PKN
    # ------------------------------------------------------------------
    print(f'Building COSMOS PKN (organism={args.organism}, subset={args.subset})...')
    bundle = _build(args)
    print(f'  Translated edges: {len(bundle.network):,}')

    # ------------------------------------------------------------------
    # 2. Format node IDs
    # ------------------------------------------------------------------
    print('Applying COSMOS node-ID formatting...')
    formatted_bundle = format_pkn(bundle, include_orphans=not args.no_orphans)
    formatted = pd.DataFrame(formatted_bundle.network)

    if formatted.empty:
        print('  No edges produced — check build options.')
        sys.exit(0)

    n_connector = (formatted['interaction_type'] == 'connector').sum()
    n_main = len(formatted) - n_connector
    print(f'  Main edges:      {n_main:,}')
    print(f'  Connector edges: {n_connector:,}')
    print(f'  Total:           {len(formatted):,}')

    # ------------------------------------------------------------------
    # 3. Optionally drop connector edges
    # ------------------------------------------------------------------
    if args.no_connector_edges:
        formatted = formatted[formatted['interaction_type'] != 'connector'].reset_index(drop=True)
        print(f'  After dropping connectors: {len(formatted):,} edges')

    # ------------------------------------------------------------------
    # 4. Select output columns and rename mor → sign
    # ------------------------------------------------------------------
    formatted['sign'] = formatted['mor']

    if args.all_columns:
        out_cols = [
            'source', 'target', 'sign',
            'interaction_type', 'resource',
            'source_type', 'target_type',
            'locations',
        ]
    else:
        out_cols = ['source', 'target', 'sign']

    out_cols = [c for c in out_cols if c in formatted.columns]
    output_df = formatted[out_cols]

    # ------------------------------------------------------------------
    # 5. Write output
    # ------------------------------------------------------------------
    sep = _separator(args)
    out_path = Path(args.output)
    output_df.to_csv(out_path, sep=sep, index=False)
    print(f'\nSaved to {out_path}  ({len(output_df):,} rows, {len(out_cols)} columns)')

    # ------------------------------------------------------------------
    # 6. Summary by resource (main edges only)
    # ------------------------------------------------------------------
    if 'resource' in formatted.columns and 'interaction_type' in formatted.columns:
        main = formatted[formatted['interaction_type'] != 'connector']
        counts = main.groupby('resource').size().sort_values(ascending=False)
        print('\nEdge counts by resource (main edges only):')
        for resource, n in counts.items():
            print(f'  {resource:<35} {n:>7,}')
        print(f'  {"TOTAL":<35} {counts.sum():>7,}')
