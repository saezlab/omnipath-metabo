#!/usr/bin/env python

"""
Build and export the COSMOS prior-knowledge network to CSV / TSV.

Pipeline::

    build() → format_pkn() → write file

Usage::

    python scripts/build_cosmos_pkn.py

    # Human PKN, transporters only, TSV output
    python scripts/build_cosmos_pkn.py \\
        --subset transporters --output cosmos_transporters.tsv

    # Mouse, no STITCH, drop connector edges
    python scripts/build_cosmos_pkn.py \\
        --organism 10090 --no-stitch --no-connector-edges

    # Full PKN with every column (for debugging)
    python scripts/build_cosmos_pkn.py --all-columns

Output columns (default)
-------------------------
``source``  Formatted COSMOS node ID (``Metab__CHEBI:...`` or ``Gene{N}__ENSG...``)
``target``  Formatted COSMOS node ID
``sign``    Mode of regulation (1 = activation/positive, -1 = inhibition/negative)

Use ``--all-columns`` to also include ``interaction_type``, ``resource``,
``locations``, and ``attrs``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='Build and export the COSMOS PKN.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

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
            'Which functional subset to build.  '
            'All subsets contain protein-metabolite interactions; '
            'they differ by biological function: '
            '"transporters" = membrane transport (TCDB, SLC, GEM transporters, Recon3D, STITCH transporters); '
            '"receptors" = ligand-receptor signalling (MRCLinksDB, STITCH receptors); '
            '"allosteric" = allosteric regulation, small molecules that activate/inhibit proteins '
            '(BRENDA, STITCH other); '
            '"enzyme_metabolite" = stoichiometric enzyme-metabolite reactions from GEMs (GEM metabolic only).  '
            '"all" is the union of all four.'
        ),
    )
    p.add_argument('--score-threshold', type=int, default=700,
                   help='STITCH combined score threshold.')
    p.add_argument('--no-stitch', action='store_true')
    p.add_argument('--no-tcdb', action='store_true')
    p.add_argument('--no-slc', action='store_true')
    p.add_argument('--no-brenda', action='store_true')
    p.add_argument('--no-mrclinksdb', action='store_true')
    p.add_argument('--no-gem', action='store_true')
    p.add_argument('--no-recon3d', action='store_true')

    # --- formatter options ---------------------------------------------------
    p.add_argument(
        '--no-orphans', action='store_true',
        help='Drop orphan pseudo-enzyme nodes (reaction_id entries with no gene).',
    )
    p.add_argument(
        '--no-connector-edges', action='store_true',
        help=(
            'Exclude connector edges (bare ID → formatted node) from the output. '
            'Connector edges link raw ChEBI / ENSG IDs to the COSMOS-formatted '
            'node IDs; they are required by the COSMOS R package but can be '
            'omitted for inspection.'
        ),
    )

    # --- output --------------------------------------------------------------
    p.add_argument(
        '--output', default='cosmos_pkn.csv',
        help=(
            'Output file path.  Extension determines the separator: '
            '.tsv / .txt → tab-separated; everything else → comma-separated.'
        ),
    )
    p.add_argument(
        '--sep', default=None,
        help='Override the column separator (e.g. "\\t" for TSV).',
    )
    p.add_argument(
        '--all-columns', action='store_true',
        help=(
            'Write all PKN columns (source, target, sign, interaction_type, '
            'resource, locations, source_type, target_type) instead of just '
            'the minimal three.'
        ),
    )

    return p.parse_args()


def _build(args: argparse.Namespace):
    """Call the appropriate build function."""
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
        'all': build,
        'transporters': build_transporters,
        'receptors': build_receptors,
        'allosteric': build_allosteric,
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
    import pandas as pd
    from omnipath_metabo.datasets.cosmos._format import format_pkn

    args = _parse_args()

    # ------------------------------------------------------------------
    # 1. Build translated PKN
    # ------------------------------------------------------------------
    print(f'Building COSMOS PKN (organism={args.organism}, subset={args.subset})...')
    df = _build(args)
    print(f'  Translated PKN: {len(df):,} edges')

    # ------------------------------------------------------------------
    # 2. Format node IDs
    # ------------------------------------------------------------------
    print('Applying COSMOS node-ID formatting...')
    formatted = format_pkn(df, include_orphans=not args.no_orphans)

    n_connector = (formatted['interaction_type'] == 'connector').sum()
    n_main = len(formatted) - n_connector
    print(f'  Main edges:      {n_main:,}')
    print(f'  Connector edges: {n_connector:,}')
    print(f'  Total:           {len(formatted):,}')

    # ------------------------------------------------------------------
    # 3. Optionally drop connector edges
    # ------------------------------------------------------------------
    if args.no_connector_edges:
        formatted = formatted[formatted['interaction_type'] != 'connector'].copy()
        print(f'  After dropping connectors: {len(formatted):,} edges')

    # ------------------------------------------------------------------
    # 4. Select output columns and rename mor → sign
    # ------------------------------------------------------------------
    formatted = formatted.copy()
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
    # 6. Summary by resource and interaction type
    # ------------------------------------------------------------------
    if 'resource' in formatted.columns and 'interaction_type' in formatted.columns:
        print('\nEdge counts by resource (main edges only):')
        main = formatted[formatted['interaction_type'] != 'connector']
        counts = (
            main.groupby('resource')
            .size()
            .sort_values(ascending=False)
        )
        for resource, n in counts.items():
            print(f'  {resource:<35} {n:>7,}')
        print(f'  {"TOTAL":<35} {counts.sum():>7,}')


if __name__ == '__main__':
    main()
