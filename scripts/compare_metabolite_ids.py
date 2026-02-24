#!/usr/bin/env python

"""
HMDB vs ChEBI metabolite ID coverage comparison.

Builds the raw (untranslated) PKN for each resource and attempts both
ChEBI and HMDB translation for every metabolite ID.  Reports per-resource
and overall success rates so you can make a data-driven decision about
whether switching the default target ID type to HMDB would reduce
information loss.

Usage::

    uv run python scripts/compare_metabolite_ids.py

    # Optional: restrict resources or organisms
    uv run python scripts/compare_metabolite_ids.py --no-brenda --organism 10090

Output:
    Printed summary table + saved as ``metabolite_id_coverage.csv`` in the
    current working directory.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict

import pandas as pd


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='Compare ChEBI vs HMDB metabolite ID coverage across PKN resources.',
    )
    p.add_argument('--organism', type=int, default=9606,
                   help='NCBI taxonomy ID (default: 9606)')
    p.add_argument('--score-threshold', type=int, default=700,
                   help='STITCH combined score threshold (default: 700)')
    p.add_argument('--no-stitch', action='store_true')
    p.add_argument('--no-tcdb', action='store_true')
    p.add_argument('--no-slc', action='store_true')
    p.add_argument('--no-brenda', action='store_true')
    p.add_argument('--no-mrclinksdb', action='store_true')
    p.add_argument('--no-gem', action='store_true')
    p.add_argument('--no-recon3d', action='store_true')
    p.add_argument('--output', default='metabolite_id_coverage.csv',
                   help='Output CSV path (default: metabolite_id_coverage.csv)')
    return p.parse_args()


def _build_raw(args: argparse.Namespace) -> pd.DataFrame:
    """Build an untranslated PKN DataFrame."""
    from omnipath_metabo.datasets.cosmos import build

    kwargs: dict = dict(
        organism=args.organism,
        translate_ids=False,
        apply_blacklist=False,
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

    return build(**kwargs)


def _extract_metabolite_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return rows where the *source* is a small molecule (the typical case).

    Also handles GEM enzyme→metabolite edges where the target is the
    small molecule.
    """
    source_met = df[df['source_type'] == 'small_molecule'].copy()
    source_met = source_met.rename(columns={
        'source': 'met_id',
        'id_type_a': 'id_type',
    })[['met_id', 'id_type', 'resource']]

    target_met = df[df['target_type'] == 'small_molecule'].copy()
    target_met = target_met.rename(columns={
        'target': 'met_id',
        'id_type_b': 'id_type',
    })[['met_id', 'id_type', 'resource']]

    return pd.concat([source_met, target_met], ignore_index=True)


def _gem_name_from_resource(resource: str) -> str:
    if resource.startswith('GEM') and ':' in resource:
        return resource.split(':', 1)[1]
    return ''


def _try_translate(met_id: str, id_type: str, resource: str) -> tuple[str | None, str | None]:
    """
    Return (chebi_result, hmdb_result) for a single metabolite ID.

    Returns ``None`` for each target type if translation fails.
    """
    from omnipath_metabo.datasets.cosmos._translate import _to_chebi, _to_hmdb

    gem = _gem_name_from_resource(resource)
    chebi = _to_chebi(met_id, id_type, gem=gem)
    hmdb = _to_hmdb(met_id, id_type, gem=gem)
    return chebi, hmdb


def _compare_resource(group: pd.DataFrame) -> dict:
    """Compute coverage statistics for one resource group."""
    total = len(group)
    chebi_ok = 0
    hmdb_ok = 0
    chebi_only = 0
    hmdb_only = 0
    both_ok = 0

    # Deduplicate: each unique (met_id, id_type) pair counted once per resource
    seen: set[tuple] = set()

    for _, row in group.iterrows():
        key = (row['met_id'], row['id_type'])
        if key in seen:
            continue
        seen.add(key)

        chebi, hmdb = _try_translate(row['met_id'], row['id_type'], row['resource'])
        c = chebi is not None
        h = hmdb is not None

        if c:
            chebi_ok += 1
        if h:
            hmdb_ok += 1
        if c and not h:
            chebi_only += 1
        if h and not c:
            hmdb_only += 1
        if c and h:
            both_ok += 1

    unique = len(seen)
    return {
        'total_edges': total,
        'unique_metabolites': unique,
        'chebi_success': chebi_ok,
        'hmdb_success': hmdb_ok,
        'chebi_pct': round(100 * chebi_ok / unique, 1) if unique else 0.0,
        'hmdb_pct': round(100 * hmdb_ok / unique, 1) if unique else 0.0,
        'chebi_only': chebi_only,
        'hmdb_only': hmdb_only,
        'both': both_ok,
    }


def _print_table(rows: list[dict], resources: list[str]) -> None:
    """Pretty-print the comparison table."""
    header = (
        f"{'Resource':<20} | {'Edges':>7} | {'Unique':>7} | "
        f"{'ChEBI':>10} | {'HMDB':>10} | "
        f"{'ChEBI-only':>11} | {'HMDB-only':>10} | {'Both':>7}"
    )
    sep = '-' * len(header)
    print('\n' + sep)
    print(header)
    print(sep)
    for r, row in zip(resources, rows):
        print(
            f"{r:<20} | {row['total_edges']:>7} | {row['unique_metabolites']:>7} | "
            f"{row['chebi_success']:>5} ({row['chebi_pct']:>4.1f}%) | "
            f"{row['hmdb_success']:>5} ({row['hmdb_pct']:>4.1f}%) | "
            f"{row['chebi_only']:>11} | {row['hmdb_only']:>10} | {row['both']:>7}"
        )
    print(sep + '\n')


def main() -> None:
    args = _parse_args()

    print(f'Building raw PKN (organism={args.organism}, translate_ids=False)...')
    df = _build_raw(args)
    print(f'  Total edges: {len(df)}')

    met_df = _extract_metabolite_rows(df)
    print(f'  Metabolite-side rows: {len(met_df)}')

    # Normalise resource names: 'GEM:Human-GEM' → 'GEM:Human-GEM' (keep as-is)
    # grouped by the full resource string so GEM variants are kept separate
    resources = sorted(met_df['resource'].unique())
    result_rows = []

    for resource in resources:
        print(f'  Translating {resource}...')
        group = met_df[met_df['resource'] == resource]
        stats = _compare_resource(group)
        result_rows.append({'resource': resource, **stats})

    # Overall row
    overall_group = met_df.copy()
    overall_stats = _compare_resource(overall_group)
    result_rows.append({'resource': 'TOTAL', **overall_stats})
    resources_for_print = resources + ['TOTAL']

    _print_table(result_rows, resources_for_print)

    out_df = pd.DataFrame(result_rows)
    out_df.to_csv(args.output, index=False)
    print(f'Saved to {args.output}')


if __name__ == '__main__':
    main()
