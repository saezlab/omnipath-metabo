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
OmniPath protein-protein and gene regulatory interactions.

Queries the legacy OmniPath web service (omnipathdb.org) for
protein-protein signaling (including ligand-receptor and kinase-
substrate) and transcription factor-target gene interactions.
These complement the metabolite-protein interactions from the other
COSMOS resources.

This is a temporary bridge module — once the new OmniPath database
build (omnipath-server) is operational, these interactions will be
sourced directly from the new database.
"""

from __future__ import annotations

__all__ = [
    'ppi_interactions',
    'grn_interactions',
]

import csv
import io
import logging
from collections.abc import Generator

from .._record import Interaction

_log = logging.getLogger(__name__)

_OMNIPATH_BASE = 'https://omnipathdb.org'


def _query_omnipath(
    endpoint: str = 'interactions',
    datasets: str | None = None,
    organism: int = 9606,
    resources: str | None = None,
    **extra_params: str,
) -> list[dict]:
    """Query the OmniPath web API and return rows as dicts.

    Args:
        endpoint: API endpoint (``'interactions'`` or ``'enzsub'``).
        datasets: Comma-separated dataset names (for ``/interactions``).
        organism: NCBI taxonomy ID.
        resources: Comma-separated resource filter (e.g. ``'SIGNOR'``).
        **extra_params: Additional query parameters.

    Returns:
        List of row dicts parsed from TSV response.
    """

    params = {
        'genesymbols': 'yes',
        'fields': 'sources,references',
        'organisms': str(organism),
    }

    if datasets:
        params['datasets'] = datasets

    if resources:
        params['resources'] = resources

    params.update(extra_params)

    query = '&'.join(f'{k}={v}' for k, v in params.items())
    url = f'{_OMNIPATH_BASE}/{endpoint}?{query}'

    _log.info('[COSMOS] Querying OmniPath: %s', url)

    try:
        from dlmachine import Download

        dl = Download(url)
        dl.download()
        text = dl.result.text

    except (ImportError, Exception):
        import ssl
        import urllib.request

        try:
            import certifi
            ctx = ssl.create_default_context(cafile=certifi.where())
        except ImportError:
            _log.warning('[COSMOS] certifi not found; disabling SSL verification for OmniPath request')
            ctx = ssl._create_unverified_context()

        with urllib.request.urlopen(url, timeout=120, context=ctx) as resp:
            text = resp.read().decode('utf-8')

    reader = csv.DictReader(io.StringIO(text), delimiter='\t')
    rows = list(reader)

    _log.info(
        '[COSMOS] OmniPath %s/%s: %d interactions',
        endpoint, datasets or '', len(rows),
    )

    return rows


def _sign_to_mor(row: dict) -> int:
    """Convert OmniPath stimulation/inhibition to MOR (mode of regulation).

    Accepts both '1'/'0' and 'True'/'False' string values — the legacy
    omnipathdb.org API returns 'True'/'False', not '1'/'0'.
    """

    def _flag(key: str) -> bool:
        return str(row.get(key, '0')).lower() in ('1', 'true')

    stim = _flag('is_stimulation')
    inhib = _flag('is_inhibition')

    if stim and not inhib:
        return 1
    if inhib and not stim:
        return -1

    return 0  # unknown or contradictory


def ppi_interactions(
    organism: int = 9606,
    datasets: str = 'omnipath,ligrecextra',
    resources: str | None = None,
) -> Generator[Interaction, None, None]:
    """Protein-protein interactions from OmniPath.

    Yields directed edges from the combined OmniPath datasets.
    All protein-protein interactions (signaling, ligand-receptor,
    kinase-substrate) are queried together to avoid duplicate edges
    from overlapping datasets.

    Args:
        organism: NCBI taxonomy ID.
        datasets: Comma-separated OmniPath dataset names.
            Default ``'omnipath,ligrecextra'``.
            Extended: ``'omnipath,ligrecextra,kinaseextra,pathwayextra'``.
        resources: Comma-separated resource filter (e.g. ``'SIGNOR'``).

    Yields:
        :class:`~.Interaction` records.
    """

    rows = _query_omnipath(datasets=datasets, organism=organism, resources=resources)

    for row in rows:
        source = row.get('source_genesymbol', '')
        target = row.get('target_genesymbol', '')

        if not source or not target:
            continue

        yield Interaction(
            source=source,
            target=target,
            source_type='protein',
            target_type='protein',
            id_type_a='genesymbol',
            id_type_b='genesymbol',
            interaction_type='signaling',
            resource=f'OmniPath:{datasets}',
            mor=_sign_to_mor(row),
            locations=(),
            attrs={
                'sources': row.get('sources', ''),
                'references': row.get('references', ''),
                'directed': row.get('is_directed', '0') == '1',
            },
        )


def grn_interactions(
    organism: int = 9606,
    datasets: str = 'collectri',
    resources: str | None = None,
    dorothea_levels: str | None = None,
) -> Generator[Interaction, None, None]:
    """Gene regulatory network from OmniPath.

    Yields TF->gene directed edges from ``collectri`` by default.
    Optionally includes ``dorothea`` with confidence level filtering.

    Args:
        organism: NCBI taxonomy ID.
        datasets: Comma-separated dataset names.
            Default ``'collectri'``.
            Use ``'collectri,dorothea'`` for combined.
        resources: Comma-separated resource filter.
        dorothea_levels: Confidence levels for DoRothEA (e.g. ``'A,B,C'``).
            Only applied when ``'dorothea'`` is in *datasets*.
            Default: ``'A,B,C'``.

    Yields:
        :class:`~.Interaction` records.
    """

    extra_params: dict[str, str] = {}

    if 'dorothea' in datasets and dorothea_levels is not None:
        extra_params['dorothea_levels'] = dorothea_levels
    elif 'dorothea' in datasets:
        extra_params['dorothea_levels'] = 'A,B,C'

    rows = _query_omnipath(
        datasets=datasets, organism=organism,
        resources=resources, **extra_params,
    )

    for row in rows:
        source = row.get('source_genesymbol', '')
        target = row.get('target_genesymbol', '')

        if not source or not target:
            continue

        yield Interaction(
            source=source,
            target=target,
            source_type='protein',
            target_type='protein',
            id_type_a='genesymbol',
            id_type_b='genesymbol',
            interaction_type='gene_regulation',
            resource=f'OmniPath:{datasets}',
            mor=_sign_to_mor(row),
            locations=(),
            attrs={
                'sources': row.get('sources', ''),
                'references': row.get('references', ''),
                'directed': row.get('is_directed', '0') == '1',
            },
        )
