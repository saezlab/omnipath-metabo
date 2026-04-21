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
OmniPath signaling, ligand-receptor, and gene regulatory interactions.

Queries the legacy OmniPath web service (omnipathdb.org) for
protein-protein signaling, ligand-receptor, and transcription
factor-target gene interactions.  These complement the metabolite-
protein interactions from the other COSMOS resources.

This is a temporary bridge module — once the new OmniPath database
build (omnipath-server) is operational, these interactions will be
sourced directly from the new database.
"""

from __future__ import annotations

__all__ = [
    'signaling_interactions',
    'ligrec_interactions',
    'grn_interactions',
]

import csv
import io
import logging
from collections.abc import Generator

from .._record import Interaction

_log = logging.getLogger(__name__)

_OMNIPATH_API = 'https://omnipathdb.org/interactions'


def _query_omnipath(
    datasets: str,
    organism: int = 9606,
    **extra_params: str,
) -> list[dict]:
    """Query the OmniPath web API and return rows as dicts.

    Args:
        datasets: Comma-separated dataset names.
        organism: NCBI taxonomy ID.
        **extra_params: Additional query parameters.

    Returns:
        List of row dicts parsed from TSV response.
    """

    params = {
        'datasets': datasets,
        'genesymbols': 'yes',
        'fields': 'sources,references',
        'organisms': str(organism),
    }
    params.update(extra_params)

    query = '&'.join(f'{k}={v}' for k, v in params.items())
    url = f'{_OMNIPATH_API}?{query}'

    _log.info('[COSMOS] Querying OmniPath: %s', url)

    try:
        from dlmachine import Download

        dl = Download(url)
        dl.download()
        text = dl.result.text

    except (ImportError, Exception):
        # Fallback to urllib
        import urllib.request

        with urllib.request.urlopen(url, timeout=120) as resp:
            text = resp.read().decode('utf-8')

    reader = csv.DictReader(io.StringIO(text), delimiter='\t')
    rows = list(reader)

    _log.info(
        '[COSMOS] OmniPath %s: %d interactions', datasets, len(rows),
    )

    return rows


def _sign_to_mor(row: dict) -> int:
    """Convert OmniPath stimulation/inhibition to MOR (mode of regulation)."""

    stim = row.get('is_stimulation', '0') == '1'
    inhib = row.get('is_inhibition', '0') == '1'

    if stim and not inhib:
        return 1
    if inhib and not stim:
        return -1

    return 0  # unknown or contradictory


def signaling_interactions(
    organism: int = 9606,
    datasets: str = 'omnipath',
) -> Generator[Interaction, None, None]:
    """Protein-protein signaling interactions from OmniPath.

    Yields directed edges from the ``omnipath`` dataset (core PPI and
    signaling). Optionally includes ``kinaseextra`` and/or
    ``pathwayextra`` datasets.

    Args:
        organism: NCBI taxonomy ID.
        datasets: Comma-separated OmniPath dataset names.
            Default ``'omnipath'``.
            Use ``'omnipath,kinaseextra,pathwayextra'`` for extended.

    Yields:
        :class:`~.Interaction` records.
    """

    rows = _query_omnipath(datasets, organism=organism)

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


def ligrec_interactions(
    organism: int = 9606,
    datasets: str = 'ligrecextra',
) -> Generator[Interaction, None, None]:
    """Ligand-receptor interactions from OmniPath.

    Yields directed edges from the ``ligrecextra`` dataset.

    Args:
        organism: NCBI taxonomy ID.
        datasets: Comma-separated dataset names. Default ``'ligrecextra'``.

    Yields:
        :class:`~.Interaction` records.
    """

    rows = _query_omnipath(datasets, organism=organism)

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
            interaction_type='ligand_receptor',
            resource=f'OmniPath:{datasets}',
            mor=_sign_to_mor(row),
            locations=(),
            attrs={
                'sources': row.get('sources', ''),
                'references': row.get('references', ''),
            },
        )


def grn_interactions(
    organism: int = 9606,
    datasets: str = 'collectri',
    dorothea_levels: str | None = None,
) -> Generator[Interaction, None, None]:
    """Gene regulatory network from OmniPath.

    Yields TF→gene directed edges from ``collectri`` by default.
    Optionally includes ``dorothea`` with confidence level filtering.

    Args:
        organism: NCBI taxonomy ID.
        datasets: Comma-separated dataset names.
            Default ``'collectri'``.
            Use ``'collectri,dorothea'`` for combined.
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

    rows = _query_omnipath(datasets, organism=organism, **extra_params)

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
