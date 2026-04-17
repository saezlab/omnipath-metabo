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

"""Dual-mode adapter: omnipath-utils DB -> omnipath-client HTTP fallback."""

from __future__ import annotations

import logging

_log = logging.getLogger(__name__)


_HTTP_CHUNK_SIZE = 500


def _chunked_translate(
    _translate_fn,
    identifiers: list,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int = 9606,
) -> dict:
    """
    Call *_translate_fn* in chunks to avoid HTTP read timeouts on large batches.

    The omnipath-client HTTP backend sends everything in a single POST with a
    120-second timeout. Gene symbol lists from Mouse-GEM (5000+) regularly
    exceed this. Splitting into chunks of :data:`_HTTP_CHUNK_SIZE` keeps each
    request well under the limit.
    """

    if len(identifiers) <= _HTTP_CHUNK_SIZE:
        return _translate_fn(identifiers, id_type, target_id_type, ncbi_tax_id)

    result: dict = {}

    for i in range(0, len(identifiers), _HTTP_CHUNK_SIZE):
        chunk = identifiers[i: i + _HTTP_CHUNK_SIZE]
        _log.debug(
            '[COSMOS] HTTP translate chunk %d-%d / %d (%s→%s)',
            i, i + len(chunk), len(identifiers), id_type, target_id_type,
        )
        result.update(_translate_fn(chunk, id_type, target_id_type, ncbi_tax_id))

    return result


def _init_backend():
    """Return (translate_fn, table_fn, mode_name).

    Selects between two backends:

    - **HTTP mode** (default): via ``omnipath-client``, queries the deployed
      service at ``utils.omnipathdb.org``.  This uses the pre-populated
      PostgreSQL database directly and is the right choice for most users.

    - **Database mode** (opt-in via ``COSMOS_MAPPING_MODE=database``): via
      ``omnipath_utils.mapping`` Python API.  Note this uses the in-memory
      backend system, not the PostgreSQL DB — it will trigger upstream
      downloads (UniChem, BioMart, etc.) on first use of each mapping pair.
      Useful for offline development with cached pypath data.
    """

    import os

    mode = os.environ.get('COSMOS_MAPPING_MODE', 'http').lower()

    if mode == 'database':
        try:
            from omnipath_utils.mapping import translate as _translate
            from omnipath_utils.mapping import translation_table as _table

            _translate(['P04637'], 'uniprot', 'genesymbol', 9606)
            _log.info(
                '[COSMOS] Using omnipath-utils in-memory mode '
                '(COSMOS_MAPPING_MODE=database).'
            )
            return _translate, _table, 'database'
        except Exception as exc:
            _log.warning(
                '[COSMOS] Database mode requested but unavailable (%s); '
                'falling back to HTTP.', exc,
            )

    from omnipath_client.utils import translate as _http_translate

    def _translate(
        identifiers: list,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int = 9606,
    ) -> dict:
        return _chunked_translate(
            _http_translate, identifiers, id_type, target_id_type, ncbi_tax_id,
        )

    def _table(id_type: str, target_id_type: str, ncbi_tax_id: int = 0) -> dict:
        """HTTP mode shim: full-table fetches are not supported without DB.

        translation_dict was removed from omnipath_client; the HTTP API only
        supports translating known IDs.  Callers that pre-load entire tables
        (e.g. pubchem→chebi) will receive an empty dict, meaning metabolite ID
        translation will be skipped in HTTP-only mode.  Use omnipath-utils with
        a reachable PostgreSQL instance for full metabolite translation support.
        """
        _log.warning(
            '[COSMOS] HTTP mode: full mapping table (%s→%s) unavailable '
            '(omnipath_client no longer provides translation_dict). '
            'Metabolite ID translation requires DB mode.',
            id_type,
            target_id_type,
        )
        return {}

    _log.info(
        '[COSMOS] Using omnipath-client HTTP mode for ID translation.'
    )
    return _translate, _table, 'http'


mapping_translate, mapping_table, mapping_mode = _init_backend()
