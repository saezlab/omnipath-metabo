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


def _init_backend():
    """Return (translate_fn, table_fn, mode_name)."""

    # Try direct DB mode first (requires omnipath-utils[db] + reachable PostgreSQL)
    try:
        from omnipath_utils.mapping import translate as _translate
        from omnipath_utils.mapping import translation_table as _table

        # Smoke test: verify both translate and translation_table work
        # against the PostgreSQL DB.  Test with a metabolite pair
        # (organism=0) — if this returns empty, the in-memory backend is
        # being used instead of the DB, which would trigger large upstream
        # downloads during the build.
        _translate(['P04637'], 'uniprot', 'genesymbol', 9606)
        _test_table = _table('pubchem', 'chebi', 0)

        if not _test_table:
            raise RuntimeError('translation_table returned empty — DB not reachable')

        _log.info('[COSMOS] Using omnipath-utils DB mode for ID translation.')
        return _translate, _table, 'database'
    except Exception:
        pass

    # Fall back to HTTP
    from omnipath_client.utils import translate as _translate

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

    _log.info('[COSMOS] Using omnipath-client HTTP mode for ID translation.')
    return _translate, _table, 'http'


mapping_translate, mapping_table, mapping_mode = _init_backend()
