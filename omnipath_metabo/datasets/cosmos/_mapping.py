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

    # Try direct DB mode first
    try:
        from omnipath_utils.mapping import translate as _translate
        from omnipath_utils.mapping import translation_table as _table

        # Smoke test: can we reach the DB?
        _translate(['P04637'], 'uniprot', 'genesymbol', 9606)
        _log.info('[COSMOS] Using omnipath-utils DB mode for ID translation.')
        return _translate, _table, 'database'
    except Exception:
        pass

    # Fall back to HTTP
    from omnipath_client.utils import translate as _translate
    from omnipath_client.utils import translation_dict as _table

    _log.info('[COSMOS] Using omnipath-client HTTP mode for ID translation.')
    return _translate, _table, 'http'


mapping_translate, mapping_table, mapping_mode = _init_backend()
