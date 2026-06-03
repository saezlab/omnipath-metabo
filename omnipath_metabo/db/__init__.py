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

"""Main-Postgres access for the metabo post-build chemistry layer."""

from __future__ import annotations

from omnipath_metabo.db._connection import connect, resolve_db_url
from omnipath_metabo.db._schema import (
    STRUCTURAL_SPECIFICITY_LEVELS,
    ensure_rdkit_extension,
    ensure_metabo_schema,
    current_build_id,
    recorded_build_id,
    record_build_id,
)

__all__ = [
    'connect',
    'resolve_db_url',
    'STRUCTURAL_SPECIFICITY_LEVELS',
    'ensure_rdkit_extension',
    'ensure_metabo_schema',
    'current_build_id',
    'recorded_build_id',
    'record_build_id',
]
