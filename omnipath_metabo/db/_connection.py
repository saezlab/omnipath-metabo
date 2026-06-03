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
Connection to the MAIN OmniPath Postgres for the metabo post-build step.

The metabo chemistry layer writes only ``metabo_*``-prefixed tables in the shared
``public`` schema; in production it should connect as a metabo-scoped role whose
grants cover those tables plus read access to ``entity`` / ``facet_entity_bitmap``
/ ``identifier_evidence`` / ``build_manifest``.
"""

from __future__ import annotations

__all__ = ['connect', 'resolve_db_url']

import os


def resolve_db_url(db_url: str | None = None) -> str:
    """Resolve the main-DB URL from the argument or environment."""
    url = (
        db_url
        or os.environ.get('OMNIPATH_DB_URL')
        or os.environ.get('DATABASE_URL')
    )
    if not url:
        raise RuntimeError(
            'No main-database URL; pass --db-url or set OMNIPATH_DB_URL '
            '(or DATABASE_URL).'
        )
    return url


def connect(db_url: str | None = None):
    """Open a psycopg2 connection to the main OmniPath Postgres."""
    import psycopg2

    return psycopg2.connect(resolve_db_url(db_url))
