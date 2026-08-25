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
Query layer for ``/sets/metsigdb``.

Selects rows from the materialized MetSigDB membership table using the five
shared filter columns of the API contract, and nothing else. It owns which rows
come back, in which order, and how many. It owns nothing about their shape.
"""

from __future__ import annotations

__all__ = [
    'DEFAULT_LIMIT',
    'MAX_LIMIT',
    'MetSigDBQuery',
    'build_query',
    'fetch',
]

from dataclasses import dataclass
from typing import Any

TABLE = 'public.metsigdb_membership'

# The same page bounds the network routes of this service already serve.
DEFAULT_LIMIT = 1000
MAX_LIMIT = 100_000

# The published row, in contract order. The query selects these columns by
# name rather than with a star, so a column added to the substrate does not
# silently widen the public response.
COLUMNS: tuple[str, ...] = (
    'metabolite_entity_id',
    'metabolite_label',
    'metabolite_entity_type',
    'inchikey',
    'smiles',
    'hmdb',
    'pubchem',
    'chebi',
    'kegg',
    'resource',
    'set_source_id',
    'set_label',
    'set_type',
    'organism',
    'set_size',
    'set_context',
    'provenance_source',
    'provenance_record',
    'build_id',
)

# Row identity, which is also the primary key, so the order is stable across
# pages and the index already serves it.
ORDER_BY = 'resource, set_source_id, metabolite_entity_id'


@dataclass(frozen=True)
class MetSigDBQuery:
    """One request's filters and page.

    Every field is a published shared filter. Adding a field here widens the
    public API, so the contract has to move first.
    """

    resource: tuple[str, ...] = ()
    set_type: tuple[str, ...] = ()
    organism: int | None = None
    set_source_id: str | None = None
    metabolite_entity_id: str | None = None
    limit: int = DEFAULT_LIMIT
    offset: int = 0


def build_query(spec: MetSigDBQuery) -> tuple[str, dict[str, Any]]:
    """The SQL and bind parameters for one request.

    Separate from execution so the shape of a filter set can be asserted
    without a database.
    """
    where: list[str] = []
    params: dict[str, Any] = {'limit': spec.limit, 'offset': spec.offset}

    # Multi-valued filters. ANY over an array keeps one bind parameter whatever
    # the caller asks for, so the plan does not change with the value count.
    for field in ('resource', 'set_type'):
        values = getattr(spec, field)
        if values:
            where.append(f'{field} = ANY(%({field})s)')
            params[field] = list(values)

    # Scalar filters. `organism` is equality, so a null-organism row never
    # matches: it stays in the substrate and out of an organism-filtered
    # response, which is what the contract says.
    for field in ('organism', 'set_source_id', 'metabolite_entity_id'):
        value = getattr(spec, field)
        if value is not None:
            where.append(f'{field} = %({field})s')
            params[field] = value

    clause = f'WHERE {" AND ".join(where)} ' if where else ''
    columns = ', '.join(COLUMNS)

    return (
        f'SELECT {columns} FROM {TABLE} {clause}'
        f'ORDER BY {ORDER_BY} LIMIT %(limit)s OFFSET %(offset)s'
    ), params


def fetch(conn, spec: MetSigDBQuery) -> list[dict[str, Any]]:
    """Run one query and return its rows as mappings.

    Reads the materialized table and nothing else. No upstream retrieval, and
    no identifier remapping.
    """
    from psycopg2.extras import RealDictCursor

    sql, params = build_query(spec)
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(sql, params)
        return [dict(row) for row in cur.fetchall()]
