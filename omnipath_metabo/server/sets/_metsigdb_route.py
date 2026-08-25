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
Route layer for ``/sets/metsigdb``.

Validates the request, bounds the page, and dispatches to the query and
projection layers. It owns neither which rows come back nor their shape.

Read-only and additive over the ``metsigdb_membership`` table that
``omnipath-build`` publishes. An absent substrate fails clearly with 503.
"""

from __future__ import annotations

__all__ = ['MetSigDBController']

import os
from typing import Any

from litestar import Controller, Request, get
from litestar.exceptions import HTTPException
from litestar.params import Parameter

from omnipath_metabo.server.sets._metsigdb_projection import project_rows
from omnipath_metabo.server.sets._metsigdb_query import (
    DEFAULT_LIMIT,
    MAX_LIMIT,
    MetSigDBQuery,
    fetch,
)

# The v1 enums, rejected at the edge rather than answered with an empty page. A
# caller asking for SMPDB has made a mistake, and an empty result would hide it.
RESOURCES = ('KEGG', 'Reactome', 'WikiPathways', 'MACdb', 'ClassyFire')
SET_TYPES = ('disease', 'pathway', 'chemical_class')

# Everything the route accepts. Litestar ignores a query parameter it does not
# know, which would answer the wrong question with a 200: `?hmdb=HMDB00077`
# reads as "the memberships of this metabolite" and would return the
# unfiltered first page. The identifier columns are published but not
# filterable in v1, so naming one has to fail loudly.
QUERY_PARAMS = (
    'resource',
    'set_type',
    'organism',
    'set_source_id',
    'metabolite_entity_id',
    'limit',
    'offset',
)


def _db_url(request: Request) -> str:
    url = (
        request.app.state.get('omnipath_db_url')
        or os.environ.get('OMNIPATH_DB_URL')
        or os.environ.get('DATABASE_URL')
    )
    if not url:
        raise HTTPException(
            status_code=503,
            detail='MetSigDB database not configured (set OMNIPATH_DB_URL).',
        )
    return url


def _connect(request: Request):
    try:
        import psycopg2
    except ImportError as exc:  # pragma: no cover
        raise HTTPException(status_code=503, detail='psycopg2 not installed') from exc
    try:
        return psycopg2.connect(_db_url(request))
    except Exception as exc:  # pragma: no cover - connection failure
        raise HTTPException(
            status_code=503, detail=f'MetSigDB database unavailable: {exc}'
        ) from exc


def _reject_unknown_parameters(request: Request) -> None:
    """Refuse a request naming anything the contract does not publish."""
    unknown = sorted(set(request.query_params) - set(QUERY_PARAMS))
    if unknown:
        raise HTTPException(
            status_code=400,
            detail=(
                f'Unsupported query parameter: {", ".join(unknown)}. '
                f'Supported: {", ".join(QUERY_PARAMS)}.'
            ),
        )


def _checked(values: list[str] | None, allowed: tuple[str, ...], name: str):
    """One multi-valued filter, rejected when it names something v1 does not."""
    if not values:
        return ()
    unknown = [value for value in values if value not in allowed]
    if unknown:
        raise HTTPException(
            status_code=400,
            detail=(
                f'Unsupported {name}: {", ".join(sorted(unknown))}. '
                f'Supported: {", ".join(allowed)}.'
            ),
        )
    return tuple(values)


class MetSigDBController(Controller):
    """Read-only membership rows from the built MetSigDB substrate."""

    path = '/sets/metsigdb'

    # The handler blocks: psycopg2 is synchronous, and a full page over a
    # 3.5-million-row table is not instant. Run it off the event loop so one
    # slow query does not stall every other request.
    @get('/', sync_to_thread=True)
    def memberships(
        self,
        request: Request,
        resource: list[str] | None = Parameter(default=None),
        set_type: list[str] | None = Parameter(default=None),
        organism: int | None = Parameter(default=None),
        set_source_id: str | None = Parameter(default=None),
        metabolite_entity_id: str | None = Parameter(default=None),
        limit: int = Parameter(default=DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
        offset: int = Parameter(default=0, ge=0),
    ) -> dict[str, Any]:
        """Membership rows in the published default projection.

        Every response is paged. No filter combination returns the whole
        substrate, and the row order is stable, so paging through a result set
        neither repeats nor skips a membership.

        An ``organism`` filter matches explicit values only. Null-organism rows
        stay in the dataset and out of an organism-filtered response.
        """
        _reject_unknown_parameters(request)

        spec = MetSigDBQuery(
            resource=_checked(resource, RESOURCES, 'resource'),
            set_type=_checked(set_type, SET_TYPES, 'set_type'),
            organism=organism,
            set_source_id=set_source_id,
            metabolite_entity_id=metabolite_entity_id,
            limit=limit,
            offset=offset,
        )

        conn = _connect(request)
        try:
            try:
                rows = fetch(conn, spec)
            except Exception as exc:
                raise HTTPException(
                    status_code=503,
                    detail=(
                        'metsigdb_membership is not present or not readable; '
                        f'run the MetSigDB build step. ({exc})'
                    ),
                ) from exc
        finally:
            conn.close()

        projected = project_rows(rows)
        return {'count': len(projected), 'rows': projected}
