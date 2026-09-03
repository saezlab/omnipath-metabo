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

from omnipath_metabo.server.sets._metsigdb_projection import (
    ALL_FIELDS,
    ALL_KEYWORD,
    MetSigDBGroupedPage,
    MetSigDBPage,
    project_rows,
    resolve_fields,
)
from omnipath_metabo.server.sets._metsigdb_query import (
    DEFAULT_LIMIT,
    MAX_LIMIT,
    MetSigDBQuery,
    count,
    count_sets,
    fetch_groups,
    fetch_page,
)

# The v1 enums, rejected at the edge rather than answered with an empty page. A
# caller asking for SMPDB has made a mistake, and an empty result would hide it.
RESOURCES = ('KEGG', 'Reactome', 'WikiPathways', 'MACdb', 'ClassyFire')
SET_TYPES = ('disease', 'pathway', 'chemical_class')

# The finer semantic, where a resource publishes one. MACdb's five trait types
# come from the source; KEGG marks its whole-metabolism overview maps. The
# other three resources leave it null, and filtering on it returns nothing for
# them, which is the honest answer.
SET_SUB_TYPES = (
    'cancer',
    'phenotype',
    'medical intervention',
    'gene abnormality',
    'genotype',
    'overview_map',
    'metabolic_map',
)

# Everything the route accepts. Litestar ignores a query parameter it does not
# know, which would answer the wrong question with a 200: `?hmdb=HMDB00077`
# reads as "the memberships of this metabolite" and would return the
# unfiltered first page. The identifier columns are published but not
# filterable in v1, so naming one has to fail loudly.
QUERY_PARAMS = (
    'resource',
    'set_type',
    'set_sub_type',
    'organism',
    'set_source_id',
    'metabolite_entity_id',
    'fields',
    'group',
    'limit',
    'offset',
    'total',
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


def _fields(requested: list[str] | None):
    """The response projection this request asks for.

    Accepts the parameter repeated and comma-separated, because both spellings
    reach a URL and refusing one is a trap rather than a rule.

    An unknown name becomes a 400 naming it. The projection layer raises the
    error and this turns it into a response: the rule belongs with the field
    registry, and the status code belongs at the edge.
    """
    if not requested:
        return resolve_fields(None)

    names = [
        name.strip()
        for value in requested
        for name in value.split(',')
        if name.strip()
    ]
    try:
        return resolve_fields(names)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


def _grouped(conn, spec, *, total: bool):
    """One page of sets, grouped by resource.

    Every count is in sets rather than rows, which is the whole point of the
    shape: a page of N sets never splits one across a boundary, so a grouped
    response is coherent on its own.
    """
    groups, has_more = fetch_groups(conn, spec)
    return {
        'count': sum(len(group['sets']) for group in groups),
        'has_more': has_more,
        'groups': groups,
        'total': count_sets(conn, spec) if total else None,
    }


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
        set_sub_type: list[str] | None = Parameter(default=None),
        organism: int | None = Parameter(default=None),
        set_source_id: str | None = Parameter(default=None),
        metabolite_entity_id: str | None = Parameter(default=None),
        fields: list[str] | None = Parameter(default=None),
        group: bool = Parameter(default=False),
        limit: int = Parameter(default=DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
        offset: int = Parameter(default=0, ge=0),
        total: bool = Parameter(default=False),
    ) -> MetSigDBPage | MetSigDBGroupedPage:
        """Membership rows in the published default projection.

        Every response is paged. No filter combination returns the whole
        substrate, and the row order is stable, so paging through a result set
        neither repeats nor skips a membership.

        An ``organism`` filter matches explicit values only. Null-organism rows
        stay in the dataset and out of an organism-filtered response.

        ``has_more`` says whether another page follows, and costs one extra
        row. ``total`` is the size of the whole result set and is computed only
        when the request asks for it, because counting a filter that matches
        three million rows is work nobody should pay for by default.
        """
        _reject_unknown_parameters(request)
        projection = _fields(fields)

        spec = MetSigDBQuery(
            resource=_checked(resource, RESOURCES, 'resource'),
            set_type=_checked(set_type, SET_TYPES, 'set_type'),
            set_sub_type=_checked(set_sub_type, SET_SUB_TYPES, 'set_sub_type'),
            organism=organism,
            set_source_id=set_source_id,
            metabolite_entity_id=metabolite_entity_id,
            fields=projection,
            limit=limit,
            offset=offset,
        )

        conn = _connect(request)
        try:
            try:
                if group:
                    return _grouped(conn, spec, total=total)
                rows, has_more = fetch_page(conn, spec)
                matched = count(conn, spec) if total else None
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

        projected = project_rows(rows, projection)
        return {
            'count': len(projected),
            'has_more': has_more,
            'rows': projected,
            'total': matched,
        }
