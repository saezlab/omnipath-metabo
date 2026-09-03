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
    RESOURCES,
    SET_SUB_TYPES,
    SET_TYPES,
    count,
    count_sets,
    fetch_groups,
    fetch_page,
    identifier_columns,
)

# The v1 enums live with the query layer, which canonicalises against them.

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
    'set',
    'entity',
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


def _listed(values: list[str] | None) -> tuple[str, ...]:
    """A repeated or comma-separated parameter as one tuple."""
    if not values:
        return ()
    return tuple(
        part.strip()
        for value in values
        for part in value.split(',')
        if part.strip()
    )


def _entity(values: list[str] | None) -> tuple[str, ...]:
    """The `entity` filter, with every value's namespace recognised at the edge.

    Recognition happens here so an unrecognisable value is a 400 naming it,
    rather than an empty page that reads as "no such metabolite".
    """
    wanted = _listed(values)
    if not wanted:
        return ()
    try:
        for value in wanted:
            identifier_columns(value)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return wanted


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
    """One multi-valued filter, rejected when it names something v1 does not.

    Matched without regard to case since cycle 012, and the canonical spelling
    is what reaches the query. An invalid value is still refused whatever its
    case: `wikipathways` is `WikiPathways`, but `smpdb` is nothing.
    """
    values = _listed(values)
    if not values:
        return ()
    canonical = {value.lower(): value for value in allowed}
    unknown = [value for value in values if value.lower() not in canonical]
    if unknown:
        raise HTTPException(
            status_code=400,
            detail=(
                f'Unsupported {name}: {", ".join(sorted(unknown))}. '
                f'Supported: {", ".join(allowed)}.'
            ),
        )
    return tuple(canonical[value.lower()] for value in values)


# The operation description the OpenAPI document publishes. A consumer who gets
# an empty result from an identifier query has to be able to find out why
# without reading the substrate, so the coverage table is part of the contract
# rather than a note in the specification.
MEMBERSHIPS_DESCRIPTION = """Membership rows in the published default projection.

Every response is paged. No filter combination returns the whole
substrate, and the row order is stable, so paging through a result set
neither repeats nor skips a membership.

An ``organism`` filter matches explicit values only. Null-organism rows
stay in the dataset and out of an organism-filtered response.

``has_more`` says whether another page follows, and costs one extra
row. ``total`` is the size of the whole result set and is computed only
when the request asks for it, because counting a filter that matches
three million rows is work nobody should pay for by default.

**Fields.** Thirteen by default. ``fields=`` adds any of the other nine
by name, and ``fields=all`` returns all twenty-two. An unknown name is
refused rather than ignored.

**Grouping.** ``group=true`` groups by resource and then by set. Under
grouping ``limit``, ``offset``, ``has_more`` and ``total`` count
**sets**, not rows, so no set straddles a page. A set is never
truncated: sets are added whole until the next would cross a 50,000-row
ceiling, and a set larger than the ceiling is returned complete and
alone.

**Identifiers.** ``entity`` accepts an internal entity id, an InChIKey,
or an HMDB, ChEBI, KEGG or PubChem identifier, recognised from the value
itself. A bare integer matches ChEBI **or** PubChem and returns the
union, because both namespaces are bare integers. SMILES is published
and reachable through ``fields``, but is not searchable: it has no
distinguishing shape.

Nothing is translated on the request path. Values are matched against
the identifier columns the substrate already publishes, and the
substrate publishes only what each source supplied. **Coverage is
therefore resource-dependent, and an empty result is often correct:**

===============  ===========  ======  =======  ======  ========
resource         metabolites  hmdb    chebi    kegg    pubchem
===============  ===========  ======  =======  ======  ========
ClassyFire       145,937      100%    17%      4%      73%
MACdb            5,389        72%     83%      36%     100%
WikiPathways     2,789        41%     74%      39%     73%
Reactome         2,191        0%      76%      3%      1%
KEGG             1,799        20%     77%      100%    98%
===============  ===========  ======  =======  ======  ========

A query by HMDB identifier cannot return a Reactome row, because no
Reactome metabolite carries one. Measured against build
``9eb5a917e9ed``; the figures move when the upstream resolution work
lands.

**Deprecated.** ``set_source_id`` and ``metabolite_entity_id`` are the
cycle 010 names for ``set`` and ``entity``. Both still work.
"""


class MetSigDBController(Controller):
    """Read-only membership rows from the built MetSigDB substrate."""

    path = '/sets/metsigdb'

    # The handler blocks: psycopg2 is synchronous, and a full page over a
    # 3.5-million-row table is not instant. Run it off the event loop so one
    # slow query does not stall every other request.
    @get('/', description=MEMBERSHIPS_DESCRIPTION, sync_to_thread=True)
    def memberships(
        self,
        request: Request,
        resource: list[str] | None = Parameter(default=None),
        set_type: list[str] | None = Parameter(default=None),
        set_sub_type: list[str] | None = Parameter(default=None),
        organism: int | None = Parameter(default=None),
        set: list[str] | None = Parameter(default=None),
        entity: list[str] | None = Parameter(default=None),
        # The cycle 010 names. Marked deprecated so a client generator and the
        # schema page both show the rename before anything enforces it; a
        # consumer should see it coming rather than meet a 400 one day.
        set_source_id: str | None = Parameter(
            default=None,
            description='Deprecated. Use `set`.',
            schema_extra={'deprecated': True},
        ),
        metabolite_entity_id: str | None = Parameter(
            default=None,
            description='Deprecated. Use `entity`.',
            schema_extra={'deprecated': True},
        ),
        fields: list[str] | None = Parameter(default=None),
        group: bool = Parameter(default=False),
        limit: int = Parameter(default=DEFAULT_LIMIT, ge=1, le=MAX_LIMIT),
        offset: int = Parameter(default=0, ge=0),
        total: bool = Parameter(default=False),
    ) -> MetSigDBPage | MetSigDBGroupedPage:
        """Membership rows in the published default projection.

        The full contract — fields, grouping, identifier coverage and the
        deprecated names — is `MEMBERSHIPS_DESCRIPTION`, which is also what the
        OpenAPI document publishes. It lives outside the docstring because this
        service does not set `use_handler_docstrings`, and turning that on would
        publish every other route's docstring too.
        """
        _reject_unknown_parameters(request)
        projection = _fields(fields)

        spec = MetSigDBQuery(
            set=_listed(set),
            entity=_entity(entity),
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
