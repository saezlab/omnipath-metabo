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
Uniform specialized-network API (Milestone G).

Reads the build-side ``network_registry`` and each network's combined-contract
matview from the MAIN OmniPath Postgres (``OMNIPATH_DB_URL`` / ``DATABASE_URL``).
Generic over networks — a new network onboarded by the build-side framework
(a definition alone) is served here with no API change. Read-only; no heavy
compute on the request path. An absent view fails clearly (503).
"""

from __future__ import annotations

__all__ = ['NetworksController']

import io
import os
from typing import Any

from litestar import Controller, Request, get
from litestar.exceptions import HTTPException, NotFoundException
from litestar.params import Parameter
from litestar.response import Response

_REGISTRY_SCHEMA = 'public'
_DEFAULT_LIMIT = 1000
_MAX_LIMIT = 100_000


def _db_url(request: Request) -> str:
    url = (
        request.app.state.get('omnipath_db_url')
        or os.environ.get('OMNIPATH_DB_URL')
        or os.environ.get('DATABASE_URL')
    )
    if not url:
        raise HTTPException(
            status_code=503,
            detail='Network database not configured (set OMNIPATH_DB_URL).',
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
            status_code=503, detail=f'Network database unavailable: {exc}'
        ) from exc


# Temporary compatibility guard, and meant to be deleted whole.
#
# A registry row that names no combined relation is a preset: a dataset served
# by filtering the main service's interaction fact table, with no matview of its
# own anywhere. These routes resolve everything through the schema/relation pair,
# so without this guard a preset would probe the relation `"None"."None"` — the
# status route would then report a live, queryable dataset as absent, which is
# the worst of the available failures because it is silent. The routes below say
# what a preset is instead, and send the caller to the main service for its rows.
# This whole shim goes when the bespoke matviews retire and these routes with
# them; nothing here should grow into serving presets from this service.
def _is_preset(entry: dict[str, Any]) -> bool:
    """Whether the registry row describes a preset rather than a matview."""
    return entry.get('combined_relation') is None


def _preset_elsewhere(name: str) -> HTTPException:
    """The 'not served here' answer, naming the service that does serve it."""
    return HTTPException(
        status_code=501,
        detail=(
            f'{name} is a preset over the interaction fact table and is not '
            f'served by this service; request it from the main OmniPath API '
            f'at /interactions/{name}.'
        ),
    )


def _registry_row(conn, name: str) -> dict[str, Any]:
    from psycopg2.extras import RealDictCursor

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(
            f'SELECT name, kind, schema_name, combined_relation, included_sources, '
            f'built_at FROM {_REGISTRY_SCHEMA}.network_registry WHERE name = %s',
            [name],
        )
        row = cur.fetchone()
    if row is None:
        raise NotFoundException(detail=f'Unknown network: {name}')
    return dict(row)


class NetworksController(Controller):
    """Uniform read API over the registered specialized networks."""

    path = '/networks'

    @get('/')
    def list_networks(self, request: Request) -> list[dict[str, Any]]:
        """All registered networks (from ``network_registry``)."""
        from psycopg2.extras import RealDictCursor

        conn = _connect(request)
        try:
            with conn.cursor(cursor_factory=RealDictCursor) as cur:
                try:
                    cur.execute(
                        f'SELECT name, kind, schema_name, combined_relation, '
                        f'included_sources, built_at '
                        f'FROM {_REGISTRY_SCHEMA}.network_registry ORDER BY name'
                    )
                except Exception as exc:
                    raise HTTPException(
                        status_code=503,
                        detail='network_registry not present; run the build.',
                    ) from exc
                return [dict(row) for row in cur.fetchall()]
        finally:
            conn.close()

    @get('/{name:str}/status')
    def network_status(self, request: Request, name: str) -> dict[str, Any]:
        """View presence + row count + the build id the data was built against."""
        conn = _connect(request)
        try:
            entry = _registry_row(conn, name)
            schema, relation = entry['schema_name'], entry['combined_relation']
            preset = _is_preset(entry)
            with conn.cursor() as cur:
                present, row_count = True, None
                if not preset:
                    cur.execute(
                        'SELECT to_regclass(%s)', [f'{schema}.{relation}']
                    )
                    present = cur.fetchone()[0] is not None
                    if present:
                        cur.execute(
                            f'SELECT count(*) FROM "{schema}"."{relation}"'
                        )
                        row_count = int(cur.fetchone()[0])
                cur.execute('SELECT build_id FROM public.build_manifest')
                manifest = cur.fetchone()
            return {
                'name': name,
                'kind': 'preset' if preset else entry['kind'],
                'present': present,
                'row_count': row_count,
                'build_id': manifest[0] if manifest else None,
            }
        finally:
            conn.close()

    @get('/{name:str}/resources')
    def network_resources(self, request: Request, name: str) -> dict[str, Any]:
        """The per-source views/sources composing the network."""
        conn = _connect(request)
        try:
            entry = _registry_row(conn, name)
            return {
                'name': name,
                'kind': entry['kind'],
                'included_sources': entry['included_sources'],
                'combined_relation': entry['combined_relation'],
            }
        finally:
            conn.close()

    @get('/{name:str}/interactions')
    def network_interactions(
        self,
        request: Request,
        name: str,
        source: str | None = Parameter(default=None),
        limit: int = Parameter(default=_DEFAULT_LIMIT, ge=1, le=_MAX_LIMIT),
        offset: int = Parameter(default=0, ge=0),
        format: str = Parameter(default='json'),
    ) -> Any:
        """Rows from the network's combined contract (JSON or streamed parquet)."""
        import pandas as pd

        conn = _connect(request)
        try:
            entry = _registry_row(conn, name)
            if _is_preset(entry):
                raise _preset_elsewhere(name)
            schema, relation = entry['schema_name'], entry['combined_relation']
            where, params = '', {'limit': limit, 'offset': offset}
            if source:
                # Generic across text / text[] `sources` columns: substring on text.
                where = (
                    'WHERE EXISTS (SELECT 1 FROM information_schema.columns c '
                    "WHERE c.table_schema = %(schema)s AND c.table_name = "
                    "%(relation)s AND c.column_name = 'sources') "
                    'AND sources::text ILIKE %(source)s'
                )
                params.update(
                    {'schema': schema, 'relation': relation, 'source': f'%{source}%'}
                )
            query = (
                f'SELECT * FROM "{schema}"."{relation}" {where} '
                f'ORDER BY 1 LIMIT %(limit)s OFFSET %(offset)s'
            )
            try:
                frame = pd.read_sql_query(query, conn, params=params)
            except Exception as exc:
                raise HTTPException(
                    status_code=503,
                    detail=f'Network view {schema}.{relation} unavailable: {exc}',
                ) from exc
        finally:
            conn.close()

        if format == 'parquet':
            buffer = io.BytesIO()
            frame.to_parquet(buffer, index=False)
            buffer.seek(0)
            return Response(
                content=buffer.getvalue(),
                media_type='application/octet-stream',
                headers={
                    'Content-Disposition': f'attachment; filename="{name}.parquet"'
                },
            )
        return {
            'name': name,
            'count': int(len(frame)),
            'rows': frame.to_dict(orient='records'),
        }
