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
Structure consistency diagnostics (spec 011 contracts/quality-control-api.md,
T104-T105).

Read-only over the precomputed ``structure_consistency_finding`` /
``structure_consistency_summary`` tables built by
:mod:`omnipath_metabo.postbuild._qc_layer` -- no structure computation ever
runs on the request path (Principle II).
"""

from __future__ import annotations

__all__ = ['QualityControlController']

import os
from typing import Any

from litestar import Controller, Request, get
from litestar.exceptions import HTTPException
from litestar.params import Parameter

_DEFAULT_LIMIT = 1000
_MAX_LIMIT = 100_000

_UNAVAILABLE_META = {
    'available': False,
    'capability': 'structure_substrate',
    'reason': (
        'the build ran without the chemistry toolkit, so no structures '
        'were compared'
    ),
}


def _db_url(request: Request) -> str:
    url = (
        request.app.state.get('omnipath_db_url')
        or os.environ.get('OMNIPATH_DB_URL')
        or os.environ.get('DATABASE_URL')
    )
    if not url:
        raise HTTPException(
            status_code=503,
            detail='OmniPath database not configured (set OMNIPATH_DB_URL).',
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
            status_code=503, detail=f'OmniPath database unavailable: {exc}'
        ) from exc


def _schema_present(cur) -> bool:
    """Whether the metabo QC tables exist at all -- absent exactly when the
    build ran without the chemistry toolkit (``ensure_rdkit_extension``
    raises before either table is created, T098).
    """
    cur.execute("SELECT to_regclass('public.structure_consistency_summary')")
    return cur.fetchone()[0] is not None


def _source_id(cur, name: str) -> int | None:
    # Called with a RealDictCursor in every caller -- row['source_id'], not
    # row[0] (RealDictRow has no integer indexing).
    cur.execute('SELECT source_id FROM public.data_source WHERE name = %s', [name])
    row = cur.fetchone()
    return int(row['source_id']) if row else None


def _checked_response(
    conn,
    check_kind: str,
    *,
    source: str | None,
    authority: str | None,
    verdict: str | None,
    limit: int,
    offset: int,
) -> dict[str, Any]:
    from psycopg2.extras import RealDictCursor

    with conn.cursor() as cur:
        if not _schema_present(cur):
            return {'summary': [], 'findings': [], 'meta': dict(_UNAVAILABLE_META)}

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        source_id = _source_id(cur, source) if source else None
        if source and source_id is None:
            raise HTTPException(status_code=404, detail=f'Unknown source: {source}')
        authority_id = _source_id(cur, authority) if authority else None
        if authority and authority_id is None:
            raise HTTPException(
                status_code=404, detail=f'Unknown authority: {authority}'
            )

        where = ['check_kind = %(check_kind)s']
        params: dict[str, Any] = {'check_kind': check_kind}
        if source_id is not None:
            where.append('summ.source_id = %(source_id)s')
            params['source_id'] = source_id
        if authority_id is not None:
            where.append('summ.authority_source_id = %(authority_id)s')
            params['authority_id'] = authority_id
        if verdict is not None:
            where.append('summ.verdict = %(verdict)s')
            params['verdict'] = verdict
        where_sql = ' AND '.join(where)

        cur.execute(
            f"""
            SELECT ds.name AS source, au.name AS authority, summ.verdict,
                   summ.finding_count
            FROM public.structure_consistency_summary summ
            JOIN public.data_source ds ON ds.source_id = summ.source_id
            LEFT JOIN public.data_source au
              ON au.source_id = summ.authority_source_id
            WHERE {where_sql}
            """,
            params,
        )
        by_pair: dict[tuple[str, str | None], dict[str, Any]] = {}
        for row in cur.fetchall():
            key = (row['source'], row['authority'])
            entry = by_pair.setdefault(
                key,
                {
                    'source': row['source'],
                    'authority': row['authority'],
                    'compared': 0,
                    'agree': 0,
                    'layer_difference': 0,
                    'different_structure': 0,
                    'unparsable': 0,
                },
            )
            entry[row['verdict']] = int(row['finding_count'])
            entry['compared'] += int(row['finding_count'])
        summary = sorted(
            by_pair.values(), key=lambda e: (e['source'], e['authority'] or '')
        )

        find_where = ['fnd.check_kind = %(check_kind)s']
        find_params: dict[str, Any] = {
            'check_kind': check_kind, 'limit': limit, 'offset': offset,
        }
        if source_id is not None:
            find_where.append('fnd.source_id = %(source_id)s')
            find_params['source_id'] = source_id
        if authority_id is not None:
            find_where.append('fnd.authority_source_id = %(authority_id)s')
            find_params['authority_id'] = authority_id
        if verdict is not None:
            find_where.append('fnd.verdict = %(verdict)s')
            find_params['verdict'] = verdict
        find_where_sql = ' AND '.join(find_where)

        cur.execute(
            f"""
            SELECT ds.name AS source, vit.name AS identifier_type,
                   fnd.value_normalized AS value, au.name AS authority,
                   fnd.structure_a, fnd.structure_b, fnd.verdict, fnd.layer
            FROM public.structure_consistency_finding fnd
            JOIN public.data_source ds ON ds.source_id = fnd.source_id
            JOIN public.vocab_identifier_type vit
              ON vit.identifier_type_id = fnd.identifier_type_id
            LEFT JOIN public.data_source au
              ON au.source_id = fnd.authority_source_id
            WHERE {find_where_sql}
            ORDER BY fnd.finding_id
            LIMIT %(limit)s OFFSET %(offset)s
            """,
            find_params,
        )
        findings = [dict(row) for row in cur.fetchall()]

    return {
        'summary': summary,
        'findings': findings,
        'meta': {'available': True},
    }


class QualityControlController(Controller):
    """Structure consistency diagnostics (three precomputed checks)."""

    path = '/qc/structure'

    @get('/authorities')
    def authorities(self, request: Request) -> dict[str, Any]:
        """Every namespace the build derived, and whether it mints structures."""
        from psycopg2.extras import RealDictCursor

        conn = _connect(request)
        try:
            with conn.cursor() as cur:
                cur.execute("SELECT to_regclass('public.metabo_entity_structure')")
                structure_table_present = cur.fetchone()[0] is not None
            with conn.cursor(cursor_factory=RealDictCursor) as cur:
                if not structure_table_present:
                    # The chemistry toolkit never ran (ensure_rdkit_extension
                    # raises before this table is created) -- still list the
                    # declared authorities, just with an honest zero rather
                    # than joining a relation that does not exist.
                    cur.execute(
                        """
                        SELECT vit.name AS identifier_type, ds.name AS source,
                               ia.is_structure_authority AS structure_authority,
                               0 AS records_with_structure
                        FROM public.identifier_authority ia
                        JOIN public.vocab_identifier_type vit
                          ON vit.identifier_type_id = ia.identifier_type_id
                        JOIN public.data_source ds ON ds.source_id = ia.source_id
                        ORDER BY ds.name
                        """
                    )
                    return {
                        'authorities': [dict(row) for row in cur.fetchall()],
                        'meta': dict(_UNAVAILABLE_META),
                    }
                cur.execute(
                    """
                    SELECT vit.name AS identifier_type, ds.name AS source,
                           ia.is_structure_authority AS structure_authority,
                           count(DISTINCT s.entity_id)
                             FILTER (WHERE ia.is_structure_authority)
                             AS records_with_structure
                    FROM public.identifier_authority ia
                    JOIN public.vocab_identifier_type vit
                      ON vit.identifier_type_id = ia.identifier_type_id
                    JOIN public.data_source ds ON ds.source_id = ia.source_id
                    LEFT JOIN public.entity_evidence_resolution eer
                      ON eer.source_id = ia.source_id
                    LEFT JOIN public.metabo_entity_structure s
                      ON s.entity_id = eer.entity_id
                    GROUP BY vit.name, ds.name, ia.is_structure_authority
                    ORDER BY ds.name
                    """
                )
                rows = [dict(row) for row in cur.fetchall()]
                for row in rows:
                    row['records_with_structure'] = int(
                        row['records_with_structure'] or 0
                    )
            return {'authorities': rows, 'meta': {'available': True}}
        finally:
            conn.close()

    @get('/internal')
    def internal(
        self,
        request: Request,
        source: str | None = Parameter(default=None),
        authority: str | None = Parameter(default=None),
        verdict: str | None = Parameter(default=None),
        limit: int = Parameter(default=_DEFAULT_LIMIT, ge=1, le=_MAX_LIMIT),
        offset: int = Parameter(default=0, ge=0),
    ) -> dict[str, Any]:
        """Do the representations on one record, from one authority, agree?"""
        conn = _connect(request)
        try:
            return _checked_response(
                conn, 'internal', source=source, authority=authority,
                verdict=verdict, limit=limit, offset=offset,
            )
        finally:
            conn.close()

    @get('/cross-reference')
    def cross_reference(
        self,
        request: Request,
        source: str | None = Parameter(default=None),
        authority: str | None = Parameter(default=None),
        verdict: str | None = Parameter(default=None),
        limit: int = Parameter(default=_DEFAULT_LIMIT, ge=1, le=_MAX_LIMIT),
        offset: int = Parameter(default=0, ge=0),
    ) -> dict[str, Any]:
        """Does a record citing an authority's identifier agree with it?"""
        conn = _connect(request)
        try:
            return _checked_response(
                conn, 'cross_reference', source=source, authority=authority,
                verdict=verdict, limit=limit, offset=offset,
            )
        finally:
            conn.close()

    @get('/cross-reference-pairs')
    def cross_reference_pairs(
        self,
        request: Request,
        source: str | None = Parameter(default=None),
        authority: str | None = Parameter(default=None),
        verdict: str | None = Parameter(default=None),
        limit: int = Parameter(default=_DEFAULT_LIMIT, ge=1, le=_MAX_LIMIT),
        offset: int = Parameter(default=0, ge=0),
    ) -> dict[str, Any]:
        """A record with no structure, citing two authorities: do they agree?"""
        conn = _connect(request)
        try:
            return _checked_response(
                conn, 'cross_reference_pair', source=source, authority=authority,
                verdict=verdict, limit=limit, offset=offset,
            )
        finally:
            conn.close()
