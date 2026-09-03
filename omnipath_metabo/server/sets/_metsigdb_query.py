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
    'count_query',
    'fetch',
    'fetch_page',
    'count',
    'ROW_CEILING',
    'IDENTIFIER_COLUMNS',
    'RESOURCES',
    'SET_SUB_TYPES',
    'SET_TYPES',
    'canonical',
    'count_sets',
    'fetch_groups',
    'identifier_columns',
]

import re
from dataclasses import dataclass
from typing import Any

from omnipath_metabo.server.sets._metsigdb_projection import (
    COLUMN_NAMES,
    DEFAULT_FIELDS,
    MEMBER_FIELDS,
    SET_FIELDS,
    columns_for,
)

TABLE = 'public.metsigdb_membership'

# The same page bounds the network routes of this service already serve.
DEFAULT_LIMIT = 1000
MAX_LIMIT = 100_000

# What a grouped response may carry in total. `limit` counts sets, and set
# sizes in this substrate span four orders of magnitude — from one member to
# 145,937 — so a page of N sets has no bounded row count without this.
#
# Sets are added whole until the next would cross it. A set is never cut, so a
# set larger than the ceiling is returned complete and alone: the ceiling bounds
# how many sets share a page, never whether a published set can be returned.
ROW_CEILING = 50_000

# The v1 enums. They live here rather than in the route because the query layer
# canonicalises against them: a filter value is matched to its published
# spelling before it reaches SQL, which is what makes case-insensitivity free.
RESOURCES = ('KEGG', 'Reactome', 'WikiPathways', 'MACdb', 'ClassyFire')
SET_TYPES = ('disease', 'pathway', 'chemical_class')
SET_SUB_TYPES = (
    'cancer',
    'phenotype',
    'medical intervention',
    'gene abnormality',
    'genotype',
    'overview_map',
    'metabolic_map',
)

VOCABULARIES: dict[str, tuple[str, ...]] = {
    'resource': RESOURCES,
    'set_type': SET_TYPES,
    'set_sub_type': SET_SUB_TYPES,
}


def canonical(field: str, value: str) -> str:
    """A vocabulary value in its published spelling, whatever case it arrived in.

    Case-insensitivity is a normalisation at the edge, **not** a `lower()` in
    the predicate. Lowering the column costs the primary key: `resource` leads
    it, so `resource = ANY(...)` becomes an index condition while
    `lower(resource) = ANY(...)` becomes a filter. Measured on this substrate,
    the difference for one five-row page of WikiPathways is 6 buffers against
    91,442, because the scan discards 3,503,370 ClassyFire rows first.
    """
    for published in VOCABULARIES.get(field, ()):
        if published.lower() == value.lower():
            return published
    return value

# The identifier columns `entity` accepts, and the shape that names each. A
# caller pastes the identifier they hold; the route works out the namespace. A
# type parameter would make the easy case harder and the hard case no easier,
# because the namespace is often what the caller is asking about.
#
# `smiles` is absent on purpose. It has no distinguishing shape — any string
# could be a SMILES — so a rule for it would swallow every value the others
# declined. It is also the one identifier column left unindexed.
IDENTIFIER_COLUMNS: tuple[str, ...] = (
    'metabolite_entity_id',
    'inchikey',
    'hmdb',
    'chebi',
    'kegg',
    'pubchem',
)

_UUID = re.compile(r'^[0-9a-f]{8}(-[0-9a-f]{4}){3}-[0-9a-f]{12}$', re.I)
_INCHIKEY = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
_HMDB = re.compile(r'^HMDB\d{5,11}$', re.I)
_CHEBI_PREFIXED = re.compile(r'^CHEBI:\d+$', re.I)
_KEGG = re.compile(r'^C\d{5}$', re.I)
_BARE_INTEGER = re.compile(r'^\d+$')


def identifier_columns(value: str) -> tuple[str, ...]:
    """Which published columns a value could be an identifier in.

    Usually one. A bare integer is two: ChEBI and PubChem identifiers are both
    bare integers, and the value alone cannot decide between them. The union is
    the honest answer — guessing one namespace would silently drop the other's
    memberships from a result the caller believes is complete.

    Raises for a value no rule claims, rather than returning nothing and
    answering with an empty page that reads as "no such metabolite".
    """
    value = value.strip()

    if _UUID.match(value):
        return ('metabolite_entity_id',)
    if _INCHIKEY.match(value):
        return ('inchikey',)
    if _HMDB.match(value):
        return ('hmdb',)
    if _CHEBI_PREFIXED.match(value):
        return ('chebi',)
    if _KEGG.match(value):
        return ('kegg',)
    if _BARE_INTEGER.match(value):
        return ('chebi', 'pubchem')

    raise ValueError(
        f'Unrecognised entity identifier: {value}. Supported: an internal '
        'entity id, an InChIKey, or an HMDB, ChEBI, KEGG or PubChem '
        'identifier. SMILES is published but not searchable.'
    )


def _entity_clause(
    values: tuple[str, ...],
    params: dict[str, Any],
) -> str:
    """One OR-group matching any supplied identifier in any column it could be.

    `CHEBI:` is stripped before binding: the substrate stores bare numbers, and
    a caller who pastes the prefixed form means the same identifier.
    """
    by_column: dict[str, list[str]] = {}
    for value in values:
        value = value.strip()
        for column in identifier_columns(value):
            stored = value.split(':', 1)[1] if column == 'chebi' and ':' in value else value
            by_column.setdefault(column, []).append(stored)

    predicates = []
    for index, (column, wanted) in enumerate(sorted(by_column.items())):
        key = f'entity_{index}'
        # The entity id is a uuid column; the others are text. Without the cast
        # Postgres has no `uuid = text` operator and the query fails rather than
        # falling back to a slower comparison.
        cast = '::uuid[]' if column == 'metabolite_entity_id' else '::text[]'
        predicates.append(f'{column} = ANY(%({key})s{cast})')
        params[key] = wanted

    return f'({" OR ".join(predicates)})' 

# The columns a query reads come from the projection's field registry, which is
# the single source of truth for what the substrate publishes. This module used
# to carry its own copy of the twenty-two column names; two lists of the same
# thing drift, and the registry is the one the response shape is built from.

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
    set_sub_type: tuple[str, ...] = ()
    organism: int | None = None
    # Cycle 010 named these after the columns that store them. Cycle 012 names
    # them after what they select, and the old names are gone rather than
    # deprecated: the route refuses an unsupported parameter, so a caller using
    # one is told, instead of being answered as though they had filtered.
    set: tuple[str, ...] = ()
    entity: tuple[str, ...] = ()
    # The response fields this request asked for, which decide the columns the
    # query reads. A page of thirteen fields does not pull twenty-two off the
    # disk and discard nine — `set_context` and `provenance_record` alone are
    # the two heaviest columns in the table.
    fields: tuple[str, ...] = DEFAULT_FIELDS
    limit: int = DEFAULT_LIMIT
    offset: int = 0

    def columns(self) -> tuple[str, ...]:
        """The substrate columns this request needs, in contract order."""
        return columns_for(self.fields)


def build_query(spec: MetSigDBQuery) -> tuple[str, dict[str, Any]]:
    """The SQL and bind parameters for one request.

    Separate from execution so the shape of a filter set can be asserted
    without a database.
    """
    # ANY over an array keeps one bind parameter whatever the caller asks for,
    # so the plan does not change with the value count. `organism` is equality,
    # so a null-organism row never matches: it stays in the substrate and out of
    # an organism-filtered response, which is what the contract says.
    where, params = _predicates(spec)
    params |= {'limit': spec.limit, 'offset': spec.offset}
    clause = f'WHERE {" AND ".join(where)} ' if where else ''
    columns = ', '.join(spec.columns())

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


def count_query(spec: MetSigDBQuery) -> tuple[str, dict[str, Any]]:
    """How many rows the filters match, ignoring the page.

    Separate from `build_query` because a caller pays for this only when it
    asks: counting a filter that matches three million rows is real work.
    """
    sql, params = build_query(spec)
    where = sql.partition(f'FROM {TABLE} ')[2].partition('ORDER BY')[0]
    return (
        f'SELECT count(*) FROM {TABLE} {where}',
        {k: v for k, v in params.items() if k not in ('limit', 'offset')},
    )


def fetch_page(conn, spec: MetSigDBQuery) -> tuple[list[dict[str, Any]], bool]:
    """One page, and whether another one follows.

    Asks for one row beyond the page and drops it. That answers "is there
    more" for the cost of a single row, where a count would re-run the filter
    over the whole substrate.
    """
    probe = MetSigDBQuery(
        **{
            **{f: getattr(spec, f) for f in spec.__dataclass_fields__},
            'limit': spec.limit + 1,
        }
    )
    rows = fetch(conn, probe)
    has_more = len(rows) > spec.limit
    return rows[: spec.limit], has_more


def count(conn, spec: MetSigDBQuery) -> int:
    """The size of the whole result set."""
    sql, params = count_query(spec)
    with conn.cursor() as cur:
        cur.execute(sql, params)
        return int(cur.fetchone()[0])


# ------------------------------------------------------------ grouped results


def _predicates(spec: MetSigDBQuery) -> tuple[list[str], dict[str, Any]]:
    """Every filter this request carries, as clauses and bind parameters.

    One builder for all four query shapes. A member-level filter has to narrow
    the flat page, the set page, the member fetch and both counts identically,
    or `returned` disagrees with what a grouped response carries.
    """
    where: list[str] = []
    params: dict[str, Any] = {}

    # The three closed vocabularies match without regard to case, by canonical
    # spelling rather than by lowering the column. See `canonical`.
    for field in ('resource', 'set_type', 'set_sub_type'):
        values = getattr(spec, field)
        if values:
            where.append(f'{field} = ANY(%({field})s)')
            params[field] = [canonical(field, value) for value in values]

    # `set` never reads as an entity identifier: MACdb set ids are bare
    # integers and collide with ChEBI ids. Cycle 010's first set-name
    # measurement reported 645 MACdb sets instead of 269 for that reason.
    if spec.set:
        where.append('set_source_id = ANY(%(set)s)')
        params['set'] = list(spec.set)

    if spec.entity:
        where.append(_entity_clause(spec.entity, params))

    if spec.organism is not None:
        where.append('organism = %(organism)s')
        params['organism'] = spec.organism

    return where, params


def _where(spec: MetSigDBQuery) -> tuple[str, dict[str, Any]]:
    """The filter clause and its parameters, shared by every query shape."""
    where, params = _predicates(spec)
    return (f'WHERE {" AND ".join(where)} ' if where else ''), params


def _set_columns(spec: MetSigDBQuery) -> tuple[str, ...]:
    """The set-side columns this request needs.

    `resource` and `set_source_id` are always read: they are the group key and
    the paging order, whether or not the caller asked to see them.
    """
    wanted = [field for field in spec.fields if field in SET_FIELDS]
    columns = [COLUMN_NAMES.get(field, field) for field in wanted]
    for key in ('resource', 'set_source_id'):
        if key not in columns:
            columns.append(key)
    return tuple(columns)


def set_page_query(spec: MetSigDBQuery) -> tuple[str, dict[str, Any]]:
    """The page's sets, with how many members each has under this filter.

    One row per set rather than per membership, so the page can be measured in
    sets before a single member is read. `count(*)` is the filtered count; the
    set's published population is `set_size`, which the same row carries.

    One row beyond the page is asked for, which is how `has_more` costs a row
    rather than a second count.
    """
    clause, params = _where(spec)
    columns = _set_columns(spec)
    grouped = ', '.join(columns)
    params |= {'limit': spec.limit + 1, 'offset': spec.offset}

    return (
        f'SELECT {grouped}, count(*) AS returned FROM {TABLE} {clause}'
        f'GROUP BY {grouped} ORDER BY resource, set_source_id '
        f'LIMIT %(limit)s OFFSET %(offset)s'
    ), params


def count_sets_query(spec: MetSigDBQuery) -> tuple[str, dict[str, Any]]:
    """How many sets the filter matches, for a grouped `total`."""
    clause, params = _where(spec)
    return (
        f'SELECT count(*) FROM (SELECT 1 FROM {TABLE} {clause}'
        f'GROUP BY resource, set_source_id) AS matched_sets'
    ), params


def members_query(
    spec: MetSigDBQuery,
    keys: list[tuple[str, str]],
) -> tuple[str, dict[str, Any]]:
    """The members of a named list of sets, under the same filter.

    The filter is applied again rather than trusted from the set query: a
    member-level filter has to narrow the members too, or `returned` would
    disagree with what the response carries.
    """
    clause, params = _where(spec)
    member_columns = [
        COLUMN_NAMES.get(field, field)
        for field in spec.fields
        if field in MEMBER_FIELDS
    ]
    columns = ['resource', 'set_source_id', *member_columns]
    joiner = 'AND' if clause else 'WHERE'
    # Two parallel arrays rather than a list of pairs: psycopg2 adapts a list of
    # lists to a Postgres array, which cannot be coerced to a row list. Unnesting
    # them together rebuilds the pairs in SQL.
    params['key_resources'] = [key[0] for key in keys]
    params['key_sets'] = [key[1] for key in keys]

    return (
        f'SELECT {", ".join(columns)} FROM {TABLE} {clause}'
        f'{joiner} (resource, set_source_id) IN ('
        f'SELECT r, s FROM unnest(%(key_resources)s::text[], %(key_sets)s::text[]) '
        f'AS pairs(r, s)) ORDER BY {ORDER_BY}'
    ), params


def _fit_the_ceiling(rows: list[dict[str, Any]], limit: int) -> tuple[list, bool]:
    """As many whole sets as the page and the ceiling allow.

    Two bounds, and the tighter one wins. `limit` counts sets. `ROW_CEILING`
    counts members, and a set is added whole or not at all — except the first,
    which is returned even when it exceeds the ceiling alone, because a
    published set must stay reachable.
    """
    has_more = len(rows) > limit
    rows = rows[:limit]

    taken: list[dict[str, Any]] = []
    total = 0
    for row in rows:
        if taken and total + row['returned'] > ROW_CEILING:
            has_more = True
            break
        taken.append(row)
        total += row['returned']

    return taken, has_more


def fetch_groups(conn, spec: MetSigDBQuery) -> tuple[list[dict[str, Any]], bool]:
    """One page of sets, grouped by resource, each with its members.

    Two queries rather than one: the sets first, so the page can be measured and
    the ceiling applied before any member is read, then the members of the sets
    that survived. Reading members first would mean fetching rows only to
    discard them, and for ClassyFire that is millions of them.
    """
    from psycopg2.extras import RealDictCursor

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        sql, params = set_page_query(spec)
        cur.execute(sql, params)
        set_rows = [dict(row) for row in cur.fetchall()]

    set_rows, has_more = _fit_the_ceiling(set_rows, spec.limit)
    if not set_rows:
        return [], has_more

    keys = [(row['resource'], row['set_source_id']) for row in set_rows]
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        sql, params = members_query(spec, keys)
        cur.execute(sql, params)
        member_rows = [dict(row) for row in cur.fetchall()]

    return _assemble(spec, set_rows, member_rows), has_more


def _assemble(
    spec: MetSigDBQuery,
    set_rows: list[dict[str, Any]],
    member_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Sets and members into resource groups, keeping the query's order."""
    from omnipath_metabo.server.sets._metsigdb_projection import (
        project_row,
        project_rows,
    )

    set_fields = tuple(f for f in spec.fields if f in SET_FIELDS)
    member_fields = tuple(f for f in spec.fields if f in MEMBER_FIELDS)

    by_key: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in member_rows:
        by_key.setdefault((row['resource'], row['set_source_id']), []).append(row)

    groups: list[dict[str, Any]] = []
    for row in set_rows:
        key = (row['resource'], row['set_source_id'])
        if not groups or groups[-1]['resource'] != row['resource']:
            groups.append({'resource': row['resource'], 'sets': []})

        one_set = project_row(row, set_fields)
        one_set.pop('resource', None)
        one_set['returned'] = row['returned']
        one_set['members'] = project_rows(by_key.get(key, ()), member_fields)
        groups[-1]['sets'].append(one_set)

    return groups


def count_sets(conn, spec: MetSigDBQuery) -> int:
    """How many sets the filter matches."""
    sql, params = count_sets_query(spec)
    with conn.cursor() as cur:
        cur.execute(sql, params)
        return int(cur.fetchone()[0])
