"""Query-layer tests for /sets/metsigdb (cycle 010).

The query layer owns the SQL: which rows come back for a filter set, in which
order, and how many. It owns nothing about their shape.

Construction is tested without a database. Filter semantics need the built
substrate::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55435/omnipath \
        uv run --with pytest --with psycopg2-binary \
        pytest tests/test_sets_metsigdb_query.py -v
"""

from __future__ import annotations

import os

import pytest

from omnipath_metabo.server.sets._metsigdb_query import (
    DEFAULT_LIMIT,
    MAX_LIMIT,
    MetSigDBQuery,
    build_query,
    fetch,
)

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')


@pytest.fixture(scope='module')
def conn():
    psycopg2 = pytest.importorskip('psycopg2')
    if not DB_URL:
        pytest.skip('no OMNIPATH_DB_URL/DATABASE_URL')
    connection = psycopg2.connect(DB_URL)
    connection.autocommit = True
    try:
        yield connection
    finally:
        connection.close()


# ----------------------------------------------------------------- construction


def test_an_empty_spec_selects_the_whole_table_bounded():
    sql, params = build_query(MetSigDBQuery())
    assert 'WHERE' not in sql
    assert params['limit'] == DEFAULT_LIMIT
    assert params['offset'] == 0


def test_order_is_the_row_identity_so_paging_is_stable():
    sql, _ = build_query(MetSigDBQuery())
    assert 'ORDER BY resource, set_source_id, metabolite_entity_id' in sql


def test_every_query_is_bounded():
    """No filter combination may return the whole substrate in one response."""
    for spec in (MetSigDBQuery(), MetSigDBQuery(resource=('KEGG',))):
        sql, params = build_query(spec)
        assert 'LIMIT' in sql and 'OFFSET' in sql
        assert params['limit'] <= MAX_LIMIT


def test_resource_and_set_type_accept_several_values():
    sql, params = build_query(
        MetSigDBQuery(resource=('KEGG', 'Reactome'), set_type=('pathway',))
    )
    assert 'resource = ANY(%(resource)s)' in sql
    assert 'set_type = ANY(%(set_type)s)' in sql
    assert params['resource'] == ['KEGG', 'Reactome']


def test_scalar_filters_use_equality():
    sql, params = build_query(
        MetSigDBQuery(set_source_id='R-HSA-1059683', organism=9606)
    )
    assert 'set_source_id = %(set_source_id)s' in sql
    assert 'organism = %(organism)s' in sql
    assert params['organism'] == 9606


def test_only_shared_columns_are_filterable():
    """The query layer offers the published filters and no others.

    `set_sub_type` joined them when MACdb's trait type and KEGG's overview maps
    became columns. The six identifier columns are still published and still
    not filterable.
    """
    filterable = {
        field
        for field in MetSigDBQuery.__dataclass_fields__
        if field not in {'limit', 'offset'}
    }
    assert filterable == {
        'resource',
        'set_type',
        'set_sub_type',
        'organism',
        'set_source_id',
        'metabolite_entity_id',
    }


def test_the_sub_type_filter_narrows_rows(conn):
    """MACdb calls every trait a disease; 116 of its 269 are not."""
    rows = fetch(conn, MetSigDBQuery(set_sub_type=('cancer',), limit=50))
    assert rows
    assert {row['resource'] for row in rows} == {'MACdb'}
    assert {row['set_sub_type'] for row in rows} == {'cancer'}


def test_a_page_knows_whether_more_follows(conn):
    from omnipath_metabo.server.sets._metsigdb_query import fetch_page

    rows, has_more = fetch_page(conn, MetSigDBQuery(resource=('KEGG',), limit=10))
    assert len(rows) == 10 and has_more

    rows, has_more = fetch_page(
        conn, MetSigDBQuery(set_source_id='rn00270', limit=1000)
    )
    assert len(rows) == 16 and not has_more


def test_the_count_ignores_the_page(conn):
    from omnipath_metabo.server.sets._metsigdb_query import count

    spec = MetSigDBQuery(resource=('KEGG',), limit=10)
    assert count(conn, spec) == 4969


# -------------------------------------------------------------------- semantics


def test_resource_filter_narrows_rows(conn):
    rows = fetch(conn, MetSigDBQuery(resource=('KEGG',), limit=50))
    assert rows
    assert {row['resource'] for row in rows} == {'KEGG'}


def test_set_type_filter_narrows_rows(conn):
    rows = fetch(conn, MetSigDBQuery(set_type=('disease',), limit=50))
    assert rows
    assert {row['set_type'] for row in rows} == {'disease'}
    assert {row['resource'] for row in rows} == {'MACdb'}


def test_organism_filter_matches_explicit_values_only(conn):
    """Null-organism rows stay in the dataset and never match the filter."""
    rows = fetch(conn, MetSigDBQuery(organism=9606, limit=50))
    assert rows
    assert all(row['organism'] == 9606 for row in rows)
    assert {row['resource'] for row in rows} == {'Reactome'}


def test_set_source_id_filter_returns_one_set(conn):
    rows = fetch(conn, MetSigDBQuery(set_source_id='R-HSA-1059683', limit=500))
    assert rows
    assert {row['set_source_id'] for row in rows} == {'R-HSA-1059683'}
    assert rows[0]['set_size'] == len(rows)


def test_metabolite_filter_can_span_resources(conn):
    """One metabolite, every set it belongs to, wherever the set came from."""
    probe = fetch(conn, MetSigDBQuery(resource=('MACdb',), limit=1))[0]
    rows = fetch(
        conn,
        MetSigDBQuery(
            metabolite_entity_id=probe['metabolite_entity_id'], limit=MAX_LIMIT
        ),
    )
    assert len(rows) >= 1
    assert {row['metabolite_entity_id'] for row in rows} == {
        probe['metabolite_entity_id']
    }


def test_an_impossible_filter_pair_returns_no_rows(conn):
    """MACdb is the only disease producer, so KEGG has no disease sets."""
    assert fetch(conn, MetSigDBQuery(resource=('KEGG',), set_type=('disease',))) == []


def test_paging_neither_repeats_nor_skips(conn):
    first = fetch(conn, MetSigDBQuery(resource=('KEGG',), limit=100, offset=0))
    second = fetch(conn, MetSigDBQuery(resource=('KEGG',), limit=100, offset=100))
    whole = fetch(conn, MetSigDBQuery(resource=('KEGG',), limit=200, offset=0))

    def identity(rows):
        return [
            (row['resource'], row['set_source_id'], row['metabolite_entity_id'])
            for row in rows
        ]

    assert identity(first) + identity(second) == identity(whole)
    assert len(set(identity(whole))) == len(whole)
