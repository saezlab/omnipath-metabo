"""Exclusion of metabolites from /sets/metsigdb (cycle 013).

``exclude_metabolites`` leaves the metabolites of named classes out of the
response. The membership table keeps every row; the exclusion is applied to the
query, against the ``exclusion_metabolite`` table that ``omnipath-build``
publishes beside it.

The data tests run against the real membership table through a scratch schema.
The schema holds a view over ``public.metsigdb_membership`` and a small
exclusion table built from real entities, so nothing in ``public`` is touched.

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55435/omnipath \
        uv run --with 'litestar[standard]' --with psycopg2-binary --with pytest \
        pytest tests/test_sets_metsigdb_exclusion.py -v
"""

from __future__ import annotations

import os
import re
import uuid

import pytest

from omnipath_metabo.server.sets import _metsigdb_query as query
from omnipath_metabo.server.sets._metsigdb_query import (
    ROW_CEILING,
    MetSigDBQuery,
    build_query,
    count_query,
    exclusion_classes,
    set_page_query,
)

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')

PATH = '/sets/metsigdb'
CLASSES = ('ions', 'atoms', 'small_molecules', 'xenobiotics', 'cofactors')

# A small KEGG set. Every member goes into `xenobiotics`, a class nothing else
# uses, so that excluding it empties the set and leaves the rest alone.
EMPTIED_SET = 'rn00270'

needs_database = pytest.mark.skipif(
    not DB_URL,
    reason='No OMNIPATH_DB_URL/DATABASE_URL; needs the built substrate',
)


# ------------------------------------------------------------ the class list


def test_no_parameter_and_none_exclude_nothing():
    assert exclusion_classes(()) == ()
    assert exclusion_classes(('none',)) == ()
    assert exclusion_classes(('NONE',)) == ()


def test_all_expands_to_every_class():
    assert exclusion_classes(('all',)) == CLASSES
    assert exclusion_classes(('All',)) == CLASSES


def test_a_class_list_matches_without_regard_to_case_and_answers_in_published_spelling():
    assert exclusion_classes(('IONS', 'Atoms')) == ('ions', 'atoms')
    assert exclusion_classes(('ions', 'ions')) == ('ions',)


@pytest.mark.parametrize(
    'values',
    [('none', 'ions'), ('all', 'ions'), ('ions', 'all'), ('all', 'none')],
)
def test_none_and_all_cannot_be_combined_with_a_class(values):
    with pytest.raises(ValueError, match='none|all'):
        exclusion_classes(values)


def test_an_unknown_class_is_refused_naming_it_and_the_valid_ones():
    with pytest.raises(ValueError) as error:
        exclusion_classes(('ions', 'lipids'))
    message = str(error.value)
    assert 'lipids' in message
    for name in ('none', 'all', *CLASSES):
        assert name in message


def test_xenobiotics_is_a_valid_class():
    assert exclusion_classes(('xenobiotics',)) == ('xenobiotics',)


def test_an_unfiltered_query_is_the_cycle_twelve_query():
    spec = MetSigDBQuery(resource=('KEGG',))
    sql, params = build_query(spec)
    assert 'exclusion' not in sql
    assert 'exclude_metabolites' not in params
    assert build_query(
        MetSigDBQuery(resource=('KEGG',), exclude_metabolites=())
    ) == (sql, params)
    assert count_query(spec)[0].count('exclusion') == 0


def test_an_exclusion_is_an_anti_join_on_the_selected_classes():
    spec = MetSigDBQuery(exclude_metabolites=('ions', 'atoms'))
    sql, params = build_query(spec)
    assert 'NOT EXISTS' in sql
    assert 'exclusion_metabolite' in sql
    assert params['exclude_metabolites'] == ['ions', 'atoms']


def test_the_count_applies_the_same_exclusion_as_the_page():
    spec = MetSigDBQuery(resource=('KEGG',), exclude_metabolites=('cofactors',))
    sql, params = count_query(spec)
    assert 'NOT EXISTS' in sql
    assert params['exclude_metabolites'] == ['cofactors']
    assert 'limit' not in params


def test_the_grouped_set_page_keeps_every_candidate_set():
    spec = MetSigDBQuery(resource=('KEGG',), exclude_metabolites=('ions',))
    sql, params = set_page_query(spec)
    # The exclusion thins the count and never the set list.
    where = sql.partition(f'FROM {query.membership_table()}')[2].partition(
        ' GROUP BY'
    )[0]
    assert 'exclusion_metabolite' not in where
    assert 'FILTER' in sql
    assert params['exclude_metabolites'] == ['ions']


def test_emptied_sets_count_toward_the_page_and_not_toward_the_row_ceiling():
    rows = [
        {'resource': 'KEGG', 'set_source_id': 'a', 'returned': ROW_CEILING},
        {'resource': 'KEGG', 'set_source_id': 'b', 'returned': 0},
        {'resource': 'KEGG', 'set_source_id': 'c', 'returned': 0},
    ]
    taken, has_more = query._fit_the_ceiling(rows, 3)
    assert [row['set_source_id'] for row in taken] == ['a', 'b', 'c']
    assert has_more is False
    taken, has_more = query._fit_the_ceiling(rows, 2)
    assert [row['set_source_id'] for row in taken] == ['a', 'b']
    assert has_more is True


# -------------------------------------------------------------- the schema


def test_the_default_schema_leaves_every_existing_query_unchanged(monkeypatch):
    monkeypatch.delenv('OMNIPATH_PG_SCHEMA', raising=False)
    spec = MetSigDBQuery(resource=('KEGG',))
    sql, _ = build_query(spec)
    assert sql == (
        f'SELECT {", ".join(spec.columns())} FROM public.metsigdb_membership '
        'WHERE resource = ANY(%(resource)s) '
        'ORDER BY resource, set_source_id, metabolite_entity_id '
        'LIMIT %(limit)s OFFSET %(offset)s'
    )


def test_one_schema_setting_names_both_tables(monkeypatch):
    monkeypatch.setenv('OMNIPATH_PG_SCHEMA', 'release_7')
    sql, _ = build_query(MetSigDBQuery(exclude_metabolites=('ions',)))
    assert 'FROM release_7.metsigdb_membership' in sql
    assert 'release_7.exclusion_metabolite' in sql
    assert 'public.' not in sql


@pytest.mark.parametrize(
    'schema', ['x; DROP TABLE y', 'a b', '', '1abc', 'a-b']
)
def test_a_schema_name_that_is_not_an_identifier_is_refused(
    monkeypatch, schema
):
    monkeypatch.setenv('OMNIPATH_PG_SCHEMA', schema)
    with pytest.raises(ValueError, match='OMNIPATH_PG_SCHEMA'):
        build_query(MetSigDBQuery())


# ------------------------------------------------------- data and the route


@pytest.fixture(scope='module')
def schemas():
    """Two scratch schemas: one with an exclusion table, one without."""
    psycopg2 = pytest.importorskip('psycopg2')
    if not DB_URL:
        pytest.skip('no OMNIPATH_DB_URL/DATABASE_URL')
    suffix = uuid.uuid4().hex[:12]
    full, bare = f'exclusion_full_{suffix}', f'exclusion_bare_{suffix}'
    connection = psycopg2.connect(DB_URL)
    connection.autocommit = True
    try:
        with connection.cursor() as cur:
            for schema in (full, bare):
                cur.execute(f'CREATE SCHEMA {schema}')
                cur.execute(
                    f'CREATE VIEW {schema}.metsigdb_membership AS '
                    'SELECT * FROM public.metsigdb_membership'
                )
            cur.execute(
                f"""
                CREATE TABLE {full}.exclusion_metabolite (
                  entity_id uuid NOT NULL,
                  exclusion_class text NOT NULL,
                  assignment_source text NOT NULL,
                  origin text NOT NULL,
                  match_kind text NOT NULL,
                  seed_namespace text,
                  seed_identifier text,
                  evidence text NOT NULL,
                  via_entity_id uuid,
                  list_version text NOT NULL,
                  build_id text NOT NULL,
                  PRIMARY KEY (entity_id, exclusion_class, assignment_source)
                );
                CREATE INDEX ON {full}.exclusion_metabolite (exclusion_class, entity_id);
                INSERT INTO {full}.exclusion_metabolite
                  (entity_id, exclusion_class, assignment_source, origin, match_kind,
                   evidence, list_version, build_id)
                SELECT DISTINCT m.metabolite_entity_id, v.cls, 'metaproviz_seed', 'curated',
                       'entity', 'fixture', 'fixture', 'fixture'
                FROM public.metsigdb_membership m
                JOIN (VALUES ('15378', 'ions'), ('15377', 'small_molecules'),
                             ('16526', 'small_molecules'), ('30616', 'cofactors'),
                             ('57540', 'cofactors')) AS v(chebi, cls) ON m.chebi = v.chebi;
                INSERT INTO {full}.exclusion_metabolite
                  (entity_id, exclusion_class, assignment_source, origin, match_kind,
                   evidence, list_version, build_id)
                SELECT DISTINCT metabolite_entity_id, 'xenobiotics', 'metaproviz_seed',
                       'curated', 'entity', 'fixture', 'fixture', 'fixture'
                FROM public.metsigdb_membership
                WHERE resource = 'KEGG' AND set_source_id = '{EMPTIED_SET}';
                ANALYZE {full}.exclusion_metabolite;
                """
            )
        yield full, bare
    finally:
        with connection.cursor() as cur:
            for schema in (full, bare):
                cur.execute(f'DROP SCHEMA {schema} CASCADE')
        connection.close()


@pytest.fixture
def full_schema(schemas, monkeypatch):
    monkeypatch.setenv('OMNIPATH_PG_SCHEMA', schemas[0])
    return schemas[0]


@pytest.fixture
def bare_schema(schemas, monkeypatch):
    monkeypatch.setenv('OMNIPATH_PG_SCHEMA', schemas[1])
    return schemas[1]


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


@pytest.fixture(scope='module')
def client():
    pytest.importorskip('litestar')
    pytest.importorskip('psycopg2')
    from litestar.testing import TestClient

    from omnipath_metabo.server._app import create_app

    with TestClient(app=create_app()) as test_client:
        yield test_client


def excluded(conn, schema, *classes):
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT DISTINCT entity_id::text FROM {schema}.exclusion_metabolite '
            'WHERE exclusion_class = ANY(%s)',
            [list(classes)],
        )
        return {row[0] for row in cur.fetchall()}


def entity_ids(conn, chebi):
    with conn.cursor() as cur:
        cur.execute(
            'SELECT DISTINCT metabolite_entity_id::text FROM public.metsigdb_membership '
            'WHERE chebi = %s',
            [chebi],
        )
        return [row[0] for row in cur.fetchall()]


def rows(client, **params):
    response = client.get(PATH, params=params)
    assert response.status_code == 200, response.text
    return response.json()


def entities(body):
    return {row['entity'] for row in body['rows']}


# ----------------------------------------------------------------- refusals


@needs_database
def test_an_unknown_class_is_refused_with_400_naming_the_value_and_the_valid_classes(
    client, full_schema
):
    response = client.get(PATH, params={'exclude_metabolites': 'lipids'})
    assert response.status_code == 400
    detail = response.json()['detail']
    assert 'lipids' in detail
    for name in ('none', 'all', *CLASSES):
        assert name in detail


@needs_database
def test_none_with_a_class_is_refused(client, full_schema):
    response = client.get(PATH, params={'exclude_metabolites': 'none,ions'})
    assert response.status_code == 400


@needs_database
def test_the_parameter_is_accepted_and_listed_in_the_unsupported_parameter_message(
    client, full_schema
):
    assert (
        client.get(
            PATH, params={'limit': 1, 'exclude_metabolites': 'none'}
        ).status_code
        == 200
    )
    response = client.get(PATH, params={'exclude_metaboltes': 'all'})
    assert response.status_code == 400
    assert 'exclude_metabolites' in response.json()['detail']


@needs_database
def test_none_answers_exactly_as_no_parameter_does(client, full_schema):
    params = {'resource': 'KEGG', 'limit': 25, 'total': True}
    assert rows(client, **params) == rows(
        client, exclude_metabolites='none', **params
    )


@needs_database
def test_xenobiotics_is_accepted_and_removes_only_the_rows_of_its_members(
    client, full_schema, conn
):
    listed = excluded(conn, full_schema, 'xenobiotics')
    with conn.cursor() as cur:
        cur.execute(
            'SELECT count(*) FROM public.metsigdb_membership '
            "WHERE resource = 'KEGG' AND metabolite_entity_id::text = ANY(%s)",
            [list(listed)],
        )
        removed = cur.fetchone()[0]
    params = {'resource': 'KEGG', 'limit': 1, 'total': True}
    body = rows(client, exclude_metabolites='xenobiotics', **params)
    assert body['total'] == rows(client, **params)['total'] - removed


# ------------------------------------------------------------ flat responses


@needs_database
def test_all_leaves_no_returned_member_that_is_in_the_table(
    client, full_schema, conn
):
    listed = excluded(conn, full_schema, *CLASSES)
    assert listed
    body = rows(
        client, resource='Reactome', exclude_metabolites='all', limit=5000
    )
    assert body['count'] > 0
    assert not entities(body) & listed


@needs_database
def test_a_named_class_removes_that_class_and_keeps_the_others(
    client, full_schema, conn
):
    atp = entity_ids(conn, '30616')
    proton = entity_ids(conn, '15378')
    assert atp and proton
    kept = rows(client, entity=atp, exclude_metabolites='ions,atoms', limit=50)
    assert kept['count'] > 0
    gone = rows(client, entity=atp, exclude_metabolites='cofactors', limit=50)
    assert gone['count'] == 0
    gone = rows(client, entity=proton, exclude_metabolites='ions', limit=50)
    assert gone['count'] == 0
    kept = rows(
        client, entity=proton, exclude_metabolites='cofactors', limit=50
    )
    assert kept['count'] > 0


@needs_database
def test_a_repeated_parameter_is_the_same_as_a_comma_list(client, full_schema):
    params = {'resource': 'Reactome', 'limit': 100, 'total': True}
    commas = rows(client, exclude_metabolites='ions,cofactors', **params)
    repeated = client.get(
        PATH, params={**params, 'exclude_metabolites': ['ions', 'cofactors']}
    ).json()
    assert commas == repeated


@needs_database
def test_the_total_and_has_more_describe_the_filtered_result(
    client, full_schema, conn
):
    listed = excluded(conn, full_schema, 'small_molecules')
    with conn.cursor() as cur:
        cur.execute(
            'SELECT count(*), count(*) FILTER (WHERE metabolite_entity_id::text = ANY(%s)) '
            "FROM public.metsigdb_membership WHERE resource = 'Reactome'",
            [list(listed)],
        )
        everything, removed = cur.fetchone()
    assert removed > 0
    body = rows(
        client,
        resource='Reactome',
        exclude_metabolites='small_molecules',
        limit=1,
        total=True,
    )
    assert body['total'] == everything - removed
    assert body['has_more'] is True

    # A set the exclusion empties has no row at all: the page ends where it should.
    emptied = rows(
        client, set=EMPTIED_SET, exclude_metabolites='xenobiotics', total=True
    )
    assert (emptied['count'], emptied['total'], emptied['has_more']) == (
        0,
        0,
        False,
    )
    assert emptied['rows'] == []


@needs_database
def test_asking_for_an_excluded_metabolite_returns_an_empty_result_and_no_error(
    client, full_schema, conn
):
    atp = entity_ids(conn, '30616')
    body = rows(client, entity=atp, exclude_metabolites='all', total=True)
    assert (body['count'], body['total'], body['rows']) == (0, 0, [])


# --------------------------------------------------------- grouped responses


def grouped(client, **params):
    return rows(client, group=True, **params)


def only_sets(body):
    return [one for group in body['groups'] for one in group['sets']]


@needs_database
def test_a_set_that_loses_every_member_is_returned_empty_with_its_population(
    client, full_schema
):
    plain = only_sets(grouped(client, set=EMPTIED_SET))
    assert len(plain) == 1 and plain[0]['returned'] == plain[0]['set_size'] > 0

    filtered = only_sets(
        grouped(client, set=EMPTIED_SET, exclude_metabolites='xenobiotics')
    )
    assert len(filtered) == 1
    [one] = filtered
    assert one['members'] == []
    assert one['returned'] == 0
    assert one['set_size'] == plain[0]['set_size']


@needs_database
def test_returned_counts_the_surviving_members_and_set_size_stays_the_population(
    client, full_schema, conn
):
    listed = excluded(conn, full_schema, 'small_molecules', 'ions')
    body = grouped(
        client,
        resource='Reactome',
        exclude_metabolites='small_molecules,ions',
        limit=40,
    )
    thinned = False
    for one in only_sets(body):
        assert one['returned'] == len(one['members'])
        assert not {m['entity'] for m in one['members']} & listed
        thinned |= one['returned'] < one['set_size']
    assert thinned, 'expected at least one Reactome set to lose a member'


@needs_database
def test_the_candidate_sets_and_the_grouped_total_do_not_depend_on_the_exclusion(
    client, full_schema
):
    plain = grouped(client, resource='MACdb', limit=30, total=True)
    filtered = grouped(
        client,
        resource='MACdb',
        limit=30,
        total=True,
        exclude_metabolites='all',
    )
    assert filtered['total'] == plain['total'] == 269
    assert [s['set'] for s in only_sets(filtered)] == [
        s['set'] for s in only_sets(plain)
    ]


@needs_database
def test_a_paged_walk_sees_every_candidate_set_exactly_once(
    client, full_schema
):
    seen, offset = [], 0
    while True:
        body = grouped(
            client,
            resource='MACdb',
            limit=60,
            offset=offset,
            exclude_metabolites='all',
        )
        seen += [s['set'] for s in only_sets(body)]
        if not body['has_more']:
            break
        offset += 60
    assert len(seen) == len(set(seen)) == 269


# --------------------------------------------------------------- absent table


@needs_database
@pytest.mark.parametrize('value', ['all', 'ions', 'ions,cofactors'])
def test_a_request_for_an_exclusion_without_the_table_gets_503(
    client, bare_schema, value
):
    response = client.get(
        PATH, params={'limit': 1, 'exclude_metabolites': value}
    )
    assert response.status_code == 503
    assert 'exclusion list' in response.json()['detail']
    assert 'not been built' in response.json()['detail']


@needs_database
@pytest.mark.parametrize('params', [{}, {'exclude_metabolites': 'none'}])
def test_a_request_without_an_exclusion_never_reads_the_table(
    client, bare_schema, params
):
    assert client.get(PATH, params={'limit': 3, **params}).status_code == 200
    body = client.get(PATH, params={'limit': 3, 'group': True, **params})
    assert body.status_code == 200


@needs_database
def test_an_invalid_class_is_a_400_even_when_the_table_is_absent(
    client, bare_schema
):
    assert (
        client.get(PATH, params={'exclude_metabolites': 'lipids'}).status_code
        == 400
    )


# -------------------------------------------------------------- no full scans


@needs_database
@pytest.mark.parametrize(
    ('label', 'spec'),
    [
        (
            'resource and every class',
            MetSigDBQuery(
                resource=('ClassyFire',), exclude_metabolites=CLASSES
            ),
        ),
        (
            'resource and two classes',
            MetSigDBQuery(
                resource=('ClassyFire',), exclude_metabolites=('ions', 'atoms')
            ),
        ),
        (
            'reactome and cofactors',
            MetSigDBQuery(
                resource=('Reactome',), exclude_metabolites=('cofactors',)
            ),
        ),
        (
            'no other filter',
            MetSigDBQuery(exclude_metabolites=('small_molecules',)),
        ),
    ],
)
def test_an_exclusion_does_not_make_the_page_read_the_substrate(
    conn, full_schema, label, spec
):
    """The anti-join must stay a cheap probe per row of the page.

    The same rule as the other filters: a plan that discards a hundred thousand
    rows to fill one page has read the table, whatever the node is called.
    """
    sql, params = build_query(spec)
    with conn.cursor() as cur:
        cur.execute(f'EXPLAIN (ANALYZE, TIMING OFF, FORMAT TEXT) {sql}', params)
        plan = '\n'.join(row[0] for row in cur.fetchall())

    assert 'Seq Scan on metsigdb_membership' not in plan, f'{label}:\n{plan}'
    discarded = sum(
        int(n) for n in re.findall(r'Rows Removed by Filter: (\d+)', plan)
    )
    assert discarded < 100_000, f'{label} discarded {discarded:,} rows:\n{plan}'
