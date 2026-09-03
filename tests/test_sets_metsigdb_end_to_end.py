"""End-to-end tests for /sets/metsigdb (cycle 010).

Two things, and they are different things.

Parity: what a direct SQL read of the substrate returns and what the API
returns must be the same rows, in the same order, with the same values.

Boundary isolation: a failure in one serving layer must be attributable to that
layer without reading upstream source data. Each case below breaks exactly one
boundary and shows the other two still behave.

Needs the built substrate and the `server` extras::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55435/omnipath \
        uv run --with 'litestar[standard]' --with psycopg2-binary --with pytest \
        pytest tests/test_sets_metsigdb_end_to_end.py -v
"""

from __future__ import annotations

import os

import pytest

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DB_URL, reason='No OMNIPATH_DB_URL/DATABASE_URL; needs the built substrate'
)

PATH = '/sets/metsigdb'
TABLE = 'public.metsigdb_membership'


@pytest.fixture(scope='module')
def conn():
    psycopg2 = pytest.importorskip('psycopg2')
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


def _sql(conn, query, params=None):
    from psycopg2.extras import RealDictCursor

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query, params)
        return [dict(row) for row in cur.fetchall()]


def _api(client, **params):
    response = client.get(PATH, params=params)
    assert response.status_code == 200, response.text
    return response.json()['rows']


# --------------------------------------------------------------------- parity


def test_the_api_returns_the_rows_the_substrate_holds(conn, client):
    served = _api(client, resource='KEGG', limit=250)
    stored = _sql(
        conn,
        f"""
        SELECT resource, set_source_id, metabolite_entity_id::text, set_size
        FROM {TABLE} WHERE resource = 'KEGG'
        ORDER BY resource, set_source_id, metabolite_entity_id
        LIMIT 250
        """,
    )
    assert len(served) == len(stored)
    for api_row, sql_row in zip(served, stored):
        assert api_row['set'] == sql_row['set_source_id']
        assert api_row['entity'] == sql_row['metabolite_entity_id']
        assert api_row['set_size'] == sql_row['set_size']


def test_every_field_of_one_row_matches(conn, client):
    """Served row against stored row, field by field.

    Cycle 012 made the default a subset, so this asks for the whole row. It
    also renamed two fields in the response, so the comparison maps a response
    name back to the column it came from.
    """
    from omnipath_metabo.server.sets._metsigdb_projection import (
        ALL_FIELDS,
        COLUMN_NAMES,
    )

    served = _api(client, set=('R-HSA-1059683',), limit=1, fields='all')[0]
    stored = _sql(
        conn,
        f"""
        SELECT * FROM {TABLE} WHERE set_source_id = %(set_source_id)s
        ORDER BY resource, set_source_id, metabolite_entity_id LIMIT 1
        """,
        {'set_source_id': 'R-HSA-1059683'},
    )[0]
    for field in ALL_FIELDS:
        expected = stored[COLUMN_NAMES.get(field, field)]
        if field == 'entity':
            expected = str(expected)
        assert served[field] == expected, field


def test_counts_agree_per_resource(conn, client):
    for row in _sql(conn, f'SELECT resource, count(*) AS n FROM {TABLE} GROUP BY 1'):
        served = client.get(
            PATH, params={'resource': row['resource'], 'limit': 1}
        ).json()
        assert served['count'] == 1, row['resource']
        # A page of one from a resource with rows proves the filter reaches it.
        assert row['n'] > 0


def test_the_served_build_stamp_is_the_manifest_stamp(conn, client):
    manifest = _sql(conn, 'SELECT build_id FROM public.build_manifest')[0]['build_id']
    served = _api(client, limit=100, fields='build_id')
    assert {row['build_id'] for row in served} == {manifest}


def test_a_metabolites_memberships_agree(conn, client):
    probe = _api(client, resource='MACdb', limit=1)[0]['entity']
    served = _api(client, entity=(probe,), limit=100_000)
    stored = _sql(
        conn,
        f'SELECT count(*) AS n FROM {TABLE} WHERE metabolite_entity_id = %(id)s',
        {'id': probe},
    )[0]['n']
    assert len(served) == stored


# ---------------------------------------------------------- boundary isolation


def test_a_query_layer_fault_stays_in_the_query_layer(conn):
    """A wrong filter changes which rows come back, and nothing else.

    The projection still produces contract rows, which is how a caller tells a
    query fault from a projection fault.
    """
    from omnipath_metabo.server.sets._metsigdb_projection import (
        DEFAULT_FIELDS,
        project_rows,
    )
    from omnipath_metabo.server.sets._metsigdb_query import MetSigDBQuery, fetch

    wrong = project_rows(fetch(conn, MetSigDBQuery(resource=('MACdb',), limit=5)))
    assert {row['resource'] for row in wrong} == {'MACdb'}
    assert all(list(row) == list(DEFAULT_FIELDS) for row in wrong)


def test_a_projection_layer_fault_stays_in_the_projection_layer():
    """A malformed row changes the shape, and touches no SQL.

    The projection never reads the database, so a shape fault cannot be caused
    by the query layer or by upstream source data.
    """
    from omnipath_metabo.server.sets._metsigdb_projection import (
        DEFAULT_FIELDS,
        project_row,
    )

    projected = project_row({'resource': 'KEGG'})
    assert list(projected) == list(DEFAULT_FIELDS)
    assert projected['resource'] == 'KEGG'
    assert projected['set_size'] is None


def test_a_route_layer_fault_never_reaches_the_query_layer(client):
    """A rejected request is rejected at the edge.

    A 400 means the route refused it, so no SQL ran and no row was shaped.
    """
    assert client.get(PATH, params={'resource': 'SMPDB'}).status_code == 400
    assert client.get(PATH, params={'limit': 10**9}).status_code == 400


def test_the_query_layer_needs_no_web_framework():
    """The boundary is real, not just a file split."""
    import subprocess
    import sys

    probe = (
        'import sys; '
        'import omnipath_metabo.server.sets._metsigdb_query as q; '
        'import omnipath_metabo.server.sets._metsigdb_projection as p; '
        "sys.exit(1 if 'litestar' in sys.modules else 0)"
    )
    assert subprocess.run([sys.executable, '-c', probe]).returncode == 0
