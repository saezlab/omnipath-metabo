"""Route-layer tests for /sets/metsigdb (cycle 010).

The route layer owns request validation, the page bounds, and dispatch. It owns
neither which rows come back nor their shape.

Needs the built substrate and the `server` extras::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55435/omnipath \
        uv run --with 'litestar[standard]' --with psycopg2-binary --with pytest \
        pytest tests/test_sets_metsigdb_route.py -v
"""

from __future__ import annotations

import os

import pytest

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DB_URL, reason='No OMNIPATH_DB_URL/DATABASE_URL; needs the built substrate'
)

PATH = '/sets/metsigdb'


@pytest.fixture(scope='module')
def client():
    pytest.importorskip('litestar')
    pytest.importorskip('psycopg2')
    from litestar.testing import TestClient

    from omnipath_metabo.server._app import create_app

    with TestClient(app=create_app()) as test_client:
        yield test_client


def _rows(client, **params):
    response = client.get(PATH, params=params)
    assert response.status_code == 200, response.text
    return response.json()['rows']


# ------------------------------------------------------------------- the family


def test_the_route_family_is_registered(client):
    assert client.get(PATH, params={'limit': 1}).status_code == 200


def test_the_response_carries_its_rows_and_a_count(client):
    body = client.get(PATH, params={'limit': 5}).json()
    assert set(body) == {'count', 'rows'}
    assert body['count'] == len(body['rows']) == 5


def test_a_row_satisfies_the_contract(client):
    from omnipath_metabo.server.sets._metsigdb_projection import ROW_FIELDS

    row = _rows(client, limit=1)[0]
    assert list(row) == list(ROW_FIELDS)


# ------------------------------------------------------------------------ paging


def test_the_default_page_size_is_the_shared_default(client):
    from omnipath_metabo.server.sets._metsigdb_query import DEFAULT_LIMIT

    assert len(_rows(client)) == DEFAULT_LIMIT


def test_a_page_beyond_the_maximum_is_rejected_not_served(client):
    from omnipath_metabo.server.sets._metsigdb_query import MAX_LIMIT

    assert client.get(PATH, params={'limit': MAX_LIMIT + 1}).status_code == 400


def test_a_nonsense_page_is_rejected(client):
    assert client.get(PATH, params={'limit': 0}).status_code == 400
    assert client.get(PATH, params={'offset': -1}).status_code == 400


def test_paging_is_stable_across_requests(client):
    first = _rows(client, resource='KEGG', limit=50, offset=0)
    second = _rows(client, resource='KEGG', limit=50, offset=50)
    whole = _rows(client, resource='KEGG', limit=100, offset=0)
    assert first + second == whole


# ----------------------------------------------------------------------- filters


def test_supported_filters_are_accepted(client):
    rows = _rows(client, resource='MACdb', set_type='disease', limit=10)
    assert {row['resource'] for row in rows} == {'MACdb'}
    assert {row['set_type'] for row in rows} == {'disease'}


def test_a_filter_takes_several_values(client):
    rows = _rows(client, resource=['KEGG', 'Reactome'], limit=200)
    assert {row['resource'] for row in rows} <= {'KEGG', 'Reactome'}


def test_an_unsupported_filter_name_is_rejected(client):
    """Silently ignoring a filter answers the wrong question with a 200.

    `?hmdb=HMDB00077` reads as "the memberships of this metabolite". Ignored,
    it returns the unfiltered first page, which looks like an answer and is
    not one. The identifier columns are published but not filterable in v1, so
    naming one has to fail loudly.
    """
    for param in ('hmdb', 'inchikey', 'chebi', 'metabolite_label', 'nonsense'):
        response = client.get(PATH, params={param: 'x'})
        assert response.status_code == 400, param
        assert param in response.text


def test_a_rejected_filter_names_what_is_supported(client):
    body = client.get(PATH, params={'hmdb': 'HMDB00077'}).text
    for supported in ('resource', 'set_type', 'organism', 'set_source_id'):
        assert supported in body


def test_an_unsupported_filter_value_is_rejected(client):
    assert client.get(PATH, params={'resource': 'SMPDB'}).status_code == 400
    assert client.get(PATH, params={'set_type': 'protein_association'}).status_code == 400


def test_the_organism_filter_matches_explicit_values_only(client):
    rows = _rows(client, organism=9606, limit=20)
    assert all(row['organism'] == 9606 for row in rows)
    assert {row['resource'] for row in rows} == {'Reactome'}


def test_null_organism_rows_stay_in_the_dataset(client):
    rows = _rows(client, resource='ClassyFire', limit=10)
    assert rows
    assert all(row['organism'] is None for row in rows)


def test_an_empty_result_keeps_the_response_schema(client):
    body = client.get(
        PATH, params={'resource': 'KEGG', 'set_type': 'disease'}
    ).json()
    assert body == {'count': 0, 'rows': []}


def test_a_set_reports_its_own_size(client):
    rows = _rows(client, set_source_id='R-HSA-1059683', limit=1000)
    assert rows
    assert rows[0]['set_size'] == len(rows)


# ------------------------------------------------------------------ the contract


def test_no_stats_or_discovery_endpoints_exist(client):
    """v1 publishes one route family and no helpers."""
    for path in (f'{PATH}/stats', f'{PATH}/resources', f'{PATH}/sets', '/sets'):
        assert client.get(path).status_code in (404, 405)


def test_named_and_unnamed_resources_both_serve(client):
    """ClassyFire names every set; KEGG names none. Both are contract-valid."""
    assert all(row['set_label'] for row in _rows(client, resource='ClassyFire', limit=50))
    assert all(row['set_label'] is None for row in _rows(client, resource='KEGG', limit=50))


def test_wikipathways_serves_many_species(client):
    rows = _rows(client, resource='WikiPathways', limit=5000)
    assert len({row['organism'] for row in rows} - {None}) > 10
