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
    assert set(body) == {'count', 'has_more', 'rows', 'total'}
    assert body['count'] == len(body['rows']) == 5
    assert body['has_more'] is True


def test_the_total_is_absent_unless_asked_for(client):
    """Counting a filter that matches three million rows is opt-in."""
    assert client.get(PATH, params={'limit': 1}).json()['total'] is None
    body = client.get(PATH, params={'resource': 'KEGG', 'limit': 1, 'total': True}).json()
    assert body['total'] == 4969
    assert body['count'] == 1


def test_has_more_ends_a_page_walk(client):
    body = client.get(PATH, params={'set_source_id': 'rn00270'}).json()
    assert body['count'] == 16
    assert body['has_more'] is False


def test_the_sub_type_filter_separates_real_diseases(client):
    rows = _rows(client, set_sub_type='cancer', limit=50)
    assert {row['resource'] for row in rows} == {'MACdb'}
    assert {row['set_sub_type'] for row in rows} == {'cancer'}

    rows = _rows(client, resource='KEGG', set_sub_type='overview_map', limit=50)
    assert {row['set'] for row in rows} <= {
        'rn01100', 'rn01110', 'rn01120', 'rn01200', 'rn01210', 'rn01212',
        'rn01220', 'rn01230', 'rn01232', 'rn01240', 'rn01250',
    }


def test_an_unsupported_sub_type_is_rejected(client):
    assert client.get(PATH, params={'set_sub_type': 'tumour'}).status_code == 400


def test_a_row_satisfies_the_contract(client):
    """The contract has two shapes since cycle 012, and an order they share.

    The default is thirteen fields; `fields=all` is the cycle 010 response. Both
    follow the contract's field order, so a consumer reading positionally is not
    at the mercy of which fields were asked for.
    """
    from omnipath_metabo.server.sets._metsigdb_projection import (
        ALL_FIELDS,
        DEFAULT_FIELDS,
    )

    assert list(_rows(client, limit=1)[0]) == list(DEFAULT_FIELDS)
    assert list(_rows(client, limit=1, fields='all')[0]) == list(ALL_FIELDS)


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
    # `organism` moved to the nine `fields` reaches in cycle 012. Filtering on a
    # field and reading it back are two separate requests now.
    rows = _rows(client, organism=9606, limit=20, fields='organism')
    assert all(row['organism'] == 9606 for row in rows)
    assert {row['resource'] for row in rows} == {'Reactome'}


def test_null_organism_rows_stay_in_the_dataset(client):
    rows = _rows(client, resource='ClassyFire', limit=10, fields='organism')
    assert rows
    assert all(row['organism'] is None for row in rows)


def test_an_empty_result_keeps_the_response_schema(client):
    body = client.get(
        PATH, params={'resource': 'KEGG', 'set_type': 'disease'}
    ).json()
    assert body == {'count': 0, 'has_more': False, 'rows': [], 'total': None}


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
    rows = _rows(client, resource='WikiPathways', limit=5000, fields='organism')
    assert len({row['organism'] for row in rows} - {None}) > 10


# ------------------------------------------------------- the projection (012)


def test_the_default_response_is_the_thirteen_fields(client):
    """Cycle 012: a response a consumer does not have to trim."""
    from omnipath_metabo.server.sets._metsigdb_projection import DEFAULT_FIELDS

    row = _rows(client, limit=1)[0]
    assert set(row) == set(DEFAULT_FIELDS)
    assert 'entity' in row and 'set' in row
    assert 'set_context' not in row
    assert 'provenance_record' not in row


def test_named_fields_are_added_over_the_wire(client):
    row = _rows(client, limit=1, fields='smiles,set_context')[0]
    assert 'smiles' in row
    assert 'set_context' in row
    assert 'entity' in row, 'fields= adds to the default, it does not replace it'


def test_all_returns_every_published_field(client):
    from omnipath_metabo.server.sets._metsigdb_projection import ALL_FIELDS

    row = _rows(client, limit=1, fields='all')[0]
    assert set(row) == set(ALL_FIELDS)
    assert len(row) == 22


def test_an_unknown_field_is_refused_by_name(client):
    """`inchi` left the contract in cycle 010. Asking for it must not pass.

    The rule cycle 010 learned the hard way: an unknown *parameter* returned a
    200 with an unfiltered page, which reads as an answer. An unknown *field*
    would be the same defect in a smaller place.
    """
    response = client.get(PATH, params={'limit': 1, 'fields': 'inchi'})
    assert response.status_code == 400
    assert 'inchi' in response.text


def test_a_column_name_is_not_a_response_name(client):
    """The two renamed fields are reachable under one name only."""
    response = client.get(
        PATH, params={'limit': 1, 'fields': 'metabolite_entity_id'}
    )
    assert response.status_code == 400


def test_fields_is_a_known_parameter(client):
    """It has to be in the allow-list, or the route refuses its own parameter."""
    assert client.get(PATH, params={'limit': 1, 'fields': 'smiles'}).status_code == 200


# ---------------------------------------------------------- grouping (US2)


def test_the_flat_shape_is_still_the_default(client):
    """Cycle 010 consumers are unaffected: grouping is asked for, not imposed."""
    body = client.get(PATH, params={'limit': 2}).json()
    assert set(body) == {'count', 'has_more', 'rows', 'total'}


def test_a_grouped_response_carries_groups_not_rows(client):
    body = client.get(PATH, params={'resource': 'MACdb', 'group': True, 'limit': 3}).json()
    assert set(body) == {'count', 'has_more', 'groups', 'total'}
    assert body['groups'][0]['resource'] == 'MACdb'
    assert body['groups'][0]['sets'][0]['members']


def test_the_grouped_page_counts_sets(client):
    """`limit` means sets here, and `total` counts the sets that matched."""
    body = client.get(
        PATH,
        params={'resource': 'MACdb', 'group': True, 'limit': 4, 'total': True},
    ).json()
    sets = [one for g in body['groups'] for one in g['sets']]
    assert body['count'] == len(sets) <= 4
    assert body['total'] == 269, 'MACdb publishes 269 sets, not 20,291 rows'


def test_a_grouped_set_reports_population_and_returned(client):
    body = client.get(
        PATH, params={'resource': 'KEGG', 'group': True, 'limit': 2}
    ).json()
    for one_set in [one for g in body['groups'] for one in g['sets']]:
        assert one_set['returned'] == len(one_set['members'])
        assert one_set['returned'] == one_set['set_size']


def test_a_grouped_member_does_not_repeat_the_set(client):
    """The set side is written once per set, never on each member."""
    body = client.get(
        PATH, params={'resource': 'KEGG', 'group': True, 'limit': 1}
    ).json()
    member = body['groups'][0]['sets'][0]['members'][0]
    for field in ('resource', 'set', 'set_label', 'set_size'):
        assert field not in member
    assert 'entity' in member
