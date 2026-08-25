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

"""Tests for the uniform /networks API (Milestone G).

Needs a built main DB (with the network views applied) reachable via
OMNIPATH_DB_URL/DATABASE_URL, and the `server` extras (litestar) + psycopg2::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with 'litestar[standard]' --with psycopg2-binary --with pytest \
        pytest tests/test_networks_api.py -v
"""

from __future__ import annotations

import os

import pytest

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DB_URL, reason='No OMNIPATH_DB_URL/DATABASE_URL; needs a built main DB'
)


@pytest.fixture(scope='module')
def client():
    pytest.importorskip('litestar')
    pytest.importorskip('psycopg2')
    from litestar.testing import TestClient

    from omnipath_metabo.server._app import create_app

    app = create_app()
    with TestClient(app=app) as test_client:
        yield test_client


def test_list_networks(client):
    resp = client.get('/networks/')
    assert resp.status_code == 200
    names = {row['name'] for row in resp.json()}
    assert {'metalinksdb', 'liana'} <= names


def _matview_backed(client):
    """The registered datasets this service can still serve rows for.

    A dataset whose registry row names a combined relation has a view here; a
    dataset that names none is a preset and its rows live in the main service.
    Both cycle-008 datasets are presets now, so this returns nothing and the
    row-serving tests skip rather than fail. They stay in the file because the
    route still exists and a dataset onboarded with a view would use it.
    """
    return [
        row['name'] for row in client.get('/networks/').json()
        if row.get('combined_relation')
    ]


@pytest.mark.parametrize('name', ['metalinksdb', 'liana'])
def test_preset_is_reported_present_and_served_elsewhere(client, name):
    """A dataset with no matview is a preset, and the status route says so.

    Its registry row names no schema and no combined relation, which used to
    make the presence probe report a live dataset as absent. Now the row count
    is unknown rather than zero, presence is not in doubt, and asking this
    service for the rows points at the service that has them.

    Both datasets answer this way since MetaLinksDB became a preset. Its views
    are still on disk and this service no longer reaches them: the registry row
    stopped naming them, which is the point — one dataset, one place it is
    served from, and no second copy on a different refresh schedule.
    """
    resp = client.get(f'/networks/{name}/status')
    assert resp.status_code == 200
    body = resp.json()
    assert body['kind'] == 'preset'
    assert body['present'] is True
    assert body['row_count'] is None
    assert body['build_id']

    rows = client.get(f'/networks/{name}/interactions', params={'limit': 10})
    assert rows.status_code == 501
    assert f'/interactions/{name}' in rows.json()['detail']


def test_matview_backed_dataset_serves_its_own_rows(client):
    """A dataset that does own a relation is still served from it."""
    names = _matview_backed(client)
    if not names:
        pytest.skip('every registered dataset is a preset; no rows are served here')
    name = names[0]
    total = client.get(f'/networks/{name}/status').json()['row_count']
    assert total > 3
    limit = 1000
    body = client.get(
        f'/networks/{name}/interactions', params={'limit': limit},
    ).json()
    assert body['count'] == min(total, limit)
    assert body['rows']
    tail = client.get(
        f'/networks/{name}/interactions',
        params={'limit': limit, 'offset': total - 3},
    ).json()
    assert tail['count'] == 3


def test_interactions_parquet_format(client):
    names = _matview_backed(client)
    if not names:
        pytest.skip('every registered dataset is a preset; no rows are served here')
    resp = client.get(
        f'/networks/{names[0]}/interactions',
        params={'limit': 10, 'format': 'parquet'},
    )
    assert resp.status_code == 200
    assert resp.headers['content-type'] == 'application/octet-stream'
    assert resp.content[:4] == b'PAR1'  # parquet magic


def test_resources_lists_sources(client):
    """The resource list is whatever the registry holds, not a fixed number.

    It was pinned at seven and the dataset has carried twelve since the
    resource expansion, so the number was wrong rather than protective. What
    the route has to do is report the registry faithfully.
    """
    resp = client.get('/networks/metalinksdb/resources')
    assert resp.status_code == 200
    served = resp.json()['included_sources']
    assert served, 'the dataset is reported as contributing from no source'
    registered = next(
        row for row in client.get('/networks/').json()
        if row['name'] == 'metalinksdb'
    )['included_sources']
    assert served == registered


def test_unknown_network_is_not_found(client):
    resp = client.get('/networks/no_such_network/status')
    assert resp.status_code == 404
