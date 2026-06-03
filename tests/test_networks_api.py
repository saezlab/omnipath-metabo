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


def test_network_status_reports_rows_and_build_id(client):
    resp = client.get('/networks/metalinksdb/status')
    assert resp.status_code == 200
    body = resp.json()
    assert body['present'] is True
    assert body['row_count'] > 0
    assert body['build_id']


def test_interactions_consistent_with_combined_contract(client):
    """API row count (unfiltered) matches the combined contract's row count."""
    status = client.get('/networks/liana/status').json()
    resp = client.get('/networks/liana/interactions', params={'limit': 100000})
    assert resp.status_code == 200
    body = resp.json()
    assert body['count'] == status['row_count']
    assert body['rows']


def test_interactions_parquet_format(client):
    resp = client.get(
        '/networks/metalinksdb/interactions',
        params={'limit': 10, 'format': 'parquet'},
    )
    assert resp.status_code == 200
    assert resp.headers['content-type'] == 'application/octet-stream'
    assert resp.content[:4] == b'PAR1'  # parquet magic


def test_resources_lists_sources(client):
    resp = client.get('/networks/metalinksdb/resources')
    assert resp.status_code == 200
    assert len(resp.json()['included_sources']) == 7


def test_unknown_network_is_not_found(client):
    resp = client.get('/networks/no_such_network/status')
    assert resp.status_code == 404
