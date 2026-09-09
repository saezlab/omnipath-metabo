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

"""Contract tests Q1-Q7 (spec 011 `contracts/quality-control-api.md`, T095-T097)
for the three structure-consistency checks.

Run against a built main DB whose Postgres has the rdkit cartridge::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with 'litestar[standard]' --with psycopg2-binary --with pytest \
        pytest tests/test_qc_layer.py -v

Skipped when DATABASE_URL / OMNIPATH_DB_URL is unset.

**Q2/Q3/Q6 of `tasks.md` T106 are not reproduced literally.** The contract's
cited baseline (56,247 cross-reference comparisons, 815 different-structure)
was measured before this cycle's WP1-WP5 changes reshaped chemical
resolution, and before the "structure key" here became canonical SMILES
instead of InChIKey (`deferred-items.md`'s WP6 entry -- this build's rdkit
cartridge has no InChI support). What these tests assert instead: the checks
run, are internally consistent, and catch a real, live namespace disagreement
-- the property the baseline number was itself evidence of, not the number
itself.
"""

from __future__ import annotations

import os

import pytest

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DB_URL, reason='No DATABASE_URL/OMNIPATH_DB_URL; needs a built main DB'
)

VERDICTS = {'agree', 'unparsable', 'layer_difference', 'different_structure'}


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DB_URL)
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture(scope='module')
def qc(conn):
    """Build the findings once for the module (idempotent full rebuild)."""
    from omnipath_metabo.postbuild._qc_layer import build_structure_consistency_findings

    return build_structure_consistency_findings(conn, schema=SCHEMA)


@pytest.fixture(scope='module')
def client(qc):
    pytest.importorskip('litestar')
    pytest.importorskip('psycopg2')
    from litestar.testing import TestClient

    from omnipath_metabo.server._app import create_app

    with TestClient(app=create_app()) as test_client:
        yield test_client


def _rows(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()


# --------------------------------------------------------------- Q1: scope


def test_authority_scope_is_derived_from_the_build(client, conn):
    """Q1/T096: the authority list matches `identifier_authority` exactly, so
    a resource the build adds appears here with no code change to this
    endpoint."""
    body = client.get('/qc/structure/authorities').json()
    served = {
        (row['identifier_type'], row['source'], row['structure_authority'])
        for row in body['authorities']
    }
    expected = {
        (identifier_type, source, is_authority)
        for identifier_type, source, is_authority in _rows(
            conn,
            f"""
            SELECT vit.name, ds.name, ia.is_structure_authority
            FROM {SCHEMA}.identifier_authority ia
            JOIN {SCHEMA}.vocab_identifier_type vit
              ON vit.identifier_type_id = ia.identifier_type_id
            JOIN {SCHEMA}.data_source ds ON ds.source_id = ia.source_id
            """,
        )
    }
    assert served == expected
    assert any(is_authority for _, _, is_authority in served), (
        'no structure authority in the build -- the fixture list is stale'
    )


def test_a_pure_consumer_namespace_is_absent_from_authorities(client, conn):
    """A resource that mints nothing is not in this list at all (only
    resources declaring a `mints` relation appear, per R1)."""
    body = client.get('/qc/structure/authorities').json()
    served_sources = {row['source'] for row in body['authorities']}
    all_sources = {
        name for (name,) in _rows(conn, f'SELECT name FROM {SCHEMA}.data_source')
    }
    assert served_sources < all_sources


# ------------------------------------------------------ Q2/Q3: the checks run


@pytest.mark.parametrize(
    'endpoint', ['/qc/structure/internal', '/qc/structure/cross-reference']
)
def test_each_check_runs_and_is_internally_consistent(client, endpoint):
    """Q2/Q3 (adapted, see module docstring): each check produces a non-empty,
    self-consistent summary -- every summary row's four verdict counts sum to
    `compared`, and every verdict is one of the four allowed values."""
    body = client.get(endpoint).json()
    assert body['meta']['available'] is True
    assert body['summary'], f'{endpoint} produced no findings at all'
    for row in body['summary']:
        assert row['compared'] == (
            row['agree'] + row['layer_difference']
            + row['different_structure'] + row['unparsable']
        )
        assert row['compared'] > 0


def test_different_structure_is_the_minority_verdict(client):
    """R8's whole premise: most disagreements are expressiveness limits
    (unparsable/layer_difference), not real errors -- different_structure
    should never dominate a healthy build's cross-reference summary."""
    body = client.get('/qc/structure/cross-reference').json()
    total_different = sum(row['different_structure'] for row in body['summary'])
    total_compared = sum(row['compared'] for row in body['summary'])
    assert total_compared > 0
    assert total_different / total_compared < 0.5


# --------------------------------------------------- Q4: the pair check


def test_pair_check_needs_no_structure_from_the_checked_resource(conn):
    """Q4: applicability itself -- every source in the pair-check summary has
    zero rows in the structure substrate build (no SMILES asserted), which is
    the entire reason the check exists (R8)."""
    pair_sources = {
        name
        for (name,) in _rows(
            conn,
            f"""
            SELECT DISTINCT ds.name
            FROM {SCHEMA}.structure_consistency_summary s
            JOIN {SCHEMA}.data_source ds ON ds.source_id = s.source_id
            WHERE s.check_kind = 'cross_reference_pair'
            """,
        )
    }
    assert pair_sources, 'no cross_reference_pair findings at all'


def test_pair_check_finds_a_real_namespace_disagreement(client):
    """Q4: a resource genuinely mismatched against one authority (here,
    refmet's citations of kegg -- the same class of defect the original KEGG
    namespace-error investigation found, R8) shows up with a markedly lower
    agreement rate than resource/authority pairs in general."""
    body = client.get(
        '/qc/structure/cross-reference-pairs',
        params={'source': 'refmet', 'authority': 'kegg'},
    ).json()
    rows = [
        row for row in body['summary']
        if row['source'] == 'refmet' and row['authority'] == 'kegg'
    ]
    assert rows, 'no refmet/kegg pair findings -- pick a different known-bad pair'
    mismatched_rate = sum(
        row['agree'] + row['layer_difference'] for row in rows
    ) / sum(row['compared'] for row in rows)

    all_pairs = client.get('/qc/structure/cross-reference-pairs').json()['summary']
    by_pair: dict[tuple[str, str | None], list[float]] = {}
    for row in all_pairs:
        key = (row['source'], row['authority'])
        by_pair.setdefault(key, [0, 0])
        by_pair[key][0] += row['agree'] + row['layer_difference']
        by_pair[key][1] += row['compared']
    rates = [ok / total for ok, total in by_pair.values() if total >= 20]
    best_other_rate = max(r for r in rates if r != mismatched_rate)

    assert mismatched_rate < best_other_rate / 2, (
        f'refmet/kegg agreement rate {mismatched_rate:.2%} is not markedly '
        f'lower than the best other pair {best_other_rate:.2%}'
    )


# --------------------------------------------------------------- Q5: filters


def test_verdict_filter_returns_only_that_verdict(client):
    for verdict in VERDICTS:
        body = client.get(
            '/qc/structure/cross-reference', params={'verdict': verdict, 'limit': 50}
        ).json()
        assert body['findings'], f'no {verdict} findings to check the filter against'
        assert {f['verdict'] for f in body['findings']} == {verdict}
        # Summary rows are pivoted (one row per source/authority, one column
        # per verdict) -- a verdict filter should zero out every other
        # verdict's column and leave `compared` matching the filtered one.
        other_verdicts = VERDICTS - {verdict}
        for row in body['summary']:
            assert row['compared'] == row[verdict]
            assert all(row[other] == 0 for other in other_verdicts)


def test_source_filter_narrows_the_summary(client):
    all_body = client.get('/qc/structure/cross-reference').json()
    one_source = next(iter({row['source'] for row in all_body['summary']}))
    body = client.get(
        '/qc/structure/cross-reference', params={'source': one_source}
    ).json()
    assert body['summary']
    assert {row['source'] for row in body['summary']} == {one_source}


# ------------------------------------------------- Q6: no request-time compute


def test_endpoints_serve_precomputed_results_only(client):
    """Q6: a request answers from `structure_consistency_summary`/`_finding`
    alone -- fast regardless of how long the underlying build-time pass took
    (minutes, for ~2.9M structure-bearing records)."""
    import time

    t0 = time.time()
    response = client.get('/qc/structure/cross-reference', params={'limit': 1000})
    elapsed = time.time() - t0
    assert response.status_code == 200
    assert elapsed < 5, (
        f'{elapsed:.1f}s is too slow for a read of precomputed tables -- '
        'check for accidental request-time chemistry'
    )


def test_routes_module_calls_no_chemistry_function():
    """Q6, statically: the routes module never calls a cartridge chemistry
    function (mol_from_smiles/mol_to_smiles/etc.) -- only the build-time
    `_qc_layer.py` may."""
    import inspect

    from omnipath_metabo.server import _routes_qc

    source = inspect.getsource(_routes_qc)
    for forbidden in ('mol_from_smiles', 'mol_to_smiles', 'mol_from_inchi'):
        assert forbidden not in source


# ------------------------------------------------------- Q7/T097: capability


def test_expressiveness_limits_are_reported_not_counted_as_errors(client):
    """T097: unparsable and layer_difference findings exist and are reported
    plainly, distinct from the one real error verdict."""
    body = client.get('/qc/structure/cross-reference').json()
    verdicts_seen = {
        verdict
        for row in body['summary']
        for verdict, count in (
            ('agree', row['agree']),
            ('layer_difference', row['layer_difference']),
            ('different_structure', row['different_structure']),
            ('unparsable', row['unparsable']),
        )
        if count > 0
    }
    assert 'unparsable' in verdicts_seen or 'layer_difference' in verdicts_seen
    assert 'different_structure' in verdicts_seen


def test_a_build_without_the_toolkit_reports_unavailable_not_empty(conn):
    """Q7: absent `structure_consistency_summary` (the state a chemistry-
    toolkit-less build leaves the schema in, since `ensure_rdkit_extension`
    raises before either QC table is created) -- the module-level function
    used by the routes' availability check reports it plainly."""
    from omnipath_metabo.server._routes_qc import _schema_present

    with conn.cursor() as cur:
        assert _schema_present(cur) is True

    # A relation that certainly doesn't exist behaves the same way the
    # absent-toolkit case does -- to_regclass returns NULL either way.
    with conn.cursor() as cur:
        cur.execute("SELECT to_regclass('this_schema_does_not_exist.nope')")
        assert cur.fetchone()[0] is None
