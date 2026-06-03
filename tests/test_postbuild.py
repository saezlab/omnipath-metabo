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

"""Integration tests for the metabo post-build chemistry layer (Milestones E/F).

Run against a built main DB whose Postgres has the rdkit cartridge, e.g.::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with psycopg2-binary --with pytest pytest tests/test_postbuild.py -v

Skipped when DATABASE_URL / OMNIPATH_DB_URL is unset.
"""

from __future__ import annotations

import os

import pytest

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')
CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'

pytestmark = pytest.mark.skipif(
    not DB_URL, reason='No DATABASE_URL/OMNIPATH_DB_URL; needs a built main DB'
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DB_URL)
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture(scope='module')
def post_build(conn):
    """Run the E layer (substrate + specificity + facet) once for the module."""
    from omnipath_metabo.postbuild import post_build_metabo

    return post_build_metabo(conn, schema=SCHEMA, conflicts=False, log=lambda *_: None)


def _scalar(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchone()[0]


def test_structure_substrate_has_molecules(conn, post_build):
    assert post_build.structures > 0
    assert _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.metabo_entity_structure') > 0
    # mol/bfp columns are populated by the cartridge.
    assert (
        _scalar(
            conn,
            f'SELECT count(*) FROM {SCHEMA}.metabo_entity_structure '
            f'WHERE mol IS NOT NULL AND mfp IS NOT NULL',
        )
        == post_build.structures
    )


def test_every_chemical_has_a_specificity_level(conn, post_build):
    """SC-006: 100% of Chemical:OM:0037 entities get a level."""
    chemicals = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity e
        JOIN {SCHEMA}.vocab_entity_type v ON v.entity_type_id = e.entity_type_id
        WHERE v.name = %s
        """,
        [CHEMICAL_ENTITY_TYPE],
    )
    classified = _scalar(
        conn,
        f'SELECT count(*) FROM {SCHEMA}.metabo_entity_structural_specificity',
    )
    assert classified == chemicals


def test_multiple_specificity_levels_present(conn, post_build):
    n_levels = _scalar(
        conn,
        f'SELECT count(DISTINCT structural_specificity_id) '
        f'FROM {SCHEMA}.metabo_entity_structural_specificity',
    )
    assert n_levels >= 3


def test_specificity_cv_seeded(conn, post_build):
    names = {
        row
        for (row,) in _rows(
            conn, f'SELECT name FROM {SCHEMA}.metabo_vocab_structural_specificity'
        )
    }
    assert names == {
        'stereospecific',
        'cis_trans_only',
        'constitution_only',
        'variable_constitution',
        'unknown_constitution',
        'no_structure',
    }


def test_structural_specificity_facet_present(conn, post_build):
    assert (
        _scalar(
            conn,
            f"SELECT count(*) FROM {SCHEMA}.facet_entity_bitmap "
            f"WHERE facet_name = 'structural_specificity'",
        )
        > 0
    )


def test_post_build_refuses_on_build_id_mismatch(conn, post_build):
    """FR-018: refuse to write when the recorded build id no longer matches."""
    from omnipath_metabo.db import record_build_id, current_build_id
    from omnipath_metabo.postbuild import post_build_metabo

    real = current_build_id(conn, SCHEMA)
    record_build_id(conn, SCHEMA, 'deadbeefdead')  # simulate a stale prior build
    try:
        with pytest.raises(RuntimeError, match='force'):
            post_build_metabo(conn, schema=SCHEMA, conflicts=False, log=lambda *_: None)
        # force=True recomputes and re-points to the real build.
        post_build_metabo(
            conn, schema=SCHEMA, force=True, conflicts=False, log=lambda *_: None
        )
        from omnipath_metabo.db import recorded_build_id

        assert recorded_build_id(conn, SCHEMA) == real
    finally:
        record_build_id(conn, SCHEMA, real)


@pytest.mark.slow
def test_ramp_conflicts_populated(conn):
    """Milestone F: RaMP multi-InChIKey conflicts get classified (needs RaMP)."""
    from omnipath_metabo.postbuild._ramp_conflicts import (
        CONFLICT_REASONS,
        populate_ramp_conflicts,
    )

    try:
        stats = populate_ramp_conflicts(conn, schema=SCHEMA, max_records=20000)
    except Exception as exc:  # RaMP download/cache unavailable in this environment
        pytest.skip(f'RaMP unavailable: {exc}')
    assert stats.ramp_rows > 0
    reasons = {
        row
        for (row,) in _rows(
            conn,
            f'SELECT DISTINCT conflict_reason '
            f'FROM {SCHEMA}.metabo_ramp_inchikey_conflict',
        )
    }
    assert reasons <= set(CONFLICT_REASONS)


def _rows(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()
