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

"""Tests for the Goslin lipid label cascade (US8 T066).

The ``parse_lipid`` unit tests need only ``pygoslin``. The integration test
additionally needs a built main DB (chemicals with names) reachable via
DATABASE_URL / OMNIPATH_DB_URL::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run --with psycopg2-binary --with pygoslin --with pytest \
        pytest tests/test_lipid_labels.py -v
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip('pygoslin')

from omnipath_metabo.postbuild._lipid_layer import (  # noqa: E402
    parse_lipid,
    resolve_lipid_labels,
)

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')
CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'


# --- parse_lipid unit tests (no DB) -------------------------------------

@pytest.mark.parametrize(
    'name',
    ['cholesterol', 'sphingosine', 'DHA', 'glucose', 'ATP', 'caffeine', ''],
)
def test_common_and_nonlipid_names_are_not_parsed(name):
    """Names without a C:D chain (or non-lipids) yield None — keep their name."""
    assert parse_lipid(name) is None


def test_molecular_species_keeps_granularity():
    """An sn / molecular name is kept at its own (highest) level, not summed."""
    sn = parse_lipid('PC 16:0/18:1')
    assert sn is not None
    assert sn['normalised_name'] == 'PC 16:0/18:1'
    assert sn['lipid_level'] == 'sn_position'

    molecular = parse_lipid('PC 16:0_18:1')
    assert molecular['normalised_name'] == 'PC 16:0_18:1'
    assert molecular['lipid_level'] == 'molecular_species'


def test_species_name_normalised():
    res = parse_lipid('PC(36:2)')
    assert res['normalised_name'] == 'PC 36:2'
    assert res['lipid_level'] == 'species'
    assert res['lipid_class'] == 'PC'
    assert res['lipid_category'] == 'GP'
    assert (res['total_carbon'], res['total_db']) == (36, 2)


def test_partial_tg_keeps_its_own_partially_specified_identity():
    """spec 011 T121/research R9: ``TG 42:3-FA18:2`` (a species total with one
    named FA, subtracted: 42:3 - 18:2 = 24:1) does NOT collapse to the
    species total ``TG 42:3`` -- the old guard's premise (that the
    subtracted rendering misrepresents a 3-acyl TG) does not hold: the
    listed parts (24:1 + 18:2) genuinely add back to the species total, and
    chains_possible/chains_listed (from the parse object, not string
    counting) correctly show 2 of TG's 3 chains listed, not 3."""
    res = parse_lipid('TG 42:3-FA18:2')
    assert res['normalised_name'] == 'TG 24:1_18:2'
    assert res['lipid_level'] == 'partially_specified'
    assert res['chains_possible'] == 3
    assert res['chains_listed'] == 2


def test_fully_specified_tg_keeps_all_chains():
    res = parse_lipid('TG 16:0/18:1/18:2')
    assert res['normalised_name'] == 'TG 16:0/18:1/18:2'
    assert res['lipid_level'] == 'sn_position'
    assert res['chains_possible'] == res['chains_listed'] == 3


@pytest.mark.parametrize(
    ('name', 'expected_chains'),
    [
        ('MG 18:1', 1),  # one-chain class: naive text-token count would be 3
        ('CE 18:1', 1),  # sterol ester: naive text-token count would be 2
    ],
)
def test_chain_counts_come_from_the_parse_not_rendered_tokens(name, expected_chains):
    """spec 011 T121: the two documented naive-token-counting bugs -- an
    empty-position-slot rendering and a sterol ester's backbone token --
    both yield the correct chain count from the parse object."""
    res = parse_lipid(name)
    assert res is not None
    assert res['chains_possible'] == expected_chains
    assert res['chains_listed'] == expected_chains


# --- DB integration -----------------------------------------------------

dbmark = pytest.mark.skipif(
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


def _scalar(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params or [])
        row = cur.fetchone()
        return row[0] if row else None


@dbmark
def test_resolve_lipid_labels_against_db(conn):
    """End-to-end: lipids get goslin_lipid labels; common names are untouched."""
    # The lipid cache table is created by ensure_metabo_schema in the real
    # pipeline; create it here so the test is independent of the rdkit cartridge.
    with conn.cursor() as cur:
        cur.execute(
            f"""
            CREATE TABLE IF NOT EXISTS {SCHEMA}.metabo_lipid_name_resolution (
              raw_name text PRIMARY KEY, normalised_name text, species_name text,
              lipid_level text, lipid_category text, lipid_class text,
              total_carbon smallint, total_db smallint, sum_formula text,
              resolver text NOT NULL, goslin_version text,
              chains_possible smallint, chains_listed smallint,
              computed_at timestamptz NOT NULL DEFAULT now())
            """
        )
    conn.commit()

    stats = resolve_lipid_labels(conn, schema=SCHEMA)
    if stats.lipids_labelled == 0:
        pytest.skip('no lipid-shorthand chemicals in this (capped) build')

    # Labelled lipids all carry the goslin_lipid rule and a non-empty label.
    bad = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE label_rule = 'goslin_lipid' AND (label IS NULL OR label = '')
        """,
    )
    assert bad == 0

    # Common lipid names are NOT turned into Goslin codes (chain guard).
    miscoded = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE lower(label) IN ('cholesterol', 'sphingosine')
          AND label_rule = 'goslin_lipid'
        """,
    )
    assert miscoded == 0

    # spec 011 T121/R9: a 2-of-3-chain TG label is no longer a bug signature
    # (the removed species-downgrade guard used to prevent it) -- it is the
    # correct rendering of a genuine partial specification. Confirm the cache
    # records it honestly: chains_listed < chains_possible, not silently
    # collapsed to a species total (chains_listed would then be 0).
    partial_tg = _scalar(
        conn,
        rf"""
        SELECT count(*) FROM {SCHEMA}.metabo_lipid_name_resolution
        WHERE normalised_name ~ '^TG [0-9]+:[0-9]+[_/][0-9]+:[0-9]+$'
          AND lipid_level = 'partially_specified'
        """,
    )
    if partial_tg:
        bad_counts = _scalar(
            conn,
            rf"""
            SELECT count(*) FROM {SCHEMA}.metabo_lipid_name_resolution
            WHERE normalised_name ~ '^TG [0-9]+:[0-9]+[_/][0-9]+:[0-9]+$'
              AND lipid_level = 'partially_specified'
              AND (chains_listed IS DISTINCT FROM 2
                OR chains_possible IS DISTINCT FROM 3)
            """,
        )
        assert bad_counts == 0
