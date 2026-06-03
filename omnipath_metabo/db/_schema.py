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

"""
``metabo_*`` schema for the in-DB chemistry layer + the build-id gate.

The chemistry runs in Postgres through the RDKit cartridge (``CREATE EXTENSION
rdkit``); :func:`ensure_metabo_schema` ensures the extension is present (and
errors clearly if not), the foundational ``metabo_entity_structure`` substrate
(``mol``/``bfp`` columns, GiST-indexed for substructure + similarity), and the
classification / conflict / build-state tables. All ``CREATE … IF NOT EXISTS`` so
the post-build is self-contained.
"""

from __future__ import annotations

__all__ = [
    'STRUCTURAL_SPECIFICITY_LEVELS',
    'ensure_rdkit_extension',
    'ensure_metabo_schema',
    'current_build_id',
    'recorded_build_id',
    'record_build_id',
]

from psycopg2 import sql
import psycopg2

# Ordered most → least specific; index (1-based) is the stored
# ``structural_specificity_id`` and the CV seeding order.
STRUCTURAL_SPECIFICITY_LEVELS: tuple[str, ...] = (
    'stereospecific',
    'cis_trans_only',
    'constitution_only',
    'variable_constitution',
    'unknown_constitution',
    'no_structure',
)


def ensure_rdkit_extension(conn) -> None:
    """``CREATE EXTENSION rdkit`` — fail with a clear message if unavailable."""
    try:
        with conn.cursor() as cur:
            cur.execute('CREATE EXTENSION IF NOT EXISTS rdkit')
        conn.commit()
    except psycopg2.Error as exc:
        conn.rollback()
        raise RuntimeError(
            'The RDKit Postgres cartridge is required for the metabo post-build '
            'but is not available (CREATE EXTENSION rdkit failed). Install the '
            'pgdg `postgresql-NN-rdkit` package into the Postgres image '
            '(omnipath-present/postgres/Dockerfile).'
        ) from exc


def ensure_metabo_schema(conn, *, schema: str = 'public') -> None:
    """Ensure the cartridge + the ``metabo_*`` tables + the specificity CV."""
    ensure_rdkit_extension(conn)
    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.metabo_vocab_structural_specificity (
                  structural_specificity_id smallint PRIMARY KEY,
                  name text NOT NULL UNIQUE
                )
                """
            ).format(schema_id)
        )
        cur.executemany(
            sql.SQL(
                """
                INSERT INTO {}.metabo_vocab_structural_specificity
                  (structural_specificity_id, name)
                VALUES (%s, %s)
                ON CONFLICT (structural_specificity_id) DO UPDATE
                SET name = EXCLUDED.name
                """
            ).format(schema_id).as_string(cur.connection),
            [
                (index, name)
                for index, name in enumerate(STRUCTURAL_SPECIFICITY_LEVELS, start=1)
            ],
        )
        # Foundational RDKit substrate: one molecule per chemical with a
        # parseable SMILES, GiST-indexed for substructure + similarity search.
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {schema}.metabo_entity_structure (
                  entity_id uuid PRIMARY KEY REFERENCES {schema}.entity(entity_id),
                  mol mol NOT NULL,
                  mfp bfp NOT NULL,
                  canonical_smiles text NOT NULL,
                  standard_inchikey text,
                  computed_at timestamptz NOT NULL DEFAULT now()
                )
                """
            ).format(schema=schema_id)
        )
        cur.execute(
            sql.SQL(
                'CREATE INDEX IF NOT EXISTS metabo_entity_structure_mol_idx '
                'ON {}.metabo_entity_structure USING gist (mol)'
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                'CREATE INDEX IF NOT EXISTS metabo_entity_structure_mfp_idx '
                'ON {}.metabo_entity_structure USING gist (mfp)'
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {schema}.metabo_entity_structural_specificity (
                  entity_id uuid PRIMARY KEY REFERENCES {schema}.entity(entity_id),
                  structural_specificity_id smallint NOT NULL
                    REFERENCES {schema}.metabo_vocab_structural_specificity(
                      structural_specificity_id
                    ),
                  standard_inchikey text,
                  computed_at timestamptz NOT NULL DEFAULT now()
                )
                """
            ).format(schema=schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.metabo_ramp_inchikey_conflict (
                  ramp_id text NOT NULL,
                  inchikey_a text NOT NULL,
                  inchikey_b text NOT NULL,
                  conflict_reason text NOT NULL,
                  PRIMARY KEY (ramp_id, inchikey_a, inchikey_b)
                )
                """
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.metabo_build_state (
                  build_id text PRIMARY KEY,
                  computed_at timestamptz NOT NULL DEFAULT now()
                )
                """
            ).format(schema_id)
        )
    conn.commit()


def current_build_id(conn, schema: str = 'public') -> str:
    """The main build's id from ``build_manifest`` (raises if the build has none)."""
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL('SELECT build_id FROM {}.build_manifest').format(
                sql.Identifier(schema)
            )
        )
        row = cur.fetchone()
    if row is None:
        raise RuntimeError(
            'No build_manifest row; run the main build derive step first '
            '(Milestone D) before the metabo post-build.'
        )
    return row[0]


def recorded_build_id(conn, schema: str = 'public') -> str | None:
    """The build id the metabo layer was last computed for, or ``None``."""
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL('SELECT build_id FROM {}.metabo_build_state').format(
                sql.Identifier(schema)
            )
        )
        row = cur.fetchone()
    return row[0] if row else None


def record_build_id(conn, schema: str, build_id: str) -> None:
    """Record (replacing) the build id the metabo layer was computed for."""
    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(sql.SQL('TRUNCATE {}.metabo_build_state').format(schema_id))
        cur.execute(
            sql.SQL(
                'INSERT INTO {}.metabo_build_state (build_id) VALUES (%s)'
            ).format(schema_id),
            [build_id],
        )
    conn.commit()
