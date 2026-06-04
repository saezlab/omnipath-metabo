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
RaMP multi-InChIKey conflict classifier (Milestone F).

RaMP rows are fetched (via pypath) into a staging table; the structural
comparison is **in-DB** with the RDKit cartridge. A ``ramp_id`` mapped to >1
distinct InChIKey is a conflict; each pair is labelled
``stereo | specificity | tautomer | similar | unrelated``.
"""

from __future__ import annotations

__all__ = ['CONFLICT_REASONS', 'populate_ramp_conflicts']

from dataclasses import dataclass

import os

from psycopg2 import sql
from psycopg2.extras import execute_values

CONFLICT_REASONS = ('stereo', 'specificity', 'tautomer', 'similar', 'unrelated')
_SIMILAR_THRESHOLD = 0.7
# Parallelism for the in-DB RDKit build (shared OMNIPATH_BUILD_* knob; see
# omnipath-build/docs/build-tuning.md).
_MAX_PARALLEL_WORKERS_PER_GATHER = os.environ.get(
    'OMNIPATH_BUILD_MAX_PARALLEL_WORKERS_PER_GATHER', '8'
)


@dataclass(frozen=True)
class RampConflictStats:
    ramp_rows: int = 0
    conflicts: int = 0


def _fetch_ramp_rows(max_records: int | None) -> list[tuple[str, str, str | None]]:
    """``(ramp_id, inchikey, iso_smiles)`` from RaMP, InChIKey-bearing rows only."""
    from pypath.inputs.ramp import ramp_omnipathmetabo

    rows: list[tuple[str, str, str | None]] = []
    for record in ramp_omnipathmetabo():
        inchikey = getattr(record, 'inchi_key', None)
        ramp_id = getattr(record, 'ramp_id', None)
        if not ramp_id or not inchikey:
            continue
        rows.append((ramp_id, inchikey, getattr(record, 'iso_smiles', None)))
        if max_records and len(rows) >= max_records:
            break
    return rows


def populate_ramp_conflicts(
    conn,
    *,
    schema: str = 'public',
    max_records: int | None = None,
) -> RampConflictStats:
    """Classify RaMP multi-InChIKey conflicts into ``metabo_ramp_inchikey_conflict``."""

    ramp_rows = _fetch_ramp_rows(max_records)
    schema_id = sql.Identifier(schema)
    stage_id = sql.Identifier(schema, '_metabo_ramp_stage')
    distinct_id = sql.Identifier(schema, '_metabo_ramp_distinct')
    mol_id = sql.Identifier(schema, '_metabo_ramp_mol')
    with conn.cursor() as cur:
        cur.execute('SET LOCAL client_min_messages = error')
        cur.execute('SET LOCAL parallel_setup_cost = 0')
        cur.execute('SET LOCAL parallel_tuple_cost = 0')
        if _MAX_PARALLEL_WORKERS_PER_GATHER:
            cur.execute(
                'SET LOCAL max_parallel_workers_per_gather = %s',
                (_MAX_PARALLEL_WORKERS_PER_GATHER,),
            )
        # Regular UNLOGGED tables, NOT TEMP: parallel workers cannot read a
        # session temp table, so a TEMP source would force the in-DB RDKit build
        # single-core.
        for table_id in (mol_id, distinct_id, stage_id):
            cur.execute(sql.SQL('DROP TABLE IF EXISTS {}').format(table_id))
        cur.execute(
            sql.SQL(
                'CREATE UNLOGGED TABLE {} '
                '(ramp_id text, inchikey text, smiles text)'
            ).format(stage_id)
        )
        if ramp_rows:
            execute_values(
                cur,
                sql.SQL(
                    'INSERT INTO {} (ramp_id, inchikey, smiles) VALUES %s'
                ).format(stage_id).as_string(cur),
                ramp_rows,
                page_size=5000,
            )
        # Dedup to one row per (ramp_id, inchikey) with a SMILES BEFORE the
        # parse (cheap), so the parse runs only once per distinct molecule.
        cur.execute(
            sql.SQL(
                """
                CREATE UNLOGGED TABLE {distinct} AS
                SELECT DISTINCT ON (ramp_id, inchikey)
                  ramp_id, inchikey, split_part(inchikey, '-', 1) AS skel, smiles
                FROM {stage}
                WHERE smiles IS NOT NULL
                ORDER BY ramp_id, inchikey
                """
            ).format(distinct=distinct_id, stage=stage_id)
        )
        cur.execute(sql.SQL('ANALYZE {}').format(distinct_id))
        # Build molecule + fingerprint + formula + stereo flag in one parallel
        # CREATE TABLE AS — the parallel seq scan distributes the RDKit work
        # across workers; the LATERAL parses each SMILES once.
        cur.execute(
            sql.SQL(
                """
                CREATE UNLOGGED TABLE {mol} AS
                SELECT
                  d.ramp_id,
                  d.inchikey,
                  d.skel,
                  m.mol,
                  morganbv_fp(m.mol) AS fp,
                  -- mol_formula returns cstring (pseudo-type); cast to text.
                  mol_formula(m.mol)::text AS formula,
                  (strpos(mol_to_smiles(m.mol)::text, '@') > 0
                   OR strpos(mol_to_smiles(m.mol)::text, '/') > 0
                   OR strpos(mol_to_smiles(m.mol)::text, chr(92)) > 0)
                    AS has_stereo
                FROM {distinct} d
                CROSS JOIN LATERAL (SELECT mol_from_smiles(d.smiles) AS mol) m
                """
            ).format(mol=mol_id, distinct=distinct_id)
        )
        # Conflicting pairs: same ramp_id, distinct InChIKey, classified in-DB.
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {schema}.metabo_ramp_inchikey_conflict
                  (ramp_id, inchikey_a, inchikey_b, conflict_reason)
                SELECT
                  a.ramp_id,
                  a.inchikey,
                  b.inchikey,
                  CASE
                    WHEN a.skel = b.skel THEN
                      CASE WHEN coalesce(a.has_stereo, false)
                                <> coalesce(b.has_stereo, false)
                           THEN 'specificity' ELSE 'stereo' END
                    WHEN a.mol IS NULL OR b.mol IS NULL THEN 'unrelated'
                    WHEN a.formula = b.formula
                         AND tanimoto_sml(a.fp, b.fp) >= %(threshold)s
                      THEN 'tautomer'
                    WHEN tanimoto_sml(a.fp, b.fp) >= %(threshold)s THEN 'similar'
                    ELSE 'unrelated'
                  END
                FROM {mol} a
                JOIN {mol} b
                  ON a.ramp_id = b.ramp_id
                 AND a.inchikey < b.inchikey
                ON CONFLICT (ramp_id, inchikey_a, inchikey_b) DO UPDATE
                SET conflict_reason = EXCLUDED.conflict_reason
                """
            ).format(schema=schema_id, mol=mol_id),
            {'threshold': _SIMILAR_THRESHOLD},
        )
        conflicts = int(cur.rowcount)
        for table_id in (mol_id, distinct_id, stage_id):
            cur.execute(sql.SQL('DROP TABLE IF EXISTS {}').format(table_id))
    conn.commit()
    return RampConflictStats(ramp_rows=len(ramp_rows), conflicts=conflicts)
