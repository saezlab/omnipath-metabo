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

from psycopg2 import sql
from psycopg2.extras import execute_values

CONFLICT_REASONS = ('stereo', 'specificity', 'tautomer', 'similar', 'unrelated')
_SIMILAR_THRESHOLD = 0.7


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
    with conn.cursor() as cur:
        cur.execute('SET LOCAL client_min_messages = error')
        cur.execute(
            'CREATE TEMP TABLE _ramp_stage '
            '(ramp_id text, inchikey text, smiles text) ON COMMIT DROP'
        )
        if ramp_rows:
            execute_values(
                cur,
                'INSERT INTO _ramp_stage (ramp_id, inchikey, smiles) VALUES %s',
                ramp_rows,
                page_size=5000,
            )
        # Per (ramp_id, inchikey): build the molecule + fingerprint + formula in-DB.
        cur.execute(
            """
            CREATE TEMP TABLE _ramp_mol ON COMMIT DROP AS
            SELECT DISTINCT ON (ramp_id, inchikey)
              ramp_id,
              inchikey,
              split_part(inchikey, '-', 1) AS skel,
              mol_from_smiles(smiles) AS mol
            FROM _ramp_stage
            WHERE smiles IS NOT NULL
            ORDER BY ramp_id, inchikey
            """
        )
        cur.execute(
            'ALTER TABLE _ramp_mol '
            'ADD COLUMN fp bfp, ADD COLUMN has_stereo boolean, '
            'ADD COLUMN formula text'
        )
        cur.execute(
            """
            UPDATE _ramp_mol SET
              fp = morganbv_fp(mol),
              formula = mol_formula(mol),
              has_stereo = (
                strpos(mol_to_smiles(mol)::text, '@') > 0
                OR strpos(mol_to_smiles(mol)::text, '/') > 0
                OR strpos(mol_to_smiles(mol)::text, chr(92)) > 0
              )
            WHERE mol IS NOT NULL
            """
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
                FROM _ramp_mol a
                JOIN _ramp_mol b
                  ON a.ramp_id = b.ramp_id
                 AND a.inchikey < b.inchikey
                ON CONFLICT (ramp_id, inchikey_a, inchikey_b) DO UPDATE
                SET conflict_reason = EXCLUDED.conflict_reason
                """
            ).format(schema=schema_id),
            {'threshold': _SIMILAR_THRESHOLD},
        )
        conflicts = int(cur.rowcount)
    conn.commit()
    return RampConflictStats(ramp_rows=len(ramp_rows), conflicts=conflicts)
