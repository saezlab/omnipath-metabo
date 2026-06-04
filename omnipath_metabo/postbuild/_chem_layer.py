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
The in-DB structure substrate + structural-specificity classification (E).

All RDKit work happens in Postgres via the cartridge: ``mol_from_smiles`` builds
``metabo_entity_structure.mol`` (+ Morgan ``bfp`` + canonical SMILES) from the
build's ``Smiles:MI:0239`` identifiers; specificity is then derived in SQL from
the isomeric canonical SMILES (stereo markers) / InChIKey presence, covering
100% of ``Chemical:OM:0037`` entities (SC-006).
"""

from __future__ import annotations

__all__ = [
    'build_structure_substrate',
    'classify_structural_specificity',
    'refresh_structural_specificity_facet',
]

import os
import threading
from dataclasses import dataclass, field

from psycopg2 import sql

from omnipath_metabo.db import STRUCTURAL_SPECIFICITY_LEVELS, connect

# Build-phase Postgres session tuning for the RDKit parse + bulk GiST rebuilds.
# Overridable per deployment via these environment variables (shared with
# omnipath-build); set one to an empty string to leave that GUC at the server
# default. Defaults assume the build's Postgres has ~10 GB of headroom (the
# lab's docker.service hard cap is 200 GB shared across all containers). These
# are applied with SET LOCAL, so they last only for the build transaction and
# never affect the global config of an instance shared with the web API/app.
_MAX_PARALLEL_WORKERS_PER_GATHER = os.environ.get(
    'OMNIPATH_BUILD_MAX_PARALLEL_WORKERS_PER_GATHER', '8'
)
_MAINTENANCE_WORK_MEM = os.environ.get(
    'OMNIPATH_BUILD_MAINTENANCE_WORK_MEM', '2GB'
)
_MAX_PARALLEL_MAINTENANCE_WORKERS = os.environ.get(
    'OMNIPATH_BUILD_MAX_PARALLEL_MAINTENANCE_WORKERS', '4'
)

CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'
SMILES_TYPE = 'Smiles:MI:0239'
INCHIKEY_TYPE = 'Standard Inchi Key:MI:1101'
_FACET_NAME = 'structural_specificity'

# Stable ids (1-based, matching the CV seeding order).
_LEVEL_ID = {name: i for i, name in enumerate(STRUCTURAL_SPECIFICITY_LEVELS, 1)}


@dataclass(frozen=True)
class SpecificityStats:
    structures: int = 0
    chemicals: int = 0
    by_level: dict[str, int] = field(default_factory=dict)


def build_structure_substrate(
    conn, *, schema: str = 'public', db_url: str | None = None,
) -> int:
    """(Re)build ``metabo_entity_structure`` from parseable SMILES; return rows.

    Performance: the RDKit cartridge functions (``mol_from_smiles`` /
    ``morganbv_fp`` / ``mol_to_smiles``) are PARALLEL SAFE + IMMUTABLE, but
    ``INSERT ... SELECT`` forces a serial plan, so the (expensive) parse ran on
    a single core. This is a two-stage ``CREATE TABLE AS``:

    1. materialise the chemical SMILES source (the entity ↔ identifier join +
       group-by) into a plain table, keeping only rows that have a SMILES;
    2. parse it in a second ``CREATE TABLE AS`` whose **parallel sequential
       scan** distributes ``mol_from_smiles`` across workers — a ``LATERAL``
       evaluates the mol once per row (it is referenced three times). Doing the
       parse in one statement on top of the join's group-by put it *above* the
       Gather (serial on the leader); splitting the stages keeps it in the
       workers.

    The two GiST indexes are dropped first and rebuilt in bulk afterwards
    (per-row GiST maintenance during the load is far slower than one bulk
    build). The session gets generous parallel/maintenance memory for the
    duration.
    """
    schema_id = sql.Identifier(schema)
    chem_src_id = sql.Identifier(schema, 'metabo_entity_structure_chem_src')
    staging_id = sql.Identifier(schema, 'metabo_entity_structure_staging')
    params = {
        'chem': CHEMICAL_ENTITY_TYPE,
        'smiles': SMILES_TYPE,
        'inchikey': INCHIKEY_TYPE,
    }
    with conn.cursor() as cur:
        # mol_from_smiles emits NOTICE/WARNING for unparseable inputs (it returns
        # NULL, not an error) — quiet those for the bulk build.
        cur.execute('SET LOCAL client_min_messages = error')
        # Encourage a parallel plan for the RDKit parse and give the bulk GiST
        # rebuilds room to work in memory (configurable per deployment; see the
        # module-level OMNIPATH_BUILD_* defaults).
        cur.execute('SET LOCAL parallel_setup_cost = 0')
        cur.execute('SET LOCAL parallel_tuple_cost = 0')
        if _MAX_PARALLEL_WORKERS_PER_GATHER:
            cur.execute(
                'SET LOCAL max_parallel_workers_per_gather = %s',
                (_MAX_PARALLEL_WORKERS_PER_GATHER,),
            )
        if _MAINTENANCE_WORK_MEM:
            cur.execute(
                'SET LOCAL maintenance_work_mem = %s', (_MAINTENANCE_WORK_MEM,)
            )
        if _MAX_PARALLEL_MAINTENANCE_WORKERS:
            cur.execute(
                'SET LOCAL max_parallel_maintenance_workers = %s',
                (_MAX_PARALLEL_MAINTENANCE_WORKERS,),
            )

        # Drop the expensive GiST indexes; rebuilt in bulk after the load.
        cur.execute(
            sql.SQL(
                'DROP INDEX IF EXISTS {}.metabo_entity_structure_mol_idx'
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                'DROP INDEX IF EXISTS {}.metabo_entity_structure_mfp_idx'
            ).format(schema_id)
        )

        # Stage 1: materialise the SMILES source (join + group-by) into a plain
        # table, keeping only chemicals that have a SMILES. Done as its own
        # statement so the expensive parse in stage 2 is NOT stacked on top of
        # this group-by (which would push it above the Gather = serial).
        cur.execute(sql.SQL('DROP TABLE IF EXISTS {}').format(chem_src_id))
        cur.execute(
            sql.SQL(
                """
                CREATE UNLOGGED TABLE {chem_src} AS
                SELECT entity_id, smiles, inchikey
                FROM (
                  SELECT
                    e.entity_id,
                    max(ie.value) FILTER (WHERE vit.name = %(smiles)s) AS smiles,
                    max(ie.value) FILTER (WHERE vit.name = %(inchikey)s) AS inchikey
                  FROM {schema}.entity e
                  JOIN {schema}.vocab_entity_type et
                    ON et.entity_type_id = e.entity_type_id
                   AND et.name = %(chem)s
                  LEFT JOIN {schema}.entity_identifier_lookup eil
                    ON eil.entity_id = e.entity_id
                  LEFT JOIN {schema}.identifier_evidence ie
                    ON ie.identifier_id = eil.identifier_id
                  LEFT JOIN {schema}.vocab_identifier_type vit
                    ON vit.identifier_type_id = ie.identifier_type_id
                   AND vit.name IN (%(smiles)s, %(inchikey)s)
                  GROUP BY e.entity_id
                ) g
                WHERE g.smiles IS NOT NULL
                """
            ).format(schema=schema_id, chem_src=chem_src_id),
            params,
        )
        # Stats so the planner sizes chem_src and picks a parallel seq scan.
        cur.execute(sql.SQL('ANALYZE {}').format(chem_src_id))

        # Stage 2: parse. A parallel sequential scan of chem_src distributes
        # mol_from_smiles across workers; the LATERAL evaluates the mol once per
        # row (referenced three times below).
        cur.execute(sql.SQL('DROP TABLE IF EXISTS {}').format(staging_id))
        cur.execute(
            sql.SQL(
                """
                CREATE UNLOGGED TABLE {staging} AS
                SELECT
                  chem.entity_id,
                  m.mol AS mol,
                  morganbv_fp(m.mol) AS mfp,
                  -- mol_to_smiles returns cstring (a pseudo-type); cast to text
                  -- so CREATE TABLE AS can give the column a real type.
                  mol_to_smiles(m.mol)::text AS canonical_smiles,
                  chem.inchikey AS standard_inchikey
                FROM {chem_src} chem
                CROSS JOIN LATERAL (
                  SELECT mol_from_smiles(chem.smiles) AS mol
                ) m
                WHERE m.mol IS NOT NULL
                """
            ).format(staging=staging_id, chem_src=chem_src_id)
        )
        cur.execute(sql.SQL('DROP TABLE IF EXISTS {}').format(chem_src_id))

        # Move the parsed rows into the (now index-free) target, then rebuild
        # the GiST indexes in one bulk pass.
        cur.execute(
            sql.SQL('TRUNCATE {}.metabo_entity_structure').format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {schema}.metabo_entity_structure
                  (entity_id, mol, mfp, canonical_smiles, standard_inchikey)
                SELECT entity_id, mol, mfp, canonical_smiles, standard_inchikey
                FROM {staging}
                """
            ).format(schema=schema_id, staging=staging_id)
        )
        rows = int(cur.rowcount)
        cur.execute(sql.SQL('DROP TABLE IF EXISTS {}').format(staging_id))
    # Commit the parsed rows + dropped indexes so the concurrent index builders
    # (separate connections) can see the data.
    conn.commit()

    _rebuild_structure_gist_indexes(conn, schema=schema, db_url=db_url)
    return rows


def _rebuild_structure_gist_indexes(
    conn, *, schema: str, db_url: str | None,
) -> None:
    """Rebuild the two ``metabo_entity_structure`` GiST indexes.

    Postgres has no parallel GiST index build, so a single ``CREATE INDEX`` is
    stuck on one core. When a ``db_url`` is available the two indexes are built
    CONCURRENTLY on their own connections — two ``CREATE INDEX`` on the same
    table take compatible ShareLocks, so this uses two cores instead of one.
    Falls back to a sequential build on the passed connection otherwise.
    """
    schema_id = sql.Identifier(schema)
    specs = (
        ('metabo_entity_structure_mol_idx', 'mol'),
        ('metabo_entity_structure_mfp_idx', 'mfp'),
    )

    def _stmt(index_name: str, column: str):
        return sql.SQL(
            'CREATE INDEX {idx} ON {schema}.metabo_entity_structure '
            'USING gist ({col})'
        ).format(
            idx=sql.Identifier(index_name),
            schema=schema_id,
            col=sql.Identifier(column),
        )

    if not db_url:
        with conn.cursor() as cur:
            for index_name, column in specs:
                cur.execute(_stmt(index_name, column))
        conn.commit()
        return

    errors: list[Exception] = []

    def _worker(index_name: str, column: str) -> None:
        try:
            index_conn = connect(db_url)
            index_conn.autocommit = True
            try:
                with index_conn.cursor() as cur:
                    if _MAINTENANCE_WORK_MEM:
                        cur.execute(
                            'SET maintenance_work_mem = %s',
                            (_MAINTENANCE_WORK_MEM,),
                        )
                    cur.execute(_stmt(index_name, column))
            finally:
                index_conn.close()
        except Exception as exc:  # noqa: BLE001 — surfaced after join
            errors.append(exc)

    threads = [
        threading.Thread(target=_worker, args=spec) for spec in specs
    ]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if errors:
        raise errors[0]


def classify_structural_specificity(
    conn,
    *,
    schema: str = 'public',
) -> SpecificityStats:
    """Assign one specificity level to every chemical (idempotent full rebuild)."""
    schema_id = sql.Identifier(schema)
    params = {
        'chem': CHEMICAL_ENTITY_TYPE,
        'inchikey': INCHIKEY_TYPE,
        'variable': _LEVEL_ID['variable_constitution'],
        'stereo': _LEVEL_ID['stereospecific'],
        'cistrans': _LEVEL_ID['cis_trans_only'],
        'constitution': _LEVEL_ID['constitution_only'],
        'unknown': _LEVEL_ID['unknown_constitution'],
        'none': _LEVEL_ID['no_structure'],
    }
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                'TRUNCATE {}.metabo_entity_structural_specificity'
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {schema}.metabo_entity_structural_specificity
                  (entity_id, structural_specificity_id, standard_inchikey)
                WITH chem AS (
                  SELECT
                    e.entity_id,
                    bool_or(vit.name = %(inchikey)s) AS has_key,
                    max(ie.value) FILTER (WHERE vit.name = %(inchikey)s) AS inchikey
                  FROM {schema}.entity e
                  JOIN {schema}.vocab_entity_type et
                    ON et.entity_type_id = e.entity_type_id
                   AND et.name = %(chem)s
                  LEFT JOIN {schema}.entity_identifier_lookup eil
                    ON eil.entity_id = e.entity_id
                  LEFT JOIN {schema}.identifier_evidence ie
                    ON ie.identifier_id = eil.identifier_id
                  LEFT JOIN {schema}.vocab_identifier_type vit
                    ON vit.identifier_type_id = ie.identifier_type_id
                   AND vit.name = %(inchikey)s
                  GROUP BY e.entity_id
                )
                SELECT
                  chem.entity_id,
                  CASE
                    WHEN strpos(s.canonical_smiles, '*') > 0 THEN %(variable)s
                    WHEN strpos(s.canonical_smiles, '@') > 0 THEN %(stereo)s
                    WHEN strpos(s.canonical_smiles, '/') > 0
                      OR strpos(s.canonical_smiles, chr(92)) > 0 THEN %(cistrans)s
                    WHEN s.canonical_smiles IS NOT NULL THEN %(constitution)s
                    WHEN chem.has_key THEN %(unknown)s
                    ELSE %(none)s
                  END,
                  COALESCE(s.standard_inchikey, chem.inchikey)
                FROM chem
                LEFT JOIN {schema}.metabo_entity_structure s
                  ON s.entity_id = chem.entity_id
                """
            ).format(schema=schema_id),
            params,
        )
        chemicals = int(cur.rowcount)
        cur.execute(
            sql.SQL(
                """
                SELECT v.name, count(*)
                FROM {schema}.metabo_entity_structural_specificity s
                JOIN {schema}.metabo_vocab_structural_specificity v
                  ON v.structural_specificity_id = s.structural_specificity_id
                GROUP BY v.name
                """
            ).format(schema=schema_id)
        )
        by_level = {name: int(count) for name, count in cur.fetchall()}
        cur.execute(
            sql.SQL(
                'SELECT count(*) FROM {}.metabo_entity_structure'
            ).format(schema_id)
        )
        structures = int(cur.fetchone()[0])
    conn.commit()
    return SpecificityStats(
        structures=structures,
        chemicals=chemicals,
        by_level=by_level,
    )


def refresh_structural_specificity_facet(conn, *, schema: str = 'public') -> int:
    """(Re)build the ``structural_specificity`` rows in ``facet_entity_bitmap``.

    Same roaringbitmap convention as omnipath-build's derive facets, so the web
    app surfaces it with no server change. Additive after derive's facet rebuild.
    """
    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                'DELETE FROM {}.facet_entity_bitmap WHERE facet_name = %s'
            ).format(schema_id),
            [_FACET_NAME],
        )
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {schema}.facet_entity_bitmap
                  (facet_name, facet_value, entity_bitmap, entity_count)
                SELECT
                  %s,
                  v.name,
                  rb_build_agg(b.bitmap_id),
                  count(*)::integer
                FROM {schema}.metabo_entity_structural_specificity s
                JOIN {schema}.metabo_vocab_structural_specificity v
                  ON v.structural_specificity_id = s.structural_specificity_id
                JOIN {schema}.entity_bitmap_id b
                  ON b.entity_id = s.entity_id
                GROUP BY v.name
                """
            ).format(schema=schema_id),
            [_FACET_NAME],
        )
        facet_rows = int(cur.rowcount)
    conn.commit()
    return facet_rows
