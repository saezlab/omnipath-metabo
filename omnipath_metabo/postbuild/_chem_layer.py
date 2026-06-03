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

from dataclasses import dataclass, field

from psycopg2 import sql

from omnipath_metabo.db import STRUCTURAL_SPECIFICITY_LEVELS

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


def build_structure_substrate(conn, *, schema: str = 'public') -> int:
    """(Re)build ``metabo_entity_structure`` from parseable SMILES; return rows."""
    schema_id = sql.Identifier(schema)
    params = {
        'chem': CHEMICAL_ENTITY_TYPE,
        'smiles': SMILES_TYPE,
        'inchikey': INCHIKEY_TYPE,
    }
    with conn.cursor() as cur:
        # mol_from_smiles emits NOTICE/WARNING for unparseable inputs (it returns
        # NULL, not an error) — quiet those for the bulk build.
        cur.execute('SET LOCAL client_min_messages = error')
        cur.execute(
            sql.SQL('TRUNCATE {}.metabo_entity_structure').format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {schema}.metabo_entity_structure
                  (entity_id, mol, mfp, canonical_smiles, standard_inchikey)
                WITH chem AS (
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
                ),
                parsed AS (
                  SELECT entity_id, mol_from_smiles(smiles) AS mol, inchikey
                  FROM chem
                  WHERE smiles IS NOT NULL
                )
                SELECT
                  entity_id,
                  mol,
                  morganbv_fp(mol),
                  mol_to_smiles(mol),
                  inchikey
                FROM parsed
                WHERE mol IS NOT NULL
                """
            ).format(schema=schema_id),
            params,
        )
        rows = int(cur.rowcount)
    conn.commit()
    return rows


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
