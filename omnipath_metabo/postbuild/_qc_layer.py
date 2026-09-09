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
Structure consistency diagnostics (spec 011 data-model.md section 6, R8):
three checks, all precomputed here and served without request-time
computation (contracts/quality-control-api.md, Principle II).

The "structure key" here is canonical SMILES, not InChIKey -- this build's
rdkit Postgres cartridge is compiled without InChI support (``mol_inchikey``
returns the literal string ``InChI not available``; ``mol_from_inchi`` does
not exist at all), and `omnipath-metabo` deliberately runs chemistry *only*
through that cartridge, never Python rdkit (pyproject.toml's own comment on
the ``postbuild`` extra: "no Python RDKit" -- a pre-existing, intentional
choice, not something this cycle should override). Two SMILES renderings from
the SAME cartridge are safely comparable by string equality:

* the full canonical (isomeric) SMILES -- ``mol_to_smiles(mol)``
* the flat (non-isomeric) SMILES -- ``mol_to_smiles(mol, false)``, which
  drops stereo/cis-trans/isotope markers, the same "skeleton" grouping WP2's
  candidate arbitration uses for InChIKey's first 14 characters

**Narrowed scope (deferred-items.md):** the internal check as specified
compares a record's SMILES, InChI and InChIKey assertions against each other;
without InChI parsing that three-way comparison isn't possible here. What
this module builds instead: when the SAME (entity, authoritative source) pair
carries more than one *distinct raw SMILES string* (a real, if less common,
situation -- e.g. a corrected or duplicate assertion), do they canonicalize to
the same molecule. Records with exactly one asserted SMILES have nothing to
compare against and correctly produce no internal-check finding.

| Check | Applies to | Compares |
|---|---|---|
| internal | an entity with >1 distinct raw SMILES from one structure authority | those SMILES, pairwise, re-parsed and canonicalized |
| cross_reference | any entity citing an authority's identifier | its own structure substrate vs. the cited authority's |
| cross_reference_pair | an entity with no structure of its own, citing two authorities | the two cited authorities against each other |

Classification (only ``different_structure`` is an error, R8): ``agree``
(identical canonical SMILES), ``layer_difference`` (identical flat SMILES,
different full form -- stereo/isotope), ``unparsable`` (a SMILES the cartridge
could not parse), ``different_structure`` (different flat SMILES).
"""

from __future__ import annotations

__all__ = [
    'QcLayerStats',
    'build_structure_consistency_findings',
    'refresh_structure_consistency_summary',
]

from dataclasses import dataclass

from psycopg2 import sql

CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'
SMILES_TYPE = 'Smiles:MI:0239'


@dataclass(frozen=True)
class QcLayerStats:
    internal: int = 0
    cross_reference: int = 0
    cross_reference_pair: int = 0
    summary_rows: int = 0


def _verdict_case(full_a: str, full_b: str, flat_a: str, flat_b: str) -> sql.SQL:
    """The shared verdict CASE, parameterised on the (already-aliased) full and
    flat SMILES expressions of the two sides. Callers guarantee both full
    expressions are non-NULL by the time this runs -- unparsable findings are
    branched on separately, before this ever executes.
    """

    return sql.SQL(
        """
        CASE
          WHEN {fa} = {fb} THEN 'agree'
          WHEN {la} = {lb} THEN 'layer_difference'
          ELSE 'different_structure'
        END
        """
    ).format(
        fa=sql.SQL(full_a), fb=sql.SQL(full_b),
        la=sql.SQL(flat_a), lb=sql.SQL(flat_b),
    )


def _layer_case(full_a: str, full_b: str, flat_a: str, flat_b: str) -> sql.SQL:
    return sql.SQL(
        """
        CASE
          WHEN {fa} = {fb} THEN NULL
          WHEN {la} = {lb} THEN 'stereo_isotope'
          ELSE NULL
        END
        """
    ).format(
        fa=sql.SQL(full_a), fb=sql.SQL(full_b),
        la=sql.SQL(flat_a), lb=sql.SQL(flat_b),
    )


def _type_id(cur, schema_id: sql.Identifier, name: str) -> int | None:
    cur.execute(
        sql.SQL(
            'SELECT identifier_type_id FROM {}.vocab_identifier_type WHERE name = %s'
        ).format(schema_id),
        [name],
    )
    row = cur.fetchone()
    return int(row[0]) if row else None


def _entity_type_id(cur, schema_id: sql.Identifier, name: str) -> int | None:
    cur.execute(
        sql.SQL(
            'SELECT entity_type_id FROM {}.vocab_entity_type WHERE name = %s'
        ).format(schema_id),
        [name],
    )
    row = cur.fetchone()
    return int(row[0]) if row else None


def build_structure_consistency_findings(
    conn, *, schema: str = 'public',
) -> QcLayerStats:
    """(Re)build ``structure_consistency_finding`` -- all three checks, one
    full rebuild per call (idempotent). Requires the structure substrate
    (:func:`_chem_layer.build_structure_substrate`) to already be current for
    this build.
    """

    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        chem_type_id = _entity_type_id(cur, schema_id, CHEMICAL_ENTITY_TYPE)
        smiles_type_id = _type_id(cur, schema_id, SMILES_TYPE)
        if chem_type_id is None or smiles_type_id is None:
            cur.execute(
                sql.SQL('TRUNCATE {}.structure_consistency_finding').format(
                    schema_id
                )
            )
            conn.commit()
            return QcLayerStats()

        cur.execute(
            sql.SQL('TRUNCATE {}.structure_consistency_finding').format(schema_id)
        )
        params = dict(chem=chem_type_id, smiles=smiles_type_id)

        internal = _build_internal(cur, schema_id, params)
        _build_authority_map(cur, schema_id, params)
        cross_reference = _build_cross_reference(cur, schema_id, params)
        cross_reference_pair = _build_cross_reference_pair(cur, schema_id, params)
        cur.execute(
            sql.SQL('DROP TABLE IF EXISTS {}').format(
                sql.Identifier(schema, 'qc_authority_identifier_entity')
            )
        )
    conn.commit()
    summary_rows = refresh_structure_consistency_summary(conn, schema=schema)
    return QcLayerStats(
        internal=internal,
        cross_reference=cross_reference,
        cross_reference_pair=cross_reference_pair,
        summary_rows=summary_rows,
    )


def _build_internal(cur, schema_id: sql.Identifier, params: dict) -> int:
    """Do a record's own (possibly several) raw SMILES assertions describe one
    molecule? Scoped to entities a structure authority resolved, where that
    authority asserted more than one *distinct* raw SMILES string for the
    same entity -- the only self-consistency question the cartridge can
    actually answer without InChI (see module docstring).
    """

    cur.execute(
        sql.SQL(
            """
            WITH raw AS (
              SELECT DISTINCT
                eer.entity_id, eer.source_id, ie.value AS raw_smiles
              FROM {schema}.entity_evidence ee
              JOIN {schema}.entity_evidence_resolution eer
                ON eer.source_id = ee.source_id
               AND eer.entity_evidence_id = ee.entity_evidence_id
              JOIN {schema}.entity_evidence_identifier eei
                ON eei.source_id = ee.source_id
               AND eei.entity_evidence_id = ee.entity_evidence_id
              JOIN {schema}.identifier_evidence ie
                ON ie.identifier_id = eei.identifier_id
               AND ie.identifier_type_id = %(smiles)s
              JOIN {schema}.identifier_authority ia
                ON ia.source_id = eer.source_id AND ia.is_structure_authority
              WHERE ee.entity_type_id = %(chem)s
                AND eer.entity_id IS NOT NULL
            ),
            multi AS (
              SELECT entity_id, source_id
              FROM raw
              GROUP BY entity_id, source_id
              HAVING count(*) > 1
            ),
            parsed AS (
              SELECT r.entity_id, r.source_id, r.raw_smiles,
                     mol_from_smiles(r.raw_smiles::cstring) AS mol
              FROM raw r
              JOIN multi m ON m.entity_id = r.entity_id AND m.source_id = r.source_id
            )
            INSERT INTO {schema}.structure_consistency_finding (
              check_kind, source_id, identifier_type_id, value_normalized,
              authority_source_id, structure_a, structure_b, verdict, layer
            )
            SELECT
              'internal', a.source_id, %(smiles)s, a.entity_id::text,
              NULL,
              mol_to_smiles(a.mol)::text, mol_to_smiles(b.mol)::text,
              CASE
                WHEN a.mol IS NULL OR b.mol IS NULL THEN 'unparsable'
                ELSE ({verdict})
              END,
              CASE
                WHEN a.mol IS NULL OR b.mol IS NULL THEN NULL
                ELSE ({layer})
              END
            FROM parsed a
            JOIN parsed b
              ON b.entity_id = a.entity_id AND b.source_id = a.source_id
             AND b.raw_smiles > a.raw_smiles
            """
        ).format(
            schema=schema_id,
            verdict=_verdict_case(
                'mol_to_smiles(a.mol)::text', 'mol_to_smiles(b.mol)::text',
                'mol_to_smiles(a.mol, false)::text',
                'mol_to_smiles(b.mol, false)::text',
            ),
            layer=_layer_case(
                'mol_to_smiles(a.mol)::text', 'mol_to_smiles(b.mol)::text',
                'mol_to_smiles(a.mol, false)::text',
                'mol_to_smiles(b.mol, false)::text',
            ),
        ),
        params,
    )
    return int(cur.rowcount)


def _build_authority_map(cur, schema_id: sql.Identifier, params: dict) -> None:
    """Materialise, once, the mapping every structure-authority identifier
    value (on the authority's OWN records) resolves to -- reused by both the
    cross-reference and cross-reference-pair checks below.

    Driven from ``identifier_evidence`` filtered to the six structure-
    authority identifier types first (~2.9M rows), not from
    ``entity_evidence`` (which would scan every chemical mention of every
    kind before filtering down to structure-bearing ones). ``ee.source_id =
    ia.source_id`` restricts the ``entity_evidence`` fan-out to only the
    authorities' own (few) partitions.
    """

    cur.execute(
        sql.SQL('DROP TABLE IF EXISTS {}').format(
            sql.Identifier(schema_id.strings[-1], 'qc_authority_identifier_entity')
        )
    )
    cur.execute(
        sql.SQL(
            """
            CREATE UNLOGGED TABLE {schema}.qc_authority_identifier_entity AS
            SELECT DISTINCT
              ia.source_id AS authority_source_id, ie.identifier_type_id,
              COALESCE(ie.value_normalized, ie.value) AS value_normalized,
              eer.entity_id
            FROM {schema}.identifier_evidence ie
            JOIN {schema}.identifier_authority ia
              ON ia.identifier_type_id = ie.identifier_type_id
             AND ia.is_structure_authority
            JOIN {schema}.entity_evidence_identifier eei
              ON eei.identifier_id = ie.identifier_id
            JOIN {schema}.entity_evidence ee
              ON ee.source_id = eei.source_id
             AND ee.entity_evidence_id = eei.entity_evidence_id
             AND ee.source_id = ia.source_id
             AND ee.entity_type_id = %(chem)s
            JOIN {schema}.entity_evidence_resolution eer
              ON eer.source_id = ee.source_id
             AND eer.entity_evidence_id = ee.entity_evidence_id
            """
        ).format(schema=schema_id),
        params,
    )
    cur.execute(
        sql.SQL(
            'CREATE INDEX ON {}.qc_authority_identifier_entity '
            '(identifier_type_id, value_normalized)'
        ).format(schema_id)
    )


def _build_cross_reference(cur, schema_id: sql.Identifier, params: dict) -> int:
    """Does a record citing an authority's identifier agree with the
    authority's own (cartridge-canonicalized) structure? Scoped to
    entity_evidence mentions citing a structure-authority namespace, whichever
    entity the citing mention itself resolved to vs. the entity the cited
    identifier itself canonicalizes to (:func:`_build_authority_map`).
    """

    cur.execute(
        sql.SQL(
            """
            WITH citation AS (
              SELECT DISTINCT
                ee.source_id AS citing_source_id,
                eer.entity_id AS citing_entity_id,
                ia.identifier_type_id,
                COALESCE(ie.value_normalized, ie.value) AS value_normalized,
                ia.source_id AS authority_source_id
              FROM {schema}.identifier_evidence ie
              JOIN {schema}.identifier_authority ia
                ON ia.identifier_type_id = ie.identifier_type_id
               AND ia.is_structure_authority
              JOIN {schema}.entity_evidence_identifier eei
                ON eei.identifier_id = ie.identifier_id
              JOIN {schema}.entity_evidence ee
                ON ee.source_id = eei.source_id
               AND ee.entity_evidence_id = eei.entity_evidence_id
               AND ee.entity_type_id = %(chem)s
               AND ee.source_id <> ia.source_id
              JOIN {schema}.entity_evidence_resolution eer
                ON eer.source_id = ee.source_id
               AND eer.entity_evidence_id = ee.entity_evidence_id
            ),
            authority_entity AS (
              -- the entity the cited identifier itself canonicalizes to, via
              -- the precomputed authority map.
              SELECT DISTINCT
                cit.citing_source_id, cit.citing_entity_id,
                cit.identifier_type_id, cit.value_normalized,
                cit.authority_source_id, am.entity_id AS authority_entity_id
              FROM citation cit
              JOIN {schema}.qc_authority_identifier_entity am
                ON am.authority_source_id = cit.authority_source_id
               AND am.identifier_type_id = cit.identifier_type_id
               AND am.value_normalized = cit.value_normalized
            )
            INSERT INTO {schema}.structure_consistency_finding (
              check_kind, source_id, identifier_type_id, value_normalized,
              authority_source_id, structure_a, structure_b, verdict, layer
            )
            SELECT DISTINCT
              'cross_reference', ae.citing_source_id, ae.identifier_type_id,
              ae.value_normalized, ae.authority_source_id,
              s_citing.canonical_smiles, s_authority.canonical_smiles,
              CASE
                WHEN s_citing.mol IS NULL OR s_authority.mol IS NULL
                  THEN 'unparsable'
                ELSE ({verdict})
              END,
              CASE
                WHEN s_citing.mol IS NULL OR s_authority.mol IS NULL THEN NULL
                ELSE ({layer})
              END
            FROM authority_entity ae
            LEFT JOIN {schema}.metabo_entity_structure s_citing
              ON s_citing.entity_id = ae.citing_entity_id
            LEFT JOIN {schema}.metabo_entity_structure s_authority
              ON s_authority.entity_id = ae.authority_entity_id
            WHERE ae.citing_entity_id <> ae.authority_entity_id
            """
        ).format(
            schema=schema_id,
            verdict=_verdict_case(
                's_citing.canonical_smiles', 's_authority.canonical_smiles',
                'mol_to_smiles(s_citing.mol, false)::text',
                'mol_to_smiles(s_authority.mol, false)::text',
            ),
            layer=_layer_case(
                's_citing.canonical_smiles', 's_authority.canonical_smiles',
                'mol_to_smiles(s_citing.mol, false)::text',
                'mol_to_smiles(s_authority.mol, false)::text',
            ),
        ),
        params,
    )
    return int(cur.rowcount)


def _build_cross_reference_pair(
    cur, schema_id: sql.Identifier, params: dict,
) -> int:
    """A record with no structure of its own, citing two authorities: do the
    two authorities agree with each other? Needs no structure from the
    resource being checked -- what makes it applicable to resources that
    have none (research R8, the check that found the KEGG defect).
    """

    cur.execute(
        sql.SQL(
            """
            WITH citation AS (
              SELECT DISTINCT
                ee.source_id AS citing_source_id,
                ee.entity_evidence_id,
                ia.identifier_type_id,
                COALESCE(ie.value_normalized, ie.value) AS value_normalized,
                ia.source_id AS authority_source_id
              FROM {schema}.identifier_evidence ie
              JOIN {schema}.identifier_authority ia
                ON ia.identifier_type_id = ie.identifier_type_id
               AND ia.is_structure_authority
              JOIN {schema}.entity_evidence_identifier eei
                ON eei.identifier_id = ie.identifier_id
              JOIN {schema}.entity_evidence ee
                ON ee.source_id = eei.source_id
               AND ee.entity_evidence_id = eei.entity_evidence_id
               AND ee.entity_type_id = %(chem)s
               AND ee.source_id <> ia.source_id
            ),
            -- only mentions with no structure of their own (SMILES
            -- identifier attached directly) -- the pair check exists
            -- precisely because these resources have none.
            no_own_structure AS (
              SELECT citing_source_id, entity_evidence_id
              FROM citation
              EXCEPT
              SELECT ee.source_id, ee.entity_evidence_id
              FROM {schema}.entity_evidence ee
              JOIN {schema}.entity_evidence_identifier eei
                ON eei.source_id = ee.source_id
               AND eei.entity_evidence_id = ee.entity_evidence_id
              JOIN {schema}.identifier_evidence ie
                ON ie.identifier_id = eei.identifier_id
               AND ie.identifier_type_id = %(smiles)s
            ),
            authority_entity AS (
              SELECT DISTINCT
                cit.citing_source_id, cit.entity_evidence_id,
                cit.identifier_type_id, cit.value_normalized,
                cit.authority_source_id, am.entity_id AS authority_entity_id
              FROM citation cit
              JOIN no_own_structure nos
                ON nos.citing_source_id = cit.citing_source_id
               AND nos.entity_evidence_id = cit.entity_evidence_id
              JOIN {schema}.qc_authority_identifier_entity am
                ON am.authority_source_id = cit.authority_source_id
               AND am.identifier_type_id = cit.identifier_type_id
               AND am.value_normalized = cit.value_normalized
            ),
            -- pair up two DISTINCT authorities cited by the same mention.
            pair AS (
              SELECT
                a.citing_source_id, a.entity_evidence_id,
                a.identifier_type_id AS type_a, a.value_normalized AS value_a,
                a.authority_source_id AS authority_a,
                a.authority_entity_id AS entity_a,
                b.identifier_type_id AS type_b, b.value_normalized AS value_b,
                b.authority_source_id AS authority_b,
                b.authority_entity_id AS entity_b
              FROM authority_entity a
              JOIN authority_entity b
                ON b.citing_source_id = a.citing_source_id
               AND b.entity_evidence_id = a.entity_evidence_id
               AND (b.authority_source_id, b.value_normalized)
                 > (a.authority_source_id, a.value_normalized)
              WHERE a.authority_entity_id <> b.authority_entity_id
            )
            INSERT INTO {schema}.structure_consistency_finding (
              check_kind, source_id, identifier_type_id, value_normalized,
              authority_source_id, structure_a, structure_b, verdict, layer
            )
            SELECT DISTINCT
              'cross_reference_pair', p.citing_source_id, p.type_a, p.value_a,
              p.authority_a, s_a.canonical_smiles, s_b.canonical_smiles,
              CASE
                WHEN s_a.mol IS NULL OR s_b.mol IS NULL THEN 'unparsable'
                ELSE ({verdict})
              END,
              CASE
                WHEN s_a.mol IS NULL OR s_b.mol IS NULL THEN NULL
                ELSE ({layer})
              END
            FROM pair p
            LEFT JOIN {schema}.metabo_entity_structure s_a
              ON s_a.entity_id = p.entity_a
            LEFT JOIN {schema}.metabo_entity_structure s_b
              ON s_b.entity_id = p.entity_b
            """
        ).format(
            schema=schema_id,
            verdict=_verdict_case(
                's_a.canonical_smiles', 's_b.canonical_smiles',
                'mol_to_smiles(s_a.mol, false)::text', 'mol_to_smiles(s_b.mol, false)::text',
            ),
            layer=_layer_case(
                's_a.canonical_smiles', 's_b.canonical_smiles',
                'mol_to_smiles(s_a.mol, false)::text', 'mol_to_smiles(s_b.mol, false)::text',
            ),
        ),
        params,
    )
    return int(cur.rowcount)


def refresh_structure_consistency_summary(conn, *, schema: str = 'public') -> int:
    """(Re)build ``structure_consistency_summary`` from the findings table --
    what the endpoints actually serve (no request-time scan of findings).
    """

    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL('TRUNCATE {}.structure_consistency_summary').format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {schema}.structure_consistency_summary (
                  check_kind, source_id, authority_source_id, verdict,
                  finding_count
                )
                SELECT check_kind, source_id, authority_source_id, verdict,
                       count(*)
                FROM {schema}.structure_consistency_finding
                GROUP BY 1, 2, 3, 4
                """
            ).format(schema=schema_id)
        )
        rows = int(cur.rowcount)
    conn.commit()
    return rows
