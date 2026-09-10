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
Goslin lipid-name labels (US8 T066, data-model D; Principle II — Goslin/RDKit
live in the metabo layer).

Goslin is a lipid **nomenclature** grammar (LIPID MAPS / SwissLipids / HMDB /
shorthand-2020 dialects); it parses lipid *names* — never SMILES — into a
normalised name + a structural ``level`` on the SwissLipids hierarchy
(category → class → species → molecular_species → sn_position →
structure_defined → full_structure → complete_structure).

Goslin can normalise a name but never invent detail a name does not carry, so
per lipid entity we feed every candidate name, keep the parses, and label the
entity with the **highest-granularity** parse available, rendered at *its own*
level (``PC 16:0/18:1`` wins over ``PC 34:1``; an E/Z-resolved name wins over a
plain sn name). We never aggregate a detailed name down to species.

A chain guard (the name must contain a ``C:D`` token) keeps this away from the
common names of lipids that Goslin would otherwise turn into a cryptic code —
e.g. ``cholesterol`` → ``ST 27:1;O`` or ``sphingosine`` → ``SPB 18:1;O2`` — so
those keep the human name the brevity-first chemical cascade (T064) chose.

Parses are cached in ``metabo_lipid_name_resolution`` keyed on the verbatim
source name, so rebuilds never re-parse.

**Enrichment only (spec 011 T122).** Lipid *identity* -- which name+level+
chain-count is one entity, and the generalization edges between them -- is
decided upstream, in the main build (``duckdb_load.py``'s
``lipid_name_evidence_resolution``, reading the same standardized-name
corpus this module's own Goslin parses ultimately feed via
``omnipath-utils``' independent ``mapping/_lipid.py``). This module never
merges or creates entities; it only sets ``entity.label``/``label_rule`` for
display (and, when the toolkit allows, could enrich with a generated
structure -- see T123's status below). That was already true before this
cycle (the function only ever ``UPDATE``s ``label``), so T122 needed no
identity-related code removed here -- only T121's chain-counting fix and
the docstring making the boundary explicit.
"""

from __future__ import annotations

__all__ = ['LipidLabelStats', 'parse_lipid', 'resolve_lipid_labels']

import re
import importlib.metadata
from dataclasses import dataclass, field

from psycopg2 import sql
from psycopg2.extras import execute_values

# Name identifier types that may carry a lipid shorthand.
_NAME_TYPES = [
    'Name:OM:0202',
    'Synonym:OM:0203',
    'Abbreviated Name:OM:0208',
    'Iupac Traditional Name:OM:0211',
]
CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'
GOSLIN_LIPID_RULE = 'goslin_lipid'

# A C:D chain token — present in every lipid shorthand, absent from common
# names (cholesterol, sphingosine, DHA, …). Mirrors the SQL prefilter.
_CHAIN_RE = re.compile(r'\d+:\d+')

# Goslin LipidLevel name (lowercased slug) → granularity rank (higher = finer).
# Matches pygoslin.domain.LipidLevel; used to pick the most-resolved candidate.
_LEVEL_RANK = {
    'complete_structure': 9,
    'full_structure': 8,
    'structure_defined': 7,
    'sn_position': 6,
    'molecular_species': 5,
    # Synthesised, not a pygoslin level (T121/R9): strictly between species
    # and molecular_species -- some chains individually listed, the rest
    # given only as a sum. Ranked accordingly, never confused with either.
    'partially_specified': 4,
    'species': 3,
    'class': 2,
    'category': 1,
    'undefined': 0,
}

_parser = None  # lazy singleton — LipidParser construction is not free


def _goslin_version() -> str:
    try:
        return importlib.metadata.version('pygoslin')
    except Exception:  # noqa: BLE001
        return ''


def _get_parser():
    global _parser
    if _parser is None:
        from pygoslin.parser.Parser import LipidParser

        _parser = LipidParser()
    return _parser


def parse_lipid(name: str) -> dict | None:
    """Parse one lipid name with Goslin → normalised fields, or ``None``.

    Returns ``None`` when the name has no ``C:D`` chain token or Goslin cannot
    parse it as a lipid (the caller records those as ``unresolved``).
    """
    if not name or not _CHAIN_RE.search(name):
        return None
    from pygoslin.domain.LipidLevel import LipidLevel

    try:
        parsed = _get_parser().parse(name)
        normalised = parsed.get_lipid_string()
    except Exception:  # noqa: BLE001 — Goslin raises a family of parse errors
        return None

    info = getattr(parsed.lipid, 'info', None)
    level = getattr(info, 'level', None)
    level_slug = level.name.lower() if level is not None else None

    def _try(fn):
        try:
            return fn()
        except Exception:  # noqa: BLE001
            return None

    species = _try(lambda: parsed.get_lipid_string(LipidLevel.SPECIES))
    category = _try(lambda: parsed.lipid.headgroup.lipid_category.name)
    lipid_class = _try(lambda: parsed.lipid.headgroup.get_lipid_string())
    sum_formula = _try(parsed.get_sum_formula)

    # Chain counts come from the parse object, not from counting C:D tokens
    # in the rendered string (research R9, spec 011 T121): a one-chain class
    # renders with empty-position-slot tokens (MG 18:1 -> MG 18:1/0:0/0:0, 3
    # tokens for 1 real chain) and a sterol ester's backbone renders its own
    # pseudo-chain token (CE 18:1 -> SE 27:1/18:1, 2 tokens for 1 real
    # chain) -- both wrong under naive counting. lipid.fa_list has one entry
    # per possible chain slot; an unresolved/empty slot has num_carbon == 0
    # (excluded), while a genuinely named chain -- whether individually
    # resolved or a same-slot combined sum standing in for several
    # unresolved chains (the TG 58:12_20:0 case) -- has num_carbon > 0
    # (counted). No species downgrade: the old guard here replaced a
    # partial specification like TG 58:12_20:0 with its species total on the
    # premise that the subtracted/summed rendering misrepresents the real
    # chain count -- it does not (R9's arithmetic, verified independently in
    # omnipath-utils' mapping/_lipid.py); the guard discarded 3,468 distinct
    # partial specifications system-wide. "Partially specified" (chains
    # listed strictly between 0 and the class total) is its own level here
    # too, matching data-model.md section 7 exactly -- mirrored
    # independently from omnipath-utils' implementation (cross-repo
    # duplication is this codebase's established pattern for an identical
    # fix needed in two services, e.g. T071's name normalization).
    poss_fa = getattr(info, 'poss_fa', None)
    fa_list = getattr(parsed.lipid, 'fa_list', None) or []
    named_chains = sum(1 for fa in fa_list if getattr(fa, 'num_carbon', 0) > 0)
    if (
        isinstance(poss_fa, int) and poss_fa > 0
        and 0 < named_chains < poss_fa
    ):
        level_slug = 'partially_specified'

    total_carbon = total_db = None
    if species:
        chain = _CHAIN_RE.search(species)
        if chain:
            total_carbon, total_db = (int(x) for x in chain.group().split(':'))

    return {
        'normalised_name': normalised,
        'species_name': species,
        'lipid_level': level_slug,
        'lipid_category': category,
        'lipid_class': lipid_class,
        'total_carbon': total_carbon,
        'total_db': total_db,
        'sum_formula': sum_formula,
        'chains_possible': poss_fa,
        'chains_listed': named_chains,
    }


def _parse_many(names: list[str]) -> list[dict | None]:
    """Goslin-parse many names, parallel across cores for large batches.

    At full scale the candidate set is ~1.9M chain-bearing names; serial Goslin
    parsing took hours. ``parse_lipid`` is a module-level function with a
    per-process lazy ``LipidParser`` singleton, so a process pool parallelises
    cleanly — each worker builds its own parser once and results stay in input
    order. Small batches stay serial (pool setup isn't worth it).
    """
    if not names:
        return []
    import os

    workers = min(16, max(1, (os.cpu_count() or 2) - 2))
    if workers == 1 or len(names) < 2000:
        return [parse_lipid(n) for n in names]
    from multiprocessing import Pool

    with Pool(processes=workers) as pool:
        return list(pool.imap(parse_lipid, names, chunksize=1000))


@dataclass(frozen=True)
class LipidLabelStats:
    names_parsed: int = 0
    names_resolved: int = 0
    names_unresolved: int = 0
    lipids_labelled: int = 0
    by_level: dict[str, int] = field(default_factory=dict)


def _level_values() -> sql.SQL:
    return sql.SQL(', ').join(
        sql.SQL('({}, {})').format(sql.Literal(slug), sql.Literal(rank))
        for slug, rank in _LEVEL_RANK.items()
    )


def resolve_lipid_labels(
    conn,
    *,
    schema: str = 'public',
    max_records: int | None = None,
) -> LipidLabelStats:
    """Parse uncached lipid names, cache them, label lipid entities.

    Idempotent: only names absent from ``metabo_lipid_name_resolution`` are
    (re)parsed; the entity-label UPDATE is a full rebuild over the cache.
    """
    schema_id = sql.Identifier(schema)
    cache_id = sql.Identifier(schema, 'metabo_lipid_name_resolution')
    version = _goslin_version()

    with conn.cursor() as cur:
        chem_type_id = _chem_type_id(cur, schema_id)
        if not chem_type_id:
            return LipidLabelStats()

        # 1) Distinct chain-bearing candidate names not yet cached.
        limit = sql.SQL('LIMIT {}').format(sql.Literal(int(max_records))) \
            if max_records else sql.SQL('')
        cur.execute(
            sql.SQL(
                """
                SELECT DISTINCT btrim(ie.value) AS raw_name
                FROM {schema}.entity e
                JOIN {schema}.entity_identifier ei
                  ON ei.entity_id = e.entity_id
                JOIN {schema}.identifier_evidence ie
                  ON ie.identifier_id = ei.identifier_id
                JOIN {schema}.vocab_identifier_type it
                  ON it.identifier_type_id = ie.identifier_type_id
                WHERE e.entity_type_id = %(chem)s
                  AND it.name = ANY(%(name_types)s)
                  AND btrim(ie.value) ~ '[0-9]+:[0-9]+'
                  AND length(btrim(ie.value)) BETWEEN 2 AND 200
                  AND NOT EXISTS (
                    SELECT 1 FROM {cache} c WHERE c.raw_name = btrim(ie.value)
                  )
                {limit}
                """
            ).format(schema=schema_id, cache=cache_id, limit=limit),
            {'chem': chem_type_id, 'name_types': _NAME_TYPES},
        )
        raw_names = [row[0] for row in cur.fetchall()]

        # 2) Parse in Python (Goslin) — parallel across cores, then accumulate
        #    cache rows. (Serial parsing of the full-scale ~1.9M candidate set
        #    took hours; the NOT-EXISTS cache filter keeps re-runs incremental.)
        rows = []
        resolved = 0
        parsed_results = _parse_many(raw_names)
        for raw_name, parsed in zip(raw_names, parsed_results):
            if parsed is None:
                rows.append((raw_name, None, None, None, None, None,
                             None, None, None, 'unresolved', version,
                             None, None))
            else:
                resolved += 1
                rows.append((
                    raw_name,
                    parsed['normalised_name'],
                    parsed['species_name'],
                    parsed['lipid_level'],
                    parsed['lipid_category'],
                    parsed['lipid_class'],
                    parsed['total_carbon'],
                    parsed['total_db'],
                    parsed['sum_formula'],
                    'goslin',
                    version,
                    parsed['chains_possible'],
                    parsed['chains_listed'],
                ))
        if rows:
            execute_values(
                cur,
                sql.SQL(
                    """
                    INSERT INTO {cache} (
                      raw_name, normalised_name, species_name, lipid_level,
                      lipid_category, lipid_class, total_carbon, total_db,
                      sum_formula, resolver, goslin_version,
                      chains_possible, chains_listed
                    ) VALUES %s
                    ON CONFLICT (raw_name) DO NOTHING
                    """
                ).format(cache=cache_id).as_string(cur.connection),
                rows,
            )

        # 3) Label each lipid entity with its highest-granularity Goslin name.
        cur.execute(
            sql.SQL(
                """
                WITH lvl_rank(lipid_level, rank) AS (VALUES {levels}),
                cand AS (
                  SELECT
                    ei.entity_id,
                    c.normalised_name,
                    lr.rank
                  FROM {schema}.entity e
                  JOIN {schema}.entity_identifier ei
                    ON ei.entity_id = e.entity_id
                  JOIN {schema}.identifier_evidence ie
                    ON ie.identifier_id = ei.identifier_id
                  JOIN {schema}.vocab_identifier_type it
                    ON it.identifier_type_id = ie.identifier_type_id
                  JOIN {cache} c
                    ON c.raw_name = btrim(ie.value)
                  JOIN lvl_rank lr ON lr.lipid_level = c.lipid_level
                  WHERE e.entity_type_id = %(chem)s
                    AND it.name = ANY(%(name_types)s)
                    AND c.resolver = 'goslin'
                    AND c.normalised_name IS NOT NULL
                ),
                ranked AS (
                  SELECT entity_id, normalised_name,
                    row_number() OVER (
                      PARTITION BY entity_id
                      ORDER BY rank DESC,
                               length(normalised_name) DESC,
                               normalised_name
                    ) AS rk
                  FROM cand
                )
                UPDATE {schema}.entity e
                SET label = ranked.normalised_name,
                    label_rule = %(rule)s
                FROM ranked
                WHERE ranked.entity_id = e.entity_id
                  AND ranked.rk = 1
                  AND e.entity_type_id = %(chem)s
                """
            ).format(schema=schema_id, cache=cache_id, levels=_level_values()),
            {
                'chem': chem_type_id,
                'name_types': _NAME_TYPES,
                'rule': GOSLIN_LIPID_RULE,
            },
        )
        lipids_labelled = cur.rowcount

        # 4) Stats.
        cur.execute(
            sql.SQL(
                "SELECT resolver, count(*) FROM {cache} GROUP BY resolver"
            ).format(cache=cache_id)
        )
        counts = {r: int(n) for r, n in cur.fetchall()}
        cur.execute(
            sql.SQL(
                """
                SELECT c.lipid_level, count(DISTINCT e.entity_id)
                FROM {schema}.entity e
                JOIN {cache} c ON c.normalised_name = e.label
                WHERE e.entity_type_id = %(chem)s
                  AND e.label_rule = %(rule)s
                GROUP BY c.lipid_level
                """
            ).format(schema=schema_id, cache=cache_id),
            {'chem': chem_type_id, 'rule': GOSLIN_LIPID_RULE},
        )
        by_level = {lvl: int(n) for lvl, n in cur.fetchall()}

    conn.commit()
    return LipidLabelStats(
        names_parsed=len(raw_names),
        names_resolved=counts.get('goslin', 0),
        names_unresolved=counts.get('unresolved', 0),
        lipids_labelled=lipids_labelled,
        by_level=by_level,
    )


def _chem_type_id(cur, schema_id: sql.Identifier):
    cur.execute(
        sql.SQL(
            'SELECT entity_type_id FROM {}.vocab_entity_type WHERE name = %s'
        ).format(schema_id),
        [CHEMICAL_ENTITY_TYPE],
    )
    row = cur.fetchone()
    return row[0] if row else None
