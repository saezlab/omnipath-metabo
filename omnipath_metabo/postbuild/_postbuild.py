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
Orchestration of the omnipath-metabo post-build chemistry layer (Milestones E/F).

Runs after the main build completes, against the MAIN OmniPath Postgres:
ensure the RDKit cartridge + ``metabo_*`` schema, gate on the build manifest,
build the structure substrate, classify structural specificity + facets, and
classify RaMP conflicts. Idempotent (full rebuild of the precomputed tables);
refuses to write against a changed build unless ``force`` is set.
"""

from __future__ import annotations

__all__ = ['PostBuildStats', 'post_build_metabo']

from dataclasses import dataclass, field

from omnipath_metabo.db import (
    connect,
    ensure_metabo_schema,
    current_build_id,
    recorded_build_id,
    record_build_id,
)
from omnipath_metabo.postbuild._chem_layer import (
    build_structure_substrate,
    classify_structural_specificity,
    refresh_structural_specificity_facet,
)
from omnipath_metabo.postbuild._ramp_conflicts import populate_ramp_conflicts


@dataclass(frozen=True)
class PostBuildStats:
    build_id: str = ''
    structures: int = 0
    chemicals: int = 0
    specificity_by_level: dict[str, int] = field(default_factory=dict)
    facet_rows: int = 0
    ramp_conflicts: int = 0


def post_build_metabo(
    conn=None,
    *,
    db_url: str | None = None,
    schema: str = 'public',
    force: bool = False,
    conflicts: bool = True,
    conflict_max_records: int | None = None,
    log=print,
) -> PostBuildStats:
    """Run the post-build chemistry layer; return what was written.

    Pass an open ``conn`` or a ``db_url`` (else ``OMNIPATH_DB_URL`` /
    ``DATABASE_URL``). With ``force=False`` the step refuses to run when the main
    build id changed since the metabo layer was last computed (FR-018).
    """
    owns_conn = conn is None
    if conn is None:
        conn = connect(db_url)
    try:
        ensure_metabo_schema(conn, schema=schema)
        build_id = current_build_id(conn, schema)
        prior = recorded_build_id(conn, schema)
        if prior is not None and prior != build_id and not force:
            raise RuntimeError(
                f'metabo chemistry layer was computed for build {prior}, but the '
                f'main build is now {build_id}; rerun with force=True to recompute '
                f'for the new build.'
            )

        log(f'[post-build-metabo] build_id={build_id} schema={schema}')
        structures = build_structure_substrate(conn, schema=schema)
        log(f'[post-build-metabo] structure substrate: {structures} molecules')
        spec = classify_structural_specificity(conn, schema=schema)
        log(
            '[post-build-metabo] specificity: '
            + ' '.join(f'{k}={v}' for k, v in sorted(spec.by_level.items()))
        )
        facet_rows = refresh_structural_specificity_facet(conn, schema=schema)
        log(f'[post-build-metabo] structural_specificity facet rows: {facet_rows}')

        ramp_conflicts = 0
        if conflicts:
            stats = populate_ramp_conflicts(
                conn, schema=schema, max_records=conflict_max_records
            )
            ramp_conflicts = stats.conflicts
            log(
                f'[post-build-metabo] RaMP rows={stats.ramp_rows} '
                f'conflicts={ramp_conflicts}'
            )

        record_build_id(conn, schema, build_id)
        log(f'[post-build-metabo] done; recorded build_id={build_id}')
        return PostBuildStats(
            build_id=build_id,
            structures=structures,
            chemicals=spec.chemicals,
            specificity_by_level=spec.by_level,
            facet_rows=facet_rows,
            ramp_conflicts=ramp_conflicts,
        )
    finally:
        if owns_conn:
            conn.close()
