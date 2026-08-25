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
Projection layer for ``/sets/metsigdb``.

Shapes a substrate row into the published row contract. It owns the shape of a
row and nothing about which rows come back.
"""

from __future__ import annotations

__all__ = ['ROW_FIELDS', 'project_row', 'project_rows']

from collections.abc import Iterable, Mapping
from typing import Any
from uuid import UUID

# The published row, in the order the row contract lists it. The metabolite
# side leads, then the set side, then provenance.
#
# `inchi` is absent by design: no InChI identifier type exists in the build, so
# the field left the contract. The six identifiers that remain are in the
# published priority order.
ROW_FIELDS: tuple[str, ...] = (
    'metabolite_entity_id',
    'metabolite_label',
    'metabolite_entity_type',
    'inchikey',
    'smiles',
    'hmdb',
    'pubchem',
    'chebi',
    'kegg',
    'resource',
    'set_source_id',
    'set_label',
    'set_type',
    'organism',
    'set_size',
    'set_context',
    'provenance_source',
    'provenance_record',
    'build_id',
)


def project_row(row: Mapping[str, Any]) -> dict[str, Any]:
    """One substrate row as one contract row.

    Every contract field appears, absent ones as null, so the response schema
    does not change with the data. Anything the substrate carries beyond the
    contract stays out: internal columns are not a public concept.
    """
    projected = {field: row.get(field) for field in ROW_FIELDS}

    # A uuid must reach the response as a string. Left as an object it depends
    # on whatever the encoder decides, which is not a contract.
    entity_id = projected['metabolite_entity_id']
    if isinstance(entity_id, UUID):
        projected['metabolite_entity_id'] = str(entity_id)

    return projected


def project_rows(rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Project a result set, keeping the query layer's order."""
    return [project_row(row) for row in rows]
