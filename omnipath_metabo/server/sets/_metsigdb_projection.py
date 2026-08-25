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

__all__ = [
    'ROW_FIELDS',
    'MetSigDBPage',
    'MetSigDBRow',
    'project_row',
    'project_rows',
]

from collections.abc import Iterable, Mapping
from typing import Any, TypedDict
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
    # Charge, stereo and tautomer variants of one molecule share this key. The
    # anchor is still the entity, so identity is unchanged; the key is what
    # lets a consumer join a membership from one resource to the same molecule
    # in another when the two picked different protonation states.
    'metabolite_structure_key',
    'inchikey',
    'smiles',
    'hmdb',
    'pubchem',
    'chebi',
    'kegg',
    'resource',
    'set_source_id',
    # The canonical entity behind the set. Names and cross-references join on
    # this, never on `set_source_id`: MACdb trait ids are bare integers that
    # collide with ChEBI ids in `entity.canonical_identifier`.
    'set_entity_id',
    'set_label',
    'set_type',
    'set_sub_type',
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
    for field in ('metabolite_entity_id', 'set_entity_id'):
        value = projected[field]
        if isinstance(value, UUID):
            projected[field] = str(value)

    return projected


def project_rows(rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Project a result set, keeping the query layer's order."""
    return [project_row(row) for row in rows]


class MetSigDBRow(TypedDict):
    """The published row, as a type.

    Declared so the OpenAPI document describes the response instead of calling
    it an object with unspecified properties. A client generator reads this and
    produces a row type; without it, every consumer hand-writes these fields
    again.

    A TypedDict rather than a dataclass on purpose: the runtime value stays a
    plain dict, so a page of 100,000 rows costs no object construction.
    """

    metabolite_entity_id: str
    metabolite_label: str
    metabolite_entity_type: str
    metabolite_structure_key: str | None
    inchikey: str | None
    smiles: str | None
    hmdb: str | None
    pubchem: str | None
    chebi: str | None
    kegg: str | None
    resource: str
    set_source_id: str
    set_entity_id: str
    set_label: str | None
    set_type: str
    set_sub_type: str | None
    organism: int | None
    set_size: int
    set_context: dict[str, Any] | None
    provenance_source: str
    provenance_record: dict[str, Any] | None
    build_id: str


class MetSigDBPage(TypedDict):
    """One page of membership rows.

    ``count`` is the number of rows in *this* page, never the size of the whole
    result. ``has_more`` answers the only question a paging client really has,
    and costs one extra row rather than a count. ``total`` is present only when
    the request asks for it, because counting a filter that matches three
    million rows is work nobody should pay for by default.
    """

    count: int
    has_more: bool
    rows: list[MetSigDBRow]
    total: int | None
