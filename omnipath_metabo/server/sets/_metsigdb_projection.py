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
    'ALL_FIELDS',
    'MEMBER_FIELDS',
    'SET_FIELDS',
    'COLUMN_NAMES',
    'DEFAULT_FIELDS',
    'OPTIONAL_FIELDS',
    'RESPONSE_NAMES',
    'ROW_FIELDS',
    'MetSigDBGroupedPage',
    'MetSigDBPage',
    'MetSigDBRow',
    'ResourceGroup',
    'SetGroup',
    'SetGroupMember',
    'columns_for',
    'project_row',
    'project_rows',
    'resolve_fields',
]

from collections.abc import Iterable, Mapping
from typing import Any, NotRequired, TypedDict
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


# Two fields reach the response under the name of the filter that selects
# them, not the name of the column that stores them. The substrate's column
# names are unchanged; this is presentation only, so a consumer filtering on
# `entity` reads `entity` back.
RESPONSE_NAMES: dict[str, str] = {
    'metabolite_entity_id': 'entity',
    'set_source_id': 'set',
}

# The reverse, for the query layer: which column a response field comes from.
COLUMN_NAMES: dict[str, str] = {
    response: column for column, response in RESPONSE_NAMES.items()
}

# Every published field, in contract order, under its response name.
ALL_FIELDS: tuple[str, ...] = tuple(
    RESPONSE_NAMES.get(column, column) for column in ROW_FIELDS
)

# What a response carries when the request says nothing. Enough to identify the
# metabolite, name the set, and join to a consumer's own data. The nine left
# out are not less valid, only less often wanted, and `fields` reaches them.
DEFAULT_FIELDS: tuple[str, ...] = (
    'entity',
    'metabolite_label',
    'inchikey',
    'hmdb',
    'pubchem',
    'chebi',
    'kegg',
    'resource',
    'set',
    'set_label',
    'set_type',
    'set_sub_type',
    'set_size',
)

OPTIONAL_FIELDS: tuple[str, ...] = tuple(
    field for field in ALL_FIELDS if field not in DEFAULT_FIELDS
)

# `fields=all` is the cycle 010 response, for a consumer that wants it back.
ALL_KEYWORD = 'all'

# A grouped response splits the row in two. The set side describes the set and
# is written once per set; the metabolite side describes one member and is
# written once per member. Repeating six set-level fields on each of a set's
# members would undo the projection work in the same response that does it —
# ClassyFire's root class alone has 145,937 of them.
SET_FIELDS: tuple[str, ...] = (
    'resource',
    'set',
    'set_entity_id',
    'set_label',
    'set_type',
    'set_sub_type',
    'organism',
    'set_size',
    'set_context',
)

MEMBER_FIELDS: tuple[str, ...] = tuple(
    field for field in ALL_FIELDS if field not in SET_FIELDS
)

# Fields whose stored value is a uuid. They must reach the response as strings:
# left as objects the shape depends on whatever the encoder decides, which is
# not a contract.
_UUID_FIELDS = frozenset({'entity', 'set_entity_id'})


def resolve_fields(requested: Iterable[str] | None = None) -> tuple[str, ...]:
    """The fields one response carries, in contract order.

    `fields` **adds to** the default rather than replacing it, so a consumer
    naming one extra column does not lose the twelve they did not mention.
    Naming a default field is accepted and changes nothing.

    An unknown name raises rather than being dropped. Cycle 010 shipped a route
    that ignored an unknown *parameter* and answered the wrong question with a
    200; the same silence about an unknown *field* would be the same defect in a
    smaller place.
    """
    if requested is None:
        return DEFAULT_FIELDS

    requested = tuple(requested)
    if not requested:
        return DEFAULT_FIELDS

    if ALL_KEYWORD in requested:
        return ALL_FIELDS

    unknown = sorted({name for name in requested if name not in ALL_FIELDS})
    if unknown:
        raise ValueError(
            f'Unsupported field: {", ".join(unknown)}. '
            f'Supported: {", ".join(ALL_FIELDS)}, or "{ALL_KEYWORD}".'
        )

    wanted = set(DEFAULT_FIELDS) | set(requested)
    return tuple(field for field in ALL_FIELDS if field in wanted)


def columns_for(fields: Iterable[str]) -> tuple[str, ...]:
    """The substrate columns a set of response fields needs.

    The query selects these by name, so a response that asks for thirteen
    fields does not read twenty-two off the disk and discard nine.
    """
    return tuple(COLUMN_NAMES.get(field, field) for field in fields)


def project_row(
    row: Mapping[str, Any],
    fields: Iterable[str] | None = None,
) -> dict[str, Any]:
    """One substrate row as one contract row.

    Every requested field appears, absent ones as null, so the response schema
    does not change with the data. Anything the substrate carries beyond the
    request stays out: internal columns are not a public concept, and neither
    are columns nobody asked for.
    """
    fields = DEFAULT_FIELDS if fields is None else tuple(fields)

    projected: dict[str, Any] = {}
    for field in fields:
        value = row.get(COLUMN_NAMES.get(field, field))
        if field in _UUID_FIELDS and isinstance(value, UUID):
            value = str(value)
        projected[field] = value

    return projected


def project_rows(
    rows: Iterable[Mapping[str, Any]],
    fields: Iterable[str] | None = None,
) -> list[dict[str, Any]]:
    """Project a result set, keeping the query layer's order."""
    fields = DEFAULT_FIELDS if fields is None else tuple(fields)
    return [project_row(row, fields) for row in rows]


class MetSigDBRow(TypedDict):
    """The published row, as a type.

    Declared so the OpenAPI document describes the response instead of calling
    it an object with unspecified properties. A client generator reads this and
    produces a row type; without it, every consumer hand-writes these fields
    again.

    A TypedDict rather than a dataclass on purpose: the runtime value stays a
    plain dict, so a page of 100,000 rows costs no object construction.
    """

    # The thirteen a response carries when the request says nothing.
    entity: str
    metabolite_label: str
    inchikey: str | None
    hmdb: str | None
    pubchem: str | None
    chebi: str | None
    kegg: str | None
    resource: str
    set: str
    set_label: str | None
    set_type: str
    set_sub_type: str | None
    set_size: int

    # The nine `fields` reaches. NotRequired rather than optional-valued: they
    # are absent from the response, not present and null, so a client can tell
    # "not asked for" from "asked for and empty".
    metabolite_entity_type: NotRequired[str]
    metabolite_structure_key: NotRequired[str | None]
    smiles: NotRequired[str | None]
    set_entity_id: NotRequired[str]
    organism: NotRequired[int | None]
    set_context: NotRequired[dict[str, Any] | None]
    provenance_source: NotRequired[str]
    provenance_record: NotRequired[dict[str, Any] | None]
    build_id: NotRequired[str]


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


class SetGroupMember(TypedDict):
    """One metabolite inside a grouped set.

    The metabolite side of the row only. The set side is written once on the
    set, not repeated on each of its members.
    """

    entity: str
    metabolite_label: str
    inchikey: NotRequired[str | None]
    hmdb: NotRequired[str | None]
    pubchem: NotRequired[str | None]
    chebi: NotRequired[str | None]
    kegg: NotRequired[str | None]
    metabolite_entity_type: NotRequired[str]
    metabolite_structure_key: NotRequired[str | None]
    smiles: NotRequired[str | None]
    provenance_source: NotRequired[str]
    provenance_record: NotRequired[dict[str, Any] | None]
    build_id: NotRequired[str]


class SetGroup(TypedDict):
    """One set, complete, with its members.

    ``set_size`` is the set's **published population**; ``returned`` is how many
    members this response carries for it. Under a member-level filter the two
    differ, and reporting only one would make a filtered view read as a shrunken
    set rather than a partial view of a whole one.
    """

    set: str
    set_label: NotRequired[str | None]
    set_type: NotRequired[str]
    set_sub_type: NotRequired[str | None]
    set_entity_id: NotRequired[str]
    organism: NotRequired[int | None]
    set_context: NotRequired[dict[str, Any] | None]
    set_size: NotRequired[int]
    returned: int
    members: list[SetGroupMember]


class ResourceGroup(TypedDict):
    """The sets one resource contributed to this page."""

    resource: str
    sets: list[SetGroup]


class MetSigDBGroupedPage(TypedDict):
    """One page of **sets**, grouped by resource.

    Every count here is in sets, not rows. ``count`` is the sets in this page,
    ``total`` the sets the filter matched, and ``has_more`` says whether further
    sets follow. The flat page counts rows for all three, which is why grouping
    is a separate response type rather than a flag on the same one.
    """

    count: int
    has_more: bool
    groups: list[ResourceGroup]
    total: int | None
