"""Projection-layer tests for /sets/metsigdb (cycles 010 and 012).

The projection layer owns the shape of a row and nothing about which rows come
back. No database: it is a pure mapping from a substrate row to the contract
row.

Cycle 012 split that shape in two. The default is thirteen fields; `fields=`
reaches the other nine, and `fields=all` is the cycle 010 response. The cycle
010 tests below are kept and, where they read a field that is now opt-in, they
ask for it.
"""

from __future__ import annotations

from uuid import UUID

import pytest

from omnipath_metabo.server.sets._metsigdb_projection import (
    ALL_FIELDS,
    DEFAULT_FIELDS,
    OPTIONAL_FIELDS,
    ROW_FIELDS,
    project_row,
    project_rows,
    resolve_fields,
)

SUBSTRATE_ROW = {
    'resource': 'Reactome',
    'set_source_id': 'R-HSA-1059683',
    'metabolite_entity_id': '971c119b-313f-f9f6-ce2e-5c7f0cec27ca',
    'metabolite_label': 'ATP',
    'metabolite_entity_type': 'Chemical:OM:0037',
    'metabolite_structure_key': 'ZKHQWZAMYRWXGA',
    'inchikey': None,
    'smiles': None,
    'hmdb': None,
    'pubchem': None,
    'chebi': '30616',
    'kegg': None,
    'set_entity_id': '8d61798c-34b5-ff31-734d-d65b8b1b0889',
    'set_label': 'Interleukin-6 signaling',
    'set_type': 'pathway',
    'organism': 9606,
    'set_size': 2,
    'set_context': None,
    'provenance_source': 'pypath.inputs_v2.reactome@6176be278088',
    'provenance_record': {'row_id': 3592, 'source_id': 32, 'dataset_id': 77},
    'build_id': '0ed5afcb26c6',
}

ALL = resolve_fields(['all'])


# --------------------------------------------------------------- the contract


def test_the_projection_publishes_the_contract_fields():
    assert set(project_row(SUBSTRATE_ROW, ALL)) == set(ROW_FIELDS) - {
        'metabolite_entity_id', 'set_source_id',
    } | {'entity', 'set'}


def test_field_order_follows_the_contract():
    assert list(project_row(SUBSTRATE_ROW, ALL)) == list(ALL_FIELDS)
    assert list(project_row(SUBSTRATE_ROW)) == list(DEFAULT_FIELDS)


def test_inchi_is_not_published():
    """No InChI identifier type exists, so the field left the contract."""
    assert 'inchi' not in ROW_FIELDS
    assert 'inchi' not in ALL_FIELDS


def test_the_default_projection_keeps_its_priority_order():
    identifiers = [f for f in ROW_FIELDS if f in {
        'inchikey', 'smiles', 'hmdb', 'pubchem', 'chebi', 'kegg',
    }]
    assert identifiers == ['inchikey', 'smiles', 'hmdb', 'pubchem', 'chebi', 'kegg']


def test_mandatory_fields_survive():
    row = project_row(SUBSTRATE_ROW, ALL)
    for field in (
        'entity', 'metabolite_label', 'metabolite_entity_type',
        'resource', 'set', 'set_type', 'set_size', 'provenance_source',
    ):
        assert row[field] is not None


def test_absent_identifiers_stay_present_and_null():
    """The schema does not change with the data, so a consumer can rely on it."""
    row = project_row(SUBSTRATE_ROW)
    assert row['inchikey'] is None
    assert 'inchikey' in row


def test_a_named_set_publishes_its_name():
    """Reactome, MACdb and ClassyFire name every set they publish."""
    assert project_row(SUBSTRATE_ROW)['set_label'] == 'Interleukin-6 signaling'


def test_an_unnamed_set_publishes_null():
    """KEGG and WikiPathways carry no name in this build."""
    assert project_row(SUBSTRATE_ROW | {'set_label': None})['set_label'] is None


def test_both_entity_ids_are_strings():
    """A uuid must not reach the response as an object the encoder guesses at."""
    row = project_row(SUBSTRATE_ROW, ALL)
    assert isinstance(row['entity'], str)
    assert isinstance(row['set_entity_id'], str)


def test_a_uuid_object_is_stringified():
    """The substrate hands over uuid objects; the response publishes strings."""
    stored = SUBSTRATE_ROW | {
        'metabolite_entity_id': UUID('971c119b-313f-f9f6-ce2e-5c7f0cec27ca'),
        'set_entity_id': UUID('8d61798c-34b5-ff31-734d-d65b8b1b0889'),
    }
    row = project_row(stored, ALL)
    assert row['entity'] == '971c119b-313f-f9f6-ce2e-5c7f0cec27ca'
    assert row['set_entity_id'] == '8d61798c-34b5-ff31-734d-d65b8b1b0889'


def test_the_structure_key_is_published():
    """It is what makes a cross-resource join on one molecule possible."""
    row = project_row(SUBSTRATE_ROW, resolve_fields(['metabolite_structure_key']))
    assert row['metabolite_structure_key'] == 'ZKHQWZAMYRWXGA'


def test_the_entity_id_is_a_string():
    """A uuid must not reach the response as an object the encoder guesses at."""
    assert isinstance(project_row(SUBSTRATE_ROW)['entity'], str)


def test_structured_metadata_stays_structured():
    row = project_row(
        SUBSTRATE_ROW | {'set_context': {'assignment': 'ancestor', 'depth': 3}},
        ALL,
    )
    assert row['set_context'] == {'assignment': 'ancestor', 'depth': 3}
    assert row['provenance_record']['source_id'] == 32


def test_a_missing_column_is_a_null_not_a_crash():
    """A substrate that gains a column must not break a running service."""
    trimmed = {k: v for k, v in SUBSTRATE_ROW.items() if k != 'smiles'}
    assert project_row(trimmed, resolve_fields(['smiles']))['smiles'] is None


def test_no_internal_column_leaks():
    """Mapping ambiguity is not a public concept, and neither is the staging."""
    row = project_row(SUBSTRATE_ROW | {'status_id': 1, 'via_reaction': 'x'}, ALL)
    assert 'status_id' not in row
    assert 'via_reaction' not in row


def test_projecting_many_rows_preserves_order():
    rows = project_rows([SUBSTRATE_ROW, SUBSTRATE_ROW | {'set_source_id': 'R-HSA-2'}])
    assert [r['set'] for r in rows] == ['R-HSA-1059683', 'R-HSA-2']


def test_an_empty_result_projects_to_an_empty_list():
    assert project_rows([]) == []


# ------------------------------------------------- the field registry (012)


def test_the_registry_partitions_the_published_row():
    """Thirteen by default, nine on request, twenty-two in total."""
    assert len(DEFAULT_FIELDS) == 13
    assert len(OPTIONAL_FIELDS) == 9
    assert len(ALL_FIELDS) == len(ROW_FIELDS) == 22
    assert set(DEFAULT_FIELDS) | set(OPTIONAL_FIELDS) == set(ALL_FIELDS)
    assert not set(DEFAULT_FIELDS) & set(OPTIONAL_FIELDS)


def test_the_default_projection_is_the_thirteen():
    """The point of the cycle: a response a consumer does not have to trim."""
    projected = project_row(SUBSTRATE_ROW)

    assert set(projected) == set(DEFAULT_FIELDS)
    # The two renames reach the response under the names their filters use.
    assert projected['entity'] == '971c119b-313f-f9f6-ce2e-5c7f0cec27ca'
    assert projected['set'] == 'R-HSA-1059683'
    assert 'metabolite_entity_id' not in projected
    assert 'set_source_id' not in projected
    # The two heaviest structured columns are absent unless asked for.
    assert 'set_context' not in projected
    assert 'provenance_record' not in projected


def test_named_fields_are_added_to_the_default():
    """`fields=` adds; it never replaces the default."""
    fields = resolve_fields(['smiles', 'set_context'])
    projected = project_row(SUBSTRATE_ROW, fields)

    assert set(projected) == set(DEFAULT_FIELDS) | {'smiles', 'set_context'}
    assert projected['resource'] == 'Reactome'


def test_naming_a_default_field_changes_nothing():
    """Accepted rather than refused, and it does not duplicate."""
    fields = resolve_fields(['resource', 'set_label'])
    assert fields == DEFAULT_FIELDS


def test_all_returns_every_published_field():
    """`fields=all` is the cycle 010 response, for a consumer that wants it."""
    projected = project_row(SUBSTRATE_ROW, resolve_fields(['all']))
    assert set(projected) == set(ALL_FIELDS)
    assert len(projected) == 22


def test_the_resolved_order_follows_the_contract():
    """A response's field order is the contract's, not the request's."""
    fields = resolve_fields(['set_context', 'smiles'])
    assert list(fields) == [f for f in ALL_FIELDS if f in set(fields)]


def test_an_unknown_field_is_refused_by_name():
    """`inchi` left the contract in cycle 010. Asking for it must not pass."""
    with pytest.raises(ValueError) as excinfo:
        resolve_fields(['inchi'])

    message = str(excinfo.value)
    assert 'inchi' in message
    # The message has to say what is allowed, or the caller guesses again.
    assert 'smiles' in message


def test_a_column_name_is_not_a_response_name():
    """The renamed fields are reachable under one name only, not both."""
    with pytest.raises(ValueError):
        resolve_fields(['metabolite_entity_id'])
    with pytest.raises(ValueError):
        resolve_fields(['set_source_id'])


def test_every_published_field_stays_reachable():
    """No published data may become unavailable through the API."""
    for field in ALL_FIELDS:
        assert field in resolve_fields([field])
