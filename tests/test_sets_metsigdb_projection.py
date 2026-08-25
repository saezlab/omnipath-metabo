"""Projection-layer tests for /sets/metsigdb (cycle 010).

The projection layer owns the shape of a row and nothing about which rows come
back. No database: it is a pure mapping from a substrate row to the contract
row.
"""

from __future__ import annotations

import pytest

from omnipath_metabo.server.sets._metsigdb_projection import (
    ROW_FIELDS,
    project_row,
    project_rows,
)

SUBSTRATE_ROW = {
    'resource': 'Reactome',
    'set_source_id': 'R-HSA-1059683',
    'metabolite_entity_id': '971c119b-313f-f9f6-ce2e-5c7f0cec27ca',
    'metabolite_label': 'ATP',
    'metabolite_entity_type': 'Chemical:OM:0037',
    'inchikey': None,
    'smiles': None,
    'hmdb': None,
    'pubchem': None,
    'chebi': '30616',
    'kegg': None,
    'set_label': None,
    'set_type': 'pathway',
    'organism': 9606,
    'set_size': 2,
    'set_context': None,
    'provenance_source': 'pypath.inputs_v2.reactome@6176be278088',
    'provenance_record': {'row_id': 3592, 'source_id': 32, 'dataset_id': 77},
    'build_id': '0ed5afcb26c6',
}


def test_the_projection_publishes_the_contract_fields():
    assert set(project_row(SUBSTRATE_ROW)) == set(ROW_FIELDS)


def test_field_order_follows_the_contract():
    assert list(project_row(SUBSTRATE_ROW)) == list(ROW_FIELDS)


def test_inchi_is_not_published():
    """No InChI identifier type exists, so the field left the contract."""
    assert 'inchi' not in ROW_FIELDS


def test_the_default_projection_keeps_its_priority_order():
    identifiers = [f for f in ROW_FIELDS if f in {
        'inchikey', 'smiles', 'hmdb', 'pubchem', 'chebi', 'kegg',
    }]
    assert identifiers == ['inchikey', 'smiles', 'hmdb', 'pubchem', 'chebi', 'kegg']


def test_mandatory_fields_survive():
    row = project_row(SUBSTRATE_ROW)
    for field in (
        'metabolite_entity_id', 'metabolite_label', 'metabolite_entity_type',
        'resource', 'set_source_id', 'set_type', 'set_size', 'provenance_source',
    ):
        assert row[field] is not None


def test_absent_identifiers_stay_present_and_null():
    """The schema does not change with the data, so a consumer can rely on it."""
    row = project_row(SUBSTRATE_ROW)
    assert row['inchikey'] is None
    assert 'inchikey' in row


def test_null_set_label_is_published_as_null():
    """ontology_terms is empty, so no set carries a readable name in v1."""
    assert project_row(SUBSTRATE_ROW)['set_label'] is None


def test_the_entity_id_is_a_string():
    """A uuid must not reach the response as an object the encoder guesses at."""
    assert isinstance(project_row(SUBSTRATE_ROW)['metabolite_entity_id'], str)


def test_structured_metadata_stays_structured():
    row = project_row(
        SUBSTRATE_ROW | {'set_context': {'assignment': 'ancestor', 'depth': 3}}
    )
    assert row['set_context'] == {'assignment': 'ancestor', 'depth': 3}
    assert row['provenance_record']['source_id'] == 32


def test_a_missing_column_is_a_null_not_a_crash():
    """A substrate that gains a column must not break a running service."""
    trimmed = {k: v for k, v in SUBSTRATE_ROW.items() if k != 'smiles'}
    assert project_row(trimmed)['smiles'] is None


def test_no_internal_column_leaks():
    """Mapping ambiguity is not a public concept, and neither is the staging."""
    row = project_row(SUBSTRATE_ROW | {'status_id': 1, 'via_reaction': 'x'})
    assert 'status_id' not in row
    assert 'via_reaction' not in row


def test_projecting_many_rows_preserves_order():
    rows = project_rows([SUBSTRATE_ROW, SUBSTRATE_ROW | {'set_source_id': 'R-HSA-2'}])
    assert [r['set_source_id'] for r in rows] == ['R-HSA-1059683', 'R-HSA-2']


def test_an_empty_result_projects_to_an_empty_list():
    assert project_rows([]) == []
