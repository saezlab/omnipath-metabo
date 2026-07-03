#!/usr/bin/env python

"""Unit tests for the CellPhoneDB v5.0 COSMOS resource (US6).

All tests are fast (no network I/O). ``pypath.inputs_v2.cellphonedb`` is
mocked via ``patch.dict(sys.modules, ...)`` (KEGG test pattern), including
``download_proteins``/``iter_csv`` for the protein-type source (CellPhoneDB's
own ``protein_input.csv`` curation -- ``receptor`` column, ``'Transporter'``
tag -- not a generic cross-database classifier).
"""

import sys
from unittest.mock import MagicMock, patch

import pytest

from omnipath_metabo.datasets.cosmos._record import Interaction
from omnipath_metabo.datasets.cosmos.resources.cellphonedb import (
    SYNTH_RE,
    _build_complex_map,
    _build_protein_type_map,
    _load_cellphonedb_data,
    _load_cellphonedb_protein_types,
    _mor_from_modulatory_effect,
    _parse_synthetic_metabolite,
    _resolve_entity,
    cellphonedb_nonpeptidic_interactions,
    cellphonedb_ppi_interactions,
)


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------

def _interaction_row(
    partner_a='Adenosine_byNT5E_and_SLC29A1',
    partner_b='P11229',
    is_ppi='False',
    modulatory_effect='',
):
    return {
        'partner_a': partner_a,
        'partner_b': partner_b,
        'is_ppi': is_ppi,
        'modulatory_effect': modulatory_effect,
    }


def _complex_row(name='integrin_a2b1_complex', members=('P05556', 'P17301')):
    row = {
        'complex_name': name,
        'uniprot_1': '',
        'uniprot_2': '',
        'uniprot_3': '',
        'uniprot_4': '',
        'uniprot_5': '',
    }
    for i, m in enumerate(members, start=1):
        row[f'uniprot_{i}'] = m
    return row


def _protein_row(uniprot='P11229', receptor='FALSE', tags=''):
    return {'uniprot': uniprot, 'receptor': receptor, 'tags': tags}


def _mock_cellphonedb_module(interaction_rows, complex_rows, protein_rows=()):
    mock_resource = MagicMock()
    mock_resource.interactions.raw.return_value = iter(interaction_rows)
    mock_resource.complexes.raw.return_value = iter(complex_rows)

    mock_module = MagicMock()
    mock_module.resource = mock_resource
    mock_module.iter_csv.return_value = iter(protein_rows)
    return mock_module


def _run_nonpeptidic(interaction_rows, complex_rows=(), organism=9606, protein_rows=()):
    mock_module = _mock_cellphonedb_module(interaction_rows, complex_rows, protein_rows)
    with patch.dict(sys.modules, {'pypath.inputs_v2.cellphonedb': mock_module}):
        return list(cellphonedb_nonpeptidic_interactions(organism=organism))


def _run_ppi(interaction_rows, complex_rows=(), organism=9606):
    mock_module = _mock_cellphonedb_module(interaction_rows, complex_rows)
    with patch.dict(sys.modules, {'pypath.inputs_v2.cellphonedb': mock_module}):
        return list(cellphonedb_ppi_interactions(organism=organism))


# ---------------------------------------------------------------------------
# T002 -- _load_cellphonedb_data
# ---------------------------------------------------------------------------

class TestLoadCellphonedbData:

    def test_loads_interactions_and_complexes(self):
        interactions = [_interaction_row(), _interaction_row(is_ppi='True')]
        complexes = [_complex_row()]
        mock_module = _mock_cellphonedb_module(interactions, complexes)

        with patch.dict(sys.modules, {'pypath.inputs_v2.cellphonedb': mock_module}):
            loaded_interactions, loaded_complexes = _load_cellphonedb_data()

        assert loaded_interactions == interactions
        assert loaded_complexes == complexes

    def test_split_matches_r1_on_fixture_data(self):
        # 1029/1878/4 split scaled down: 3 non-peptidic matches, 2 PPI, 1 out-of-scope.
        interactions = [
            _interaction_row(partner_a='Adenosine_byNT5E_and_SLC29A1'),
            _interaction_row(partner_a='Acetylcholine_byCHAT'),
            _interaction_row(partner_a='Cholesterol_byDHCR7'),
            _interaction_row(partner_a='P12830', partner_b='Q96E93', is_ppi='True'),
            _interaction_row(partner_a='P19022', partner_b='P06734', is_ppi='True'),
            _interaction_row(partner_a='hCGB1_complex', partner_b='P22888', is_ppi='False'),
        ]
        mock_module = _mock_cellphonedb_module(interactions, [])

        with patch.dict(sys.modules, {'pypath.inputs_v2.cellphonedb': mock_module}):
            loaded, _ = _load_cellphonedb_data()

        non_peptidic = [
            r for r in loaded
            if r['is_ppi'] == 'False' and SYNTH_RE.match(r['partner_a'])
        ]
        ppi = [r for r in loaded if r['is_ppi'] == 'True']
        out_of_scope = [
            r for r in loaded
            if r['is_ppi'] == 'False' and not SYNTH_RE.match(r['partner_a'])
        ]
        assert len(non_peptidic) == 3
        assert len(ppi) == 2
        assert len(out_of_scope) == 1


# ---------------------------------------------------------------------------
# T003 -- _build_complex_map
# ---------------------------------------------------------------------------

class TestBuildComplexMap:

    def test_normal_complex_row(self):
        rows = [_complex_row('integrin_a2b1_complex', ('P05556', 'P17301'))]
        result = _build_complex_map(rows)
        assert result == {'integrin_a2b1_complex': ['P05556', 'P17301']}

    def test_synthetic_machinery_row_excluded(self):
        rows = [_complex_row('Dehydroepiandrosterone_bySTS', ('Q8TF42',))]
        result = _build_complex_map(rows)
        assert result == {}

    def test_null_padded_row_keeps_only_nonnull_members(self):
        rows = [_complex_row('complex_x', ('P00001',))]
        result = _build_complex_map(rows)
        assert result == {'complex_x': ['P00001']}

    def test_multiple_rows(self):
        rows = [
            _complex_row('complex_a', ('P00001', 'P00002')),
            _complex_row('complex_b', ('P00003',)),
        ]
        result = _build_complex_map(rows)
        assert result == {
            'complex_a': ['P00001', 'P00002'],
            'complex_b': ['P00003'],
        }


# ---------------------------------------------------------------------------
# _load_cellphonedb_protein_types / _build_protein_type_map
#
# CellPhoneDB's own protein_input.csv curation (receptor column,
# 'Transporter' tag) -- replaces the generic OmniPath Intercell/TCDB/G2P
# classifier previously reused from stitch.py's _multidb_uniprot_types.
# ---------------------------------------------------------------------------

class TestLoadCellphonedbProteinTypes:

    def test_loads_protein_rows(self):
        protein_rows = [_protein_row('P03372', receptor='TRUE')]
        mock_module = _mock_cellphonedb_module([], [], protein_rows)

        with patch.dict(sys.modules, {'pypath.inputs_v2.cellphonedb': mock_module}):
            loaded = _load_cellphonedb_protein_types()

        assert loaded == protein_rows


class TestBuildProteinTypeMap:

    def test_receptor_true_maps_to_receptor(self):
        rows = [_protein_row('P03372', receptor='TRUE')]
        assert _build_protein_type_map(rows) == {'P03372': 'receptor'}

    def test_transporter_tag_maps_to_transporter(self):
        rows = [_protein_row('P43003', receptor='FALSE', tags='Transporter')]
        assert _build_protein_type_map(rows) == {'P43003': 'transporter'}

    def test_compound_tag_containing_transporter_still_matches(self):
        rows = [_protein_row('P43003', receptor='FALSE', tags='Transporter|To_comment')]
        assert _build_protein_type_map(rows) == {'P43003': 'transporter'}

    def test_neither_flag_absent_from_map(self):
        rows = [_protein_row('O00341', receptor='FALSE', tags='')]
        assert _build_protein_type_map(rows) == {}

    def test_unrelated_tag_absent_from_map(self):
        rows = [_protein_row('P12345', receptor='FALSE', tags='Glycoprotein')]
        assert _build_protein_type_map(rows) == {}

    def test_multiple_rows(self):
        rows = [
            _protein_row('P03372', receptor='TRUE'),
            _protein_row('P43003', receptor='FALSE', tags='Transporter'),
            _protein_row('O00341', receptor='FALSE', tags=''),
        ]
        assert _build_protein_type_map(rows) == {
            'P03372': 'receptor',
            'P43003': 'transporter',
        }


# ---------------------------------------------------------------------------
# T004 -- SYNTH_RE / _parse_synthetic_metabolite
# ---------------------------------------------------------------------------

class TestSynthRe:

    def test_extracts_name_and_single_machinery_gene(self):
        result = _parse_synthetic_metabolite('Acetylcholine_byCHAT')
        assert result == ('Acetylcholine', ['CHAT'])

    def test_extracts_name_and_multiple_machinery_genes(self):
        result = _parse_synthetic_metabolite('Adenosine_byNT5E_and_SLC29A1')
        assert result == ('Adenosine', ['NT5E', 'SLC29A1'])

    def test_hcgb_complex_does_not_match(self):
        assert SYNTH_RE.match('hCGB1_complex') is None
        assert _parse_synthetic_metabolite('hCGB1_complex') is None


# ---------------------------------------------------------------------------
# T005's counterpart for _resolve_entity / _mor_from_modulatory_effect
# (shared helpers used by both slices)
# ---------------------------------------------------------------------------

class TestResolveEntity:

    def test_known_complex_resolves_to_members(self):
        complex_map = {'integrin_a2b1_complex': ['P05556', 'P17301']}
        assert _resolve_entity('integrin_a2b1_complex', complex_map) == ['P05556', 'P17301']

    def test_uniprot_ac_resolves_to_itself(self):
        assert _resolve_entity('P11229', {}) == ['P11229']

    def test_unknown_value_resolves_to_empty(self):
        assert _resolve_entity('not_a_thing', {}) == []


class TestMorFromModulatoryEffect:

    def test_inhibitory_is_minus_one(self):
        assert _mor_from_modulatory_effect('Inhibitory') == -1

    def test_empty_is_plus_one(self):
        assert _mor_from_modulatory_effect('') == 1

    def test_other_value_is_plus_one(self):
        assert _mor_from_modulatory_effect('Something else') == 1


# ---------------------------------------------------------------------------
# T010-T014 -- cellphonedb_nonpeptidic_interactions (Slice A)
# ---------------------------------------------------------------------------

class TestNonpeptidicInteractions:

    def test_yields_expected_fields(self):
        rows = [_interaction_row(partner_a='Adenosine_byNT5E_and_SLC29A1', partner_b='P11229')]
        recs = _run_nonpeptidic(rows)
        assert len(recs) == 1
        rec = recs[0]
        assert rec.source_type == 'small_molecule'
        assert rec.id_type_a == 'name'
        assert rec.source == 'Adenosine'
        assert rec.target == 'P11229'
        assert rec.attrs['producing_machinery'] == ['NT5E', 'SLC29A1']

    def test_no_extra_enzyme_to_metabolite_rows(self):
        rows = [_interaction_row(partner_a='Adenosine_byNT5E_and_SLC29A1', partner_b='P11229')]
        recs = _run_nonpeptidic(rows)
        # Only the metabolite->receptor edge; no edges involving NT5E/SLC29A1 directly.
        assert len(recs) == 1
        assert all(r.source not in ('NT5E', 'SLC29A1') for r in recs)
        assert all(r.target not in ('NT5E', 'SLC29A1') for r in recs)

    def test_ppi_rows_are_ignored(self):
        rows = [_interaction_row(is_ppi='True', partner_a='P12830', partner_b='P06734')]
        assert _run_nonpeptidic(rows) == []

    def test_transporter_routing(self):
        rows = [_interaction_row(partner_a='Adenosine_byNT5E_and_SLC29A1', partner_b='P11229')]
        recs = _run_nonpeptidic(rows, protein_rows=[_protein_row('P11229', tags='Transporter')])
        assert recs[0].interaction_type == 'transport'

    def test_receptor_routing(self):
        rows = [_interaction_row(partner_a='Adenosine_byNT5E_and_SLC29A1', partner_b='P11229')]
        recs = _run_nonpeptidic(rows, protein_rows=[_protein_row('P11229', receptor='TRUE')])
        assert recs[0].interaction_type == 'ligand_receptor'

    def test_unclassified_defaults_to_ligand_receptor(self):
        rows = [_interaction_row(partner_a='Adenosine_byNT5E_and_SLC29A1', partner_b='P11229')]
        recs = _run_nonpeptidic(rows, protein_rows=[])
        assert recs[0].interaction_type == 'ligand_receptor'

    def test_organism_scoping_nonhuman_yields_nothing(self):
        rows = [_interaction_row()]
        assert _run_nonpeptidic(rows, organism=10090) == []

    def test_complex_partner_b_expands_to_one_row_per_member(self):
        rows = [_interaction_row(
            partner_a='Adenosine_byNT5E_and_SLC29A1',
            partner_b='integrin_a2b1_complex',
        )]
        complexes = [_complex_row('integrin_a2b1_complex', ('P05556', 'P17301'))]
        recs = _run_nonpeptidic(rows, complex_rows=complexes)
        assert {r.target for r in recs} == {'P05556', 'P17301'}
        assert len(recs) == 2

    def test_unknown_complex_name_logged_and_skipped(self, caplog):
        import logging
        rows = [_interaction_row(
            partner_a='Adenosine_byNT5E_and_SLC29A1',
            partner_b='unknown_complex_xyz',
        )]
        with caplog.at_level(logging.INFO):
            recs = _run_nonpeptidic(rows)
        assert recs == []
        assert 'unknown_complex_xyz' in caplog.text

    def test_four_hcgb_rows_logged_and_skipped(self, caplog):
        import logging
        rows = [
            _interaction_row(partner_a='hCGB1_complex', partner_b='P22888'),
            _interaction_row(partner_a='hCGB2_complex', partner_b='P22888'),
            _interaction_row(partner_a='hCGB3_complex', partner_b='P22888'),
            _interaction_row(partner_a='hCGB7_complex', partner_b='P22888'),
        ]
        with caplog.at_level(logging.INFO):
            recs = _run_nonpeptidic(rows)
        assert recs == []
        for name in ('hCGB1_complex', 'hCGB2_complex', 'hCGB3_complex', 'hCGB7_complex'):
            assert name in caplog.text

    def test_result_are_interaction_instances(self):
        rows = [_interaction_row()]
        recs = _run_nonpeptidic(rows)
        assert all(isinstance(r, Interaction) for r in recs)


# ---------------------------------------------------------------------------
# T018-T022 -- cellphonedb_ppi_interactions (Slice B)
# ---------------------------------------------------------------------------

class TestPpiInteractions:

    def test_yields_expected_fields(self):
        rows = [_interaction_row(is_ppi='True', partner_a='P12830', partner_b='P06734')]
        recs = _run_ppi(rows)
        assert len(recs) == 1
        rec = recs[0]
        assert rec.source_type == 'protein'
        assert rec.target_type == 'protein'
        assert rec.id_type_a == 'uniprot'
        assert rec.id_type_b == 'uniprot'
        assert rec.interaction_type == 'ligand_receptor'
        assert rec.source == 'P12830'
        assert rec.target == 'P06734'
        assert rec.resource == 'CellPhoneDB'

    def test_nonpeptidic_rows_are_ignored(self):
        rows = [_interaction_row(is_ppi='False')]
        assert _run_ppi(rows) == []

    def test_mor_inhibitory(self):
        rows = [_interaction_row(
            is_ppi='True', partner_a='P12830', partner_b='P06734',
            modulatory_effect='Inhibitory',
        )]
        recs = _run_ppi(rows)
        assert recs[0].mor == -1

    def test_mor_absent_is_plus_one(self):
        rows = [_interaction_row(
            is_ppi='True', partner_a='P12830', partner_b='P06734',
            modulatory_effect='',
        )]
        recs = _run_ppi(rows)
        assert recs[0].mor == 1

    def test_organism_scoping_nonhuman_yields_nothing(self):
        rows = [_interaction_row(is_ppi='True', partner_a='P12830', partner_b='P06734')]
        assert _run_ppi(rows, organism=10116) == []

    def test_complex_expansion_reuses_complex_map(self):
        rows = [_interaction_row(
            is_ppi='True', partner_a='P12830', partner_b='integrin_a2b1_complex',
        )]
        complexes = [_complex_row('integrin_a2b1_complex', ('P05556', 'P17301'))]
        recs = _run_ppi(rows, complex_rows=complexes)
        assert {r.target for r in recs} == {'P05556', 'P17301'}
        assert len(recs) == 2

    def test_no_complex_name_string_in_source_or_target(self):
        rows = [
            _interaction_row(
                is_ppi='True', partner_a='P12830', partner_b='integrin_a2b1_complex',
            ),
            _interaction_row(
                partner_a='Adenosine_byNT5E_and_SLC29A1', partner_b='integrin_a2b1_complex',
                is_ppi='False',
            ),
        ]
        complexes = [_complex_row('integrin_a2b1_complex', ('P05556', 'P17301'))]
        ppi_recs = _run_ppi(rows, complex_rows=complexes)
        nonpep_recs = _run_nonpeptidic(rows, complex_rows=complexes)
        all_recs = ppi_recs + nonpep_recs
        assert all('integrin_a2b1_complex' not in (r.source, r.target) for r in all_recs)

    def test_unresolved_partner_logged_and_skipped(self, caplog):
        import logging
        rows = [_interaction_row(is_ppi='True', partner_a='HLAA', partner_b='P06734')]
        with caplog.at_level(logging.INFO):
            recs = _run_ppi(rows)
        assert recs == []
        assert 'HLAA' in caplog.text

    def test_result_are_interaction_instances(self):
        rows = [_interaction_row(is_ppi='True', partner_a='P12830', partner_b='P06734')]
        recs = _run_ppi(rows)
        assert all(isinstance(r, Interaction) for r in recs)
