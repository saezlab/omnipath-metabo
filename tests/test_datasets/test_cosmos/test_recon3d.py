#!/usr/bin/env python

"""Unit tests for the Recon3D transporter resource.

All tests here are fast (no network I/O) — pypath's
``recon3d_transporter_network`` is mocked with synthetic GemInteraction
records.  Integration tests that download real Recon3D data live in
``test_resources.py`` under the ``slow`` mark.
"""

from unittest.mock import patch

import pytest

from pypath.inputs.recon3d._gem import _parse_gene_rule, _strip_compartment
from pypath.inputs.metatlas._records import GemInteraction

from omnipath_metabo.datasets.cosmos.resources.recon3d import (
    recon3d_transporter_interactions,
)
from omnipath_metabo.datasets.cosmos._record import Interaction


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _gi(src, tgt, src_type='metabolite', tgt_type='protein',
        src_comp='', tgt_comp='', rxn_id='RXN1', reverse=False):
    return GemInteraction(
        source=src, target=tgt,
        source_type=src_type, target_type=tgt_type,
        source_compartment=src_comp, target_compartment=tgt_comp,
        reaction_id=rxn_id, reverse=reverse,
    )


def _run(gi_records, **kwargs):
    """Run recon3d_transporter_interactions with mocked pypath transport network."""

    with patch(
        'pypath.inputs.recon3d._gem.recon3d_transporter_network',
        return_value=iter(gi_records),
    ):
        return list(recon3d_transporter_interactions(**kwargs))


# ---------------------------------------------------------------------------
# _parse_gene_rule  (now in pypath.inputs.recon3d._gem)
# ---------------------------------------------------------------------------

class TestParseGeneRule:

    def test_empty_rule_returns_empty(self):
        assert _parse_gene_rule('') == []

    def test_whitespace_only_returns_empty(self):
        assert _parse_gene_rule('   ') == []

    def test_none_returns_empty(self):
        assert _parse_gene_rule(None) == []

    def test_single_gene_strips_suffix(self):
        result = _parse_gene_rule('1234_AT1')
        assert result == ['1234']

    def test_multiple_atn_suffixes_stripped(self):
        assert _parse_gene_rule('1234_AT2') == ['1234']

    def test_or_rule_gives_multiple_enzymes(self):
        result = _parse_gene_rule('1234_AT1 or 5678_AT1')
        assert set(result) == {'1234', '5678'}

    def test_and_rule_gives_complex(self):
        result = _parse_gene_rule('1234_AT1 and 5678_AT1')
        assert result == ['1234_5678']

    def test_parentheses_stripped(self):
        assert len(_parse_gene_rule('(1234_AT1 or 5678_AT1)')) == 2

    def test_mixed_or_and(self):
        result = _parse_gene_rule('1111_AT1 or (2222_AT1 and 3333_AT1)')
        assert set(result) == {'1111', '2222_3333'}

    def test_complex_subunits_sorted(self):
        result = _parse_gene_rule('9999_AT1 and 2222_AT1')
        assert result == ['2222_9999']

    def test_no_suffix_unchanged(self):
        assert _parse_gene_rule('1234') == ['1234']

    def test_strip_isoforms_false_preserves_suffix(self):
        result = _parse_gene_rule('1234_AT1', strip_isoforms=False)
        assert result == ['1234_AT1']


# ---------------------------------------------------------------------------
# _strip_compartment  (now in pypath.inputs.recon3d._gem)
# ---------------------------------------------------------------------------

class TestStripCompartment:

    def test_standard_bigg_id(self):
        assert _strip_compartment('atp_c') == ('atp', 'c')

    def test_extracellular(self):
        assert _strip_compartment('atp_e') == ('atp', 'e')

    def test_no_compartment_suffix(self):
        base, comp = _strip_compartment('atp')
        assert base == 'atp'
        assert comp == ''

    def test_multichar_suffix_not_stripped(self):
        base, comp = _strip_compartment('glc_e1')
        assert base == 'glc_e1'
        assert comp == ''

    def test_compound_id(self):
        base, comp = _strip_compartment('M00001_c')
        assert base == 'M00001'
        assert comp == 'c'


# ---------------------------------------------------------------------------
# Transport detection (via mocked recon3d_transporter_network)
# ---------------------------------------------------------------------------

class TestTransportDetection:

    def test_transport_reaction_yields_edges(self):
        gis = [
            _gi('atp', '1234', src_comp='c', rxn_id='EX_atp'),
            _gi('1234', 'atp', src_type='protein', tgt_type='small_molecule',
                tgt_comp='e', rxn_id='EX_atp'),
        ]
        recs = _run(gis)
        assert len(recs) > 0

    def test_no_input_yields_no_edges(self):
        assert _run([]) == []

    def test_met_to_enzyme_edge(self):
        gis = [
            _gi('atp', '1234', src_type='metabolite', tgt_type='protein',
                src_comp='c', rxn_id='EX_atp'),
        ]
        recs = _run(gis)
        met_edges = [r for r in recs if r.source_type == 'small_molecule']
        assert len(met_edges) == 1
        assert met_edges[0].id_type_a == 'bigg'
        assert met_edges[0].id_type_b == 'entrez'
        assert met_edges[0].source == 'atp'
        assert met_edges[0].target == '1234'

    def test_enzyme_to_met_edge(self):
        gis = [
            _gi('1234', 'atp', src_type='protein', tgt_type='metabolite',
                tgt_comp='e', rxn_id='EX_atp'),
        ]
        recs = _run(gis)
        enz_edges = [r for r in recs if r.source_type == 'protein']
        assert len(enz_edges) == 1
        assert enz_edges[0].id_type_a == 'entrez'
        assert enz_edges[0].id_type_b == 'bigg'

    def test_interaction_type_is_transport(self):
        gis = [_gi('atp', '1234', src_comp='c', rxn_id='EX_atp')]
        recs = _run(gis)
        assert all(r.interaction_type == 'transport' for r in recs)

    def test_resource_is_recon3d(self):
        gis = [_gi('atp', '1234', src_comp='c', rxn_id='EX_atp')]
        recs = _run(gis)
        assert all(r.resource == 'Recon3D' for r in recs)

    def test_transport_from_to_in_attrs(self):
        gis = [
            _gi('atp', '1234', src_type='metabolite', tgt_type='protein',
                src_comp='c', rxn_id='EX_atp'),
            _gi('1234', 'atp', src_type='protein', tgt_type='metabolite',
                tgt_comp='e', rxn_id='EX_atp'),
        ]
        recs = _run(gis)
        for r in recs:
            assert 'transport_from' in r.attrs
            assert 'transport_to' in r.attrs

    def test_compartment_in_locations(self):
        gis = [
            _gi('atp', '1234', src_type='metabolite', src_comp='c', rxn_id='EX_atp'),
        ]
        recs = _run(gis)
        for r in recs:
            assert isinstance(r.locations, tuple)

    def test_enzyme_complex_detected(self):
        gis = [
            _gi('glc', '5678_9012', src_type='metabolite', tgt_type='protein',
                src_comp='c', rxn_id='GLC_TRANS'),
        ]
        recs = _run(gis)
        complex_edges = [r for r in recs if r.attrs.get('enzyme_complex')]
        assert len(complex_edges) > 0


# ---------------------------------------------------------------------------
# Reverse edges
# ---------------------------------------------------------------------------

class TestReverseEdges:

    def test_reverse_attr_present(self):
        gis = [
            _gi('glc', '5678', src_comp='c', rxn_id='GLC_TRANS', reverse=False),
            _gi('5678', 'glc', src_type='protein', tgt_type='metabolite',
                tgt_comp='e', rxn_id='GLC_TRANS', reverse=False),
            _gi('glc', '5678', src_comp='e', rxn_id='GLC_TRANS', reverse=True),
            _gi('5678', 'glc', src_type='protein', tgt_type='metabolite',
                tgt_comp='c', rxn_id='GLC_TRANS', reverse=True),
        ]
        recs = _run(gis)
        reverse_flags = {r.attrs['reverse'] for r in recs}
        assert True in reverse_flags
        assert False in reverse_flags

    def test_include_reverse_false_no_reverse_edges(self):
        # When include_reverse=False the mock returns no reverse records.
        gis = [
            _gi('glc', '5678', src_comp='c', rxn_id='GLC_TRANS', reverse=False),
            _gi('5678', 'glc', src_type='protein', tgt_type='metabolite',
                tgt_comp='e', rxn_id='GLC_TRANS', reverse=False),
        ]
        recs = _run(gis, include_reverse=False)
        assert all(not r.attrs['reverse'] for r in recs)


# ---------------------------------------------------------------------------
# Orphan reactions
# ---------------------------------------------------------------------------

class TestOrphanReactions:

    def _orphan_gis(self, rxn_id='ORPHAN_NA'):
        return [
            _gi('na', rxn_id, src_type='metabolite', tgt_type='reaction',
                src_comp='c', rxn_id=rxn_id),
            _gi(rxn_id, 'na', src_type='reaction', tgt_type='metabolite',
                tgt_comp='e', rxn_id=rxn_id),
        ]

    def test_orphan_included(self):
        recs = _run(self._orphan_gis())
        assert len(recs) > 0

    def test_orphan_attrs_flag(self):
        recs = _run(self._orphan_gis())
        assert all(r.attrs.get('orphan') for r in recs)

    def test_orphan_uses_reaction_id_as_enzyme(self):
        recs = _run(self._orphan_gis('ORPHAN_NA'))
        protein_sides = [
            r.target if r.source_type == 'small_molecule' else r.source
            for r in recs
        ]
        assert all(p == 'ORPHAN_NA' for p in protein_sides)

    def test_orphan_id_type_is_reaction_id(self):
        recs = _run(self._orphan_gis())
        for r in recs:
            if r.source_type == 'small_molecule':
                assert r.id_type_b == 'reaction_id'
            else:
                assert r.id_type_a == 'reaction_id'

    def test_normal_edges_have_no_orphan_flag(self):
        gis = [_gi('atp', '1234', src_comp='c', rxn_id='EX_atp')]
        recs = _run(gis)
        assert all(not r.attrs.get('orphan') for r in recs)


# ---------------------------------------------------------------------------
# Organism filter
# ---------------------------------------------------------------------------

class TestOrganism:

    def test_non_human_yields_nothing(self):
        gis = [_gi('atp', '1234', src_comp='c')]
        assert _run(gis, organism=10090) == []

    def test_human_proceeds_normally(self):
        gis = [_gi('glc', '1234', src_comp='e')]
        assert len(_run(gis, organism=9606)) > 0


# ---------------------------------------------------------------------------
# cell_surface_only
# ---------------------------------------------------------------------------

class TestCellSurfaceOnly:

    def _plasma_gis(self):
        return [
            _gi('glc', '1234', src_comp='c', rxn_id='PLASMA'),
            _gi('1234', 'glc', src_type='protein', tgt_type='metabolite',
                tgt_comp='e', rxn_id='PLASMA'),
        ]

    def _mito_gis(self):
        return [
            _gi('atp', '5678', src_comp='c', rxn_id='MITO'),
            _gi('5678', 'atp', src_type='protein', tgt_type='metabolite',
                tgt_comp='m', rxn_id='MITO'),
        ]

    def test_plasma_membrane_kept(self):
        recs = _run(self._plasma_gis(), cell_surface_only=True)
        assert len(recs) > 0

    def test_intracellular_dropped(self):
        recs = _run(self._mito_gis(), cell_surface_only=True)
        assert recs == []

    def test_mixed_only_plasma_survives(self):
        recs = _run(self._plasma_gis() + self._mito_gis(), cell_surface_only=True)
        rxn_ids = {r.attrs['reaction_id'] for r in recs}
        assert 'PLASMA' in rxn_ids
        assert 'MITO' not in rxn_ids

    def test_false_keeps_intracellular(self):
        recs = _run(self._mito_gis(), cell_surface_only=False)
        assert len(recs) > 0
