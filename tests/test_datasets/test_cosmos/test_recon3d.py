#!/usr/bin/env python

"""Unit tests for the Recon3D transporter resource.

All tests here are fast (no network I/O) — pypath's ``recon3d_reactions``
is mocked with synthetic reaction dicts that mirror the BiGG JSON format.
Integration tests that download real Recon3D data live in
``test_resources.py`` under the ``slow`` mark.
"""

from unittest.mock import patch

import pytest

from omnipath_metabo.datasets.cosmos.resources.recon3d import (
    _parse_gene_rule,
    _strip_compartment,
    recon3d_transporter_interactions,
)
from omnipath_metabo.datasets.cosmos._record import Interaction


# ---------------------------------------------------------------------------
# Synthetic reaction fixtures
# ---------------------------------------------------------------------------

# A directed ATP transport (cytoplasm → extracellular).
_RXN_ATP_TRANSPORT = {
    'id': 'EX_atp',
    'gene_reaction_rule': '1234_AT1',
    'metabolites': {'atp_c': -1, 'atp_e': 1},
    'lower_bound': 0,
    'upper_bound': 1000,
    'reversible': False,
}

# A reversible glucose transport with an enzyme complex.
_RXN_GLC_TRANSPORT_REV = {
    'id': 'GLC_TRANS',
    'gene_reaction_rule': '5678_AT1 and 9012_AT1',
    'metabolites': {'glc_c': -1, 'glc_e': 1},
    'lower_bound': -1000,
    'upper_bound': 1000,
    'reversible': True,
}

# An intracellular reaction — same compartment on both sides, NOT a transport.
_RXN_INTRACELL = {
    'id': 'ATP_HYDROLYSIS',
    'gene_reaction_rule': '1234_AT1',
    'metabolites': {'atp_c': -1, 'adp_c': 1},
    'lower_bound': 0,
    'upper_bound': 1000,
    'reversible': False,
}

# An orphan transport reaction — no gene rule.
_RXN_ORPHAN_TRANSPORT = {
    'id': 'ORPHAN_NA_TRANS',
    'gene_reaction_rule': '',
    'metabolites': {'na_c': -1, 'na_e': 1},
    'lower_bound': 0,
    'upper_bound': 1000,
    'reversible': False,
}

# A high-degree cofactor reaction (used to test cofactor filtering).
_RXN_COFACTOR = {
    'id': 'COFACTOR_TRANS',
    'gene_reaction_rule': '9999_AT1',
    'metabolites': {'ubiquitous_c': -1, 'ubiquitous_e': 1},
    'lower_bound': 0,
    'upper_bound': 1000,
    'reversible': False,
}


def _run(reactions, **kwargs):
    """Run recon3d_transporter_interactions with mocked pypath data."""

    with patch(
        'pypath.inputs.recon3d._gem.recon3d_reactions',
        return_value=reactions,
    ):
        return list(recon3d_transporter_interactions(**kwargs))


# ---------------------------------------------------------------------------
# _parse_gene_rule
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
        result = _parse_gene_rule('1234_AT2')
        assert result == ['1234']

    def test_or_rule_gives_multiple_enzymes(self):
        result = _parse_gene_rule('1234_AT1 or 5678_AT1')
        assert len(result) == 2
        assert '1234' in result
        assert '5678' in result

    def test_and_rule_gives_complex(self):
        result = _parse_gene_rule('1234_AT1 and 5678_AT1')
        assert len(result) == 1
        # complex subunits joined with '_', sorted
        assert result[0] == '1234_5678'

    def test_parentheses_stripped(self):
        result = _parse_gene_rule('(1234_AT1 or 5678_AT1)')
        assert len(result) == 2

    def test_mixed_or_and(self):
        result = _parse_gene_rule('1111_AT1 or (2222_AT1 and 3333_AT1)')
        assert '1111' in result
        assert '2222_3333' in result

    def test_complex_subunits_sorted(self):
        # sorted order: '2222' before '9999'
        result = _parse_gene_rule('9999_AT1 and 2222_AT1')
        assert result == ['2222_9999']


# ---------------------------------------------------------------------------
# _strip_compartment
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
        # 'glc_e1' — last segment is 'e1', not single letter
        base, comp = _strip_compartment('glc_e1')
        assert base == 'glc_e1'
        assert comp == ''

    def test_compound_id(self):
        # IDs like 'M00001_c' — only the trailing compartment is stripped
        base, comp = _strip_compartment('M00001_c')
        assert base == 'M00001'
        assert comp == 'c'


# ---------------------------------------------------------------------------
# Transport detection
# ---------------------------------------------------------------------------

class TestTransportDetection:

    def test_transport_reaction_yields_edges(self):
        recs = _run([_RXN_ATP_TRANSPORT])
        assert len(recs) > 0

    def test_intracell_reaction_yields_no_edges(self):
        # Same compartment on both sides → not a transport event.
        recs = _run([_RXN_INTRACELL])
        assert recs == []

    def test_transport_yields_two_edges_per_enzyme(self):
        # One enzyme, one transported metabolite → 2 edges (met→enz, enz→met).
        recs = _run([_RXN_ATP_TRANSPORT])
        assert len(recs) == 2

    def test_edge_directions(self):
        recs = _run([_RXN_ATP_TRANSPORT])
        sources = {r.source for r in recs}
        targets = {r.target for r in recs}
        # metabolite should appear as both source and target
        assert 'atp' in sources or 'atp' in targets

    def test_met_to_enzyme_edge(self):
        recs = _run([_RXN_ATP_TRANSPORT])
        met_to_enz = [r for r in recs if r.source_type == 'small_molecule']
        assert len(met_to_enz) == 1
        assert met_to_enz[0].id_type_a == 'bigg'
        assert met_to_enz[0].id_type_b == 'entrez'
        assert met_to_enz[0].source == 'atp'
        assert met_to_enz[0].target == '1234'

    def test_enzyme_to_met_edge(self):
        recs = _run([_RXN_ATP_TRANSPORT])
        enz_to_met = [r for r in recs if r.source_type == 'protein']
        assert len(enz_to_met) == 1
        assert enz_to_met[0].id_type_a == 'entrez'
        assert enz_to_met[0].id_type_b == 'bigg'
        assert enz_to_met[0].source == '1234'
        assert enz_to_met[0].target == 'atp'

    def test_interaction_type_is_transport(self):
        recs = _run([_RXN_ATP_TRANSPORT])
        assert all(r.interaction_type == 'transport' for r in recs)

    def test_resource_is_recon3d(self):
        recs = _run([_RXN_ATP_TRANSPORT])
        assert all(r.resource == 'Recon3D' for r in recs)

    def test_transport_from_to_in_attrs(self):
        recs = _run([_RXN_ATP_TRANSPORT])
        for r in recs:
            assert 'transport_from' in r.attrs
            assert 'transport_to' in r.attrs

    def test_compartment_in_locations(self):
        recs = _run([_RXN_ATP_TRANSPORT])
        # Each edge should carry the metabolite compartment in locations
        for r in recs:
            assert isinstance(r.locations, tuple)
            assert len(r.locations) == 1

    def test_enzyme_complex_detected(self):
        recs = _run([_RXN_GLC_TRANSPORT_REV], include_reverse=False)
        enz_ids = {r.source if r.source_type == 'protein' else r.target for r in recs}
        assert any('_' in e for e in enz_ids), 'AND-rule complex not joined'
        complex_edges = [r for r in recs if r.attrs.get('enzyme_complex')]
        assert len(complex_edges) > 0


# ---------------------------------------------------------------------------
# Reverse edges
# ---------------------------------------------------------------------------

class TestReverseEdges:

    def test_reversible_reaction_yields_double_edges(self):
        recs = _run([_RXN_GLC_TRANSPORT_REV], include_reverse=True)
        forward = [r for r in recs if not r.attrs['reverse']]
        reverse = [r for r in recs if r.attrs['reverse']]
        assert len(forward) == 2
        assert len(reverse) == 2

    def test_include_reverse_false_suppresses_reverse(self):
        recs = _run([_RXN_GLC_TRANSPORT_REV], include_reverse=False)
        assert all(not r.attrs['reverse'] for r in recs)

    def test_non_reversible_has_no_reverse_edges(self):
        recs = _run([_RXN_ATP_TRANSPORT], include_reverse=True)
        assert all(not r.attrs['reverse'] for r in recs)

    def test_reverse_swaps_transport_from_to(self):
        recs = _run([_RXN_GLC_TRANSPORT_REV], include_reverse=True)
        fwd = next(r for r in recs if not r.attrs['reverse'])
        rev = next(r for r in recs if r.attrs['reverse'])
        assert fwd.attrs['transport_from'] == rev.attrs['transport_to']
        assert fwd.attrs['transport_to'] == rev.attrs['transport_from']


# ---------------------------------------------------------------------------
# Orphan reactions
# ---------------------------------------------------------------------------

class TestOrphanReactions:

    def test_orphan_included_by_default(self):
        recs = _run([_RXN_ORPHAN_TRANSPORT])
        assert len(recs) > 0

    def test_orphan_excluded_when_flag_false(self):
        recs = _run([_RXN_ORPHAN_TRANSPORT], include_orphans=False)
        assert recs == []

    def test_orphan_attrs_flag(self):
        recs = _run([_RXN_ORPHAN_TRANSPORT])
        assert all(r.attrs.get('orphan') for r in recs)

    def test_orphan_uses_reaction_id_as_enzyme(self):
        recs = _run([_RXN_ORPHAN_TRANSPORT])
        protein_sides = [
            r.target if r.source_type == 'small_molecule' else r.source
            for r in recs
        ]
        assert all(p == 'ORPHAN_NA_TRANS' for p in protein_sides)

    def test_orphan_id_type_is_reaction_id(self):
        recs = _run([_RXN_ORPHAN_TRANSPORT])
        for r in recs:
            if r.source_type == 'small_molecule':
                assert r.id_type_b == 'reaction_id'
            else:
                assert r.id_type_a == 'reaction_id'

    def test_normal_edges_have_no_orphan_flag(self):
        recs = _run([_RXN_ATP_TRANSPORT])
        assert all(not r.attrs.get('orphan') for r in recs)


# ---------------------------------------------------------------------------
# Cofactor filtering
# ---------------------------------------------------------------------------

class TestCofactorFilter:

    def _make_many_transport_reactions(self, n, met_id='ubiq'):
        """Create n transport reactions all sharing the same metabolite."""

        return [
            {
                'id': f'TRANS_{i}',
                'gene_reaction_rule': f'{1000 + i}_AT1',
                'metabolites': {f'{met_id}_c': -1, f'{met_id}_e': 1},
                'lower_bound': 0,
                'upper_bound': 1000,
                'reversible': False,
            }
            for i in range(n)
        ]

    def test_high_degree_metabolite_filtered(self):
        # 10 reactions each contributing 2 edges → degree = 20 per metabolite
        reactions = self._make_many_transport_reactions(10, 'ubiq')
        recs = _run(reactions, metab_max_degree=5)
        assert recs == []

    def test_low_degree_metabolite_kept(self):
        reactions = self._make_many_transport_reactions(2, 'rare')
        recs = _run(reactions, metab_max_degree=100)
        assert len(recs) > 0

    def test_threshold_is_exclusive(self):
        # degree == threshold should be kept; degree > threshold dropped
        reactions = self._make_many_transport_reactions(3, 'borderline')
        # degree = 6 (3 reactions × 2 edges each)
        assert _run(reactions, metab_max_degree=6) != []
        assert _run(reactions, metab_max_degree=5) == []
