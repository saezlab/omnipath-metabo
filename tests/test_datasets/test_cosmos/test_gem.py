#!/usr/bin/env python

"""Fast unit tests for omnipath_metabo.datasets.cosmos.resources.gem."""

from unittest.mock import patch

import pytest

from pypath.inputs.metatlas._records import GemInteraction, GemMetabolite, GemReaction


# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------

def _reaction(
    rid='MAR001',
    subsystem='Metabolism',
    gene_reaction_rule='ENSG001',
    mets=None,
    lb=-1000.0,
    ub=1000.0,
):
    """Build a minimal GemReaction."""

    return GemReaction(
        id=rid,
        name=rid,
        metabolites=mets or {'MAM001c': -1, 'MAM002c': 1},
        lower_bound=lb,
        upper_bound=ub,
        gene_reaction_rule=gene_reaction_rule,
        subsystem=subsystem,
        eccodes=(),
    )


def _gi(src='MAM001c', tgt='ENSG001',
        src_type='metabolite', tgt_type='protein',
        src_comp='c', tgt_comp='',
        rxn_id='MAR001', reverse=False):
    """Build a GemInteraction record."""

    return GemInteraction(
        source=src, target=tgt,
        source_type=src_type, target_type=tgt_type,
        source_compartment=src_comp, target_compartment=tgt_comp,
        reaction_id=rxn_id, reverse=reverse,
    )


# Keep old alias for backward compat in tests that still use _interaction.
_interaction = _gi


# ---------------------------------------------------------------------------
# Mock helpers
# ---------------------------------------------------------------------------

def _patch_gem(gem_name, metabolic_ints, transport_ints=None, reactions=None):
    """Return context managers patching the three pypath functions used by gem_interactions.

    - metatlas_gem_network       → yields metabolic GemInteraction records
    - metatlas_gem_transport_network → yields transport GemInteraction records
    - metatlas_gem_yaml_reactions → yields GemReaction records (for transport_ids + gene_id_type)
    """

    transport_ints = transport_ints or []
    reactions = reactions or [_reaction()]

    def _network(**kwargs):
        yield from metabolic_ints

    def _transport_network(**kwargs):
        yield from transport_ints

    def _yaml_reactions(**kwargs):
        yield from reactions

    return (
        patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
              side_effect=_network),
        patch('pypath.inputs.metatlas._gem.metatlas_gem_transport_network',
              side_effect=_transport_network),
        patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
              side_effect=_yaml_reactions),
    )


# ---------------------------------------------------------------------------
# Tests: provenance (attrs['gems'])
# ---------------------------------------------------------------------------

class TestGemProvenance:

    def test_gems_attr_present(self):
        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        ints = [_gi()]
        p0, p1, p2 = _patch_gem('TestGEM', ints)
        with p0, p1, p2:
            records = list(gem_interactions(gem='TestGEM'))

        for rec in records:
            assert 'gems' in rec.attrs

    def test_gems_attr_single_gem_value(self):
        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        ints = [_gi()]
        p0, p1, p2 = _patch_gem('MyGEM', ints)
        with p0, p1, p2:
            records = list(gem_interactions(gem='MyGEM'))

        assert len(records) > 0
        for rec in records:
            assert rec.attrs['gems'] == ['MyGEM']

    def test_gems_attr_is_list(self):
        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        ints = [_gi()]
        p0, p1, p2 = _patch_gem('G', ints)
        with p0, p1, p2:
            records = list(gem_interactions(gem='G'))

        for rec in records:
            assert isinstance(rec.attrs['gems'], list)


# ---------------------------------------------------------------------------
# Tests: deduplication across GEMs
# ---------------------------------------------------------------------------

class TestGemDeduplication:

    def _run_two_gems(self, ints_a, ints_b, rxns=None):
        """Run gem_interactions with two mocked GEMs, return records."""

        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        rxns = rxns or [_reaction()]
        call_count = [0]

        def _network(**kwargs):
            i = call_count[0]
            call_count[0] += 1
            yield from (ints_a if i == 0 else ints_b)

        def _transport_network(**kwargs):
            return iter([])

        def _yaml_reactions(**kwargs):
            yield from rxns

        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_transport_network',
                  side_effect=_transport_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
        ):
            return list(gem_interactions(gem=['GemA', 'GemB']))

    def test_identical_edges_collapsed(self):
        """Same (source, target, reaction_id, reverse) from two GEMs → one edge."""

        edge = _gi()
        records = self._run_two_gems([edge], [edge])
        keys = [
            (r.source, r.target,
             r.attrs.get('reaction_id'), r.attrs.get('reverse'))
            for r in records
        ]
        assert len(keys) == len(set(keys))

    def test_identical_edges_provenance(self):
        """Deduplicated edge carries both GEM names."""

        edge = _gi()
        records = self._run_two_gems([edge], [edge])

        for rec in records:
            if rec.source == 'MAM001':  # metabolite→enzyme edge (compartment stripped)
                assert rec.attrs['gems'] == ['GemA', 'GemB']

    def test_distinct_edges_not_collapsed(self):
        """Different metabolites → both edges preserved."""

        edge_a = _gi(src='MAM001c', tgt='ENSG001', rxn_id='MAR001', src_comp='c')
        edge_b = _gi(src='MAM003c', tgt='ENSG002', rxn_id='MAR002', src_comp='c')

        call_count = [0]

        def _network(**kwargs):
            i = call_count[0]
            call_count[0] += 1
            yield from ([edge_a] if i == 0 else [edge_b])

        def _transport_network(**kwargs):
            return iter([])

        def _yaml_reactions(**kwargs):
            yield _reaction()

        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions
        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_transport_network',
                  side_effect=_transport_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
        ):
            records = list(gem_interactions(gem=['GemA', 'GemB']))

        met_sources = {r.source for r in records if r.source_type == 'small_molecule'}
        assert 'MAM001' in met_sources
        assert 'MAM003' in met_sources

    def test_distinct_edge_each_has_one_gem(self):
        """Non-overlapping edges each carry exactly their own GEM name."""

        edge_a = _gi(src='MAM001c', tgt='ENSG001', rxn_id='MAR001', src_comp='c')
        edge_b = _gi(src='MAM001c', tgt='ENSG001', rxn_id='MAR002', src_comp='c')

        call_count = [0]

        def _network(**kwargs):
            i = call_count[0]
            call_count[0] += 1
            yield from ([edge_a] if i == 0 else [edge_b])

        def _transport_network(**kwargs):
            return iter([])

        def _yaml_reactions(**kwargs):
            yield _reaction()

        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions
        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_transport_network',
                  side_effect=_transport_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
        ):
            records = list(gem_interactions(gem=['GemA', 'GemB']))

        met_edges = [r for r in records if r.source_type == 'small_molecule']
        rxn_ids = {r.attrs['reaction_id'] for r in met_edges}
        assert 'MAR001' in rxn_ids
        assert 'MAR002' in rxn_ids
        for rec in met_edges:
            assert len(rec.attrs['gems']) == 1

    def test_gem_names_sorted(self):
        """attrs['gems'] is always sorted."""

        edge = _gi()
        records = self._run_two_gems([edge], [edge])
        for rec in records:
            assert rec.attrs['gems'] == sorted(rec.attrs['gems'])


# ---------------------------------------------------------------------------
# Tests: orphan reactions
# ---------------------------------------------------------------------------

class TestGemOrphans:

    def _run(self, orphan_gis, include_orphans=True):
        """Run gem_interactions with orphan GemInteraction records from metatlas_gem_network."""

        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        def _network(**kwargs):
            yield from orphan_gis

        def _transport_network(**kwargs):
            return iter([])

        def _yaml_reactions(**kwargs):
            yield _reaction()

        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_transport_network',
                  side_effect=_transport_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
        ):
            return list(gem_interactions(
                gem='TestGEM',
                include_orphans=include_orphans,
            ))

    def _make_orphan_gis(self, rid='MAR_O1'):
        """Build GemInteraction records for an orphan reaction (reaction as pseudo-enzyme)."""

        return [
            _gi(src='MAM001c', tgt=rid, src_type='metabolite', tgt_type='reaction',
                src_comp='c', tgt_comp='', rxn_id=rid),
            _gi(src=rid, tgt='MAM002c', src_type='reaction', tgt_type='metabolite',
                src_comp='', tgt_comp='c', rxn_id=rid),
        ]

    def test_orphan_included_by_default(self):
        records = self._run(self._make_orphan_gis())
        assert len(records) > 0

    def test_orphan_excluded_when_flag_false(self):
        # With include_orphans=False, metatlas_gem_network would not yield orphan records.
        # We simulate this by passing empty list.
        records = self._run([], include_orphans=False)
        assert records == []

    def test_orphan_attrs_flag(self):
        records = self._run(self._make_orphan_gis())
        assert all(r.attrs.get('orphan') for r in records)

    def test_orphan_uses_reaction_id_as_enzyme(self):
        records = self._run(self._make_orphan_gis())
        enzyme_nodes = {
            r.target if r.source_type == 'small_molecule' else r.source
            for r in records
        }
        assert 'MAR_O1' in enzyme_nodes

    def test_orphan_id_type_is_reaction_id(self):
        records = self._run(self._make_orphan_gis())
        for rec in records:
            if rec.source_type == 'small_molecule':
                assert rec.id_type_b == 'reaction_id'
            else:
                assert rec.id_type_a == 'reaction_id'

    def test_orphan_has_gems_provenance(self):
        records = self._run(self._make_orphan_gis())
        for rec in records:
            assert rec.attrs['gems'] == ['TestGEM']

    def test_no_duplicate_forward_reverse_keys(self):
        """Each (source, target, reaction_id, reverse) key is unique."""

        # Forward + reverse orphan edges for a reversible reaction.
        rid = 'MAR_O1'
        gis = [
            _gi(src='MAM001c', tgt=rid, src_type='metabolite', tgt_type='reaction',
                src_comp='c', tgt_comp='', rxn_id=rid, reverse=False),
            _gi(src=rid, tgt='MAM002c', src_type='reaction', tgt_type='metabolite',
                src_comp='', tgt_comp='c', rxn_id=rid, reverse=False),
            _gi(src='MAM002c', tgt=rid, src_type='metabolite', tgt_type='reaction',
                src_comp='c', tgt_comp='', rxn_id=rid, reverse=True),
            _gi(src=rid, tgt='MAM001c', src_type='reaction', tgt_type='metabolite',
                src_comp='', tgt_comp='c', rxn_id=rid, reverse=True),
        ]
        records = self._run(gis)
        full_keys = [
            (r.source, r.target, r.attrs.get('reaction_id'), r.attrs.get('reverse'))
            for r in records
        ]
        assert len(full_keys) == len(set(full_keys))

    def test_non_orphan_has_no_orphan_flag(self):
        normal_gi = _gi(src='MAM001c', tgt='ENSG001', src_type='metabolite',
                        tgt_type='protein', src_comp='c', rxn_id='MAR001')

        def _network(**kwargs):
            yield normal_gi

        def _transport_network(**kwargs):
            return iter([])

        def _yaml_reactions(**kwargs):
            yield _reaction()

        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions
        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_transport_network',
                  side_effect=_transport_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
        ):
            records = list(gem_interactions(gem='TestGEM'))

        assert all(not r.attrs.get('orphan') for r in records)


