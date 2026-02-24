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


def _metabolite(mid='MAM001c', compartment='c'):
    return GemMetabolite(id=mid, name=mid, compartment=compartment,
                         formula='', charge=0)


def _interaction(src='MAM001c', tgt='ENSG001',
                 src_type='metabolite', tgt_type='gene',
                 src_comp='c', tgt_comp='',
                 rxn_id='MAR001', reverse=False):
    return GemInteraction(
        source=src, target=tgt,
        source_type=src_type, target_type=tgt_type,
        source_compartment=src_comp, target_compartment=tgt_comp,
        reaction_id=rxn_id, reverse=reverse,
    )


# ---------------------------------------------------------------------------
# Mock helpers
# ---------------------------------------------------------------------------

_METS = [_metabolite('MAM001c', 'c'), _metabolite('MAM002c', 'c')]
_METS_LOW = [_metabolite('MAM003c', 'c'), _metabolite('MAM004c', 'c')]


def _patch_gem(gem_name, interactions, reactions=None, metabolites=None):
    """Return a context manager that patches all three metatlas functions."""

    reactions = reactions or [_reaction()]
    metabolites = metabolites or _METS

    def _network(**kwargs):
        yield from interactions

    def _yaml_reactions(**kwargs):
        yield from reactions

    def _yaml_metabolites(**kwargs):
        yield from metabolites

    return (
        patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
              side_effect=_network),
        patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
              side_effect=_yaml_reactions),
        patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_metabolites',
              side_effect=_yaml_metabolites),
    )


# ---------------------------------------------------------------------------
# Tests: provenance (attrs['gems'])
# ---------------------------------------------------------------------------

class TestGemProvenance:

    def test_gems_attr_present(self):
        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        ints = [_interaction()]
        with _patch_gem('TestGEM', ints)[0], \
             _patch_gem('TestGEM', ints)[1], \
             _patch_gem('TestGEM', ints)[2]:
            records = list(gem_interactions(gem='TestGEM'))

        for rec in records:
            assert 'gems' in rec.attrs

    def test_gems_attr_single_gem_value(self):
        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        ints = [_interaction()]
        p0, p1, p2 = _patch_gem('MyGEM', ints)
        with p0, p1, p2:
            records = list(gem_interactions(gem='MyGEM'))

        assert len(records) > 0
        for rec in records:
            assert rec.attrs['gems'] == ['MyGEM']

    def test_gems_attr_is_list(self):
        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        ints = [_interaction()]
        p0, p1, p2 = _patch_gem('G', ints)
        with p0, p1, p2:
            records = list(gem_interactions(gem='G'))

        for rec in records:
            assert isinstance(rec.attrs['gems'], list)


# ---------------------------------------------------------------------------
# Tests: deduplication across GEMs
# ---------------------------------------------------------------------------

class TestGemDeduplication:

    def _run_two_gems(self, ints_a, ints_b, rxns_a=None, rxns_b=None,
                     mets_a=None, mets_b=None):
        """Run gem_interactions with two mocked GEMs, return records."""

        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        rxns_a = rxns_a or [_reaction()]
        rxns_b = rxns_b or [_reaction()]
        mets_a = mets_a or _METS
        mets_b = mets_b or _METS

        call_count = [0]

        def _network(**kwargs):
            i = call_count[0]
            call_count[0] += 1
            yield from (ints_a if i == 0 else ints_b)

        def _yaml_reactions(**kwargs):
            yield from rxns_a  # transport IDs don't differ between gems in tests

        def _yaml_metabolites(**kwargs):
            yield from mets_a

        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_metabolites',
                  side_effect=_yaml_metabolites),
        ):
            return list(gem_interactions(gem=['GemA', 'GemB']))

    def test_identical_edges_collapsed(self):
        """Same (source, target, reaction_id, reverse) from two GEMs → one edge."""

        edge = _interaction()
        records = self._run_two_gems([edge], [edge])
        # met→enzyme and enzyme→met: each direction should be deduplicated
        sources = [r.source for r in records]
        # No duplicates in (source, target, reaction_id, reverse) space
        keys = [
            (r.source, r.target,
             r.attrs.get('reaction_id'), r.attrs.get('reverse'))
            for r in records
        ]
        assert len(keys) == len(set(keys))

    def test_identical_edges_provenance(self):
        """Deduplicated edge carries both GEM names."""

        edge = _interaction()
        records = self._run_two_gems([edge], [edge])

        for rec in records:
            if rec.source == 'MAM001':  # metabolite→enzyme edge (compartment stripped)
                assert rec.attrs['gems'] == ['GemA', 'GemB']

    def test_distinct_edges_not_collapsed(self):
        """Different metabolites → both edges preserved."""

        edge_a = _interaction(src='MAM001c', tgt='ENSG001', rxn_id='MAR001',
                              src_comp='c')
        edge_b = _interaction(src='MAM003c', tgt='ENSG002', rxn_id='MAR002',
                              src_comp='c')

        mets_b = [_metabolite('MAM003c', 'c'), _metabolite('MAM004c', 'c')]

        call_count = [0]
        rxns = [_reaction()]
        all_mets = [_METS, mets_b]

        def _network(**kwargs):
            i = call_count[0]
            call_count[0] += 1
            yield from ([edge_a] if i == 0 else [edge_b])

        def _yaml_reactions(**kwargs):
            yield from rxns

        met_call = [0]

        def _yaml_metabolites(**kwargs):
            i = met_call[0]
            met_call[0] += 1
            yield from all_mets[min(i, 1)]

        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions
        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_metabolites',
                  side_effect=_yaml_metabolites),
        ):
            records = list(gem_interactions(gem=['GemA', 'GemB']))

        # met→enzyme edges for MAM001 and MAM003 should both be present
        met_sources = {r.source for r in records if r.source_type == 'small_molecule'}
        assert 'MAM001' in met_sources
        assert 'MAM003' in met_sources

    def test_distinct_edge_each_has_one_gem(self):
        """Non-overlapping edges each carry exactly their own GEM name."""

        edge_a = _interaction(src='MAM001c', tgt='ENSG001', rxn_id='MAR001',
                              src_comp='c')
        edge_b = _interaction(src='MAM001c', tgt='ENSG001', rxn_id='MAR002',
                              src_comp='c')  # same met+enzyme, different reaction

        call_count = [0]

        def _network(**kwargs):
            i = call_count[0]
            call_count[0] += 1
            yield from ([edge_a] if i == 0 else [edge_b])

        def _yaml_reactions(**kwargs):
            yield from [_reaction()]

        def _yaml_metabolites(**kwargs):
            yield from _METS

        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions
        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_metabolites',
                  side_effect=_yaml_metabolites),
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

        edge = _interaction()
        records = self._run_two_gems([edge], [edge])
        for rec in records:
            assert rec.attrs['gems'] == sorted(rec.attrs['gems'])


# ---------------------------------------------------------------------------
# Tests: orphan reactions (gap from previous session)
# ---------------------------------------------------------------------------

class TestGemOrphans:

    def _run(self, rxns, mets, include_orphans=True):
        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions

        def _network(**kwargs):
            return iter([])  # no normal edges

        def _yaml_reactions(**kwargs):
            yield from rxns

        def _yaml_metabolites(**kwargs):
            yield from mets

        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_metabolites',
                  side_effect=_yaml_metabolites),
        ):
            return list(gem_interactions(
                gem='TestGEM',
                include_orphans=include_orphans,
            ))

    def test_orphan_included_by_default(self):
        orphan = _reaction(rid='MAR_O1', gene_reaction_rule='',
                           mets={'MAM001c': -1, 'MAM002c': 1})
        mets = [_metabolite('MAM001c', 'c'), _metabolite('MAM002c', 'c')]
        records = self._run([orphan], mets)
        assert len(records) > 0

    def test_orphan_excluded_when_flag_false(self):
        orphan = _reaction(rid='MAR_O1', gene_reaction_rule='',
                           mets={'MAM001c': -1, 'MAM002c': 1})
        mets = [_metabolite('MAM001c', 'c'), _metabolite('MAM002c', 'c')]
        records = self._run([orphan], mets, include_orphans=False)
        assert records == []

    def test_orphan_attrs_flag(self):
        orphan = _reaction(rid='MAR_O1', gene_reaction_rule='',
                           mets={'MAM001c': -1, 'MAM002c': 1})
        mets = [_metabolite('MAM001c', 'c'), _metabolite('MAM002c', 'c')]
        records = self._run([orphan], mets)
        assert all(r.attrs.get('orphan') for r in records)

    def test_orphan_uses_reaction_id_as_enzyme(self):
        orphan = _reaction(rid='MAR_O1', gene_reaction_rule='',
                           mets={'MAM001c': -1, 'MAM002c': 1})
        mets = [_metabolite('MAM001c', 'c'), _metabolite('MAM002c', 'c')]
        records = self._run([orphan], mets)
        enzyme_nodes = {
            r.target if r.source_type == 'small_molecule' else r.source
            for r in records
        }
        assert 'MAR_O1' in enzyme_nodes

    def test_orphan_id_type_is_reaction_id(self):
        orphan = _reaction(rid='MAR_O1', gene_reaction_rule='',
                           mets={'MAM001c': -1, 'MAM002c': 1})
        mets = [_metabolite('MAM001c', 'c'), _metabolite('MAM002c', 'c')]
        records = self._run([orphan], mets)
        for rec in records:
            if rec.source_type == 'small_molecule':
                assert rec.id_type_b == 'reaction_id'
            else:
                assert rec.id_type_a == 'reaction_id'

    def test_orphan_has_gems_provenance(self):
        orphan = _reaction(rid='MAR_O1', gene_reaction_rule='',
                           mets={'MAM001c': -1, 'MAM002c': 1})
        mets = [_metabolite('MAM001c', 'c'), _metabolite('MAM002c', 'c')]
        records = self._run([orphan], mets)
        for rec in records:
            assert rec.attrs['gems'] == ['TestGEM']

    def test_non_orphan_has_no_orphan_flag(self):
        normal = _reaction(rid='MAR001', gene_reaction_rule='ENSG001',
                           mets={'MAM001c': -1, 'MAM002c': 1})
        mets = [_metabolite('MAM001c', 'c'), _metabolite('MAM002c', 'c')]

        def _network(**kwargs):
            yield _interaction(src='MAM001c', tgt='ENSG001', rxn_id='MAR001',
                               src_comp='c')

        def _yaml_reactions(**kwargs):
            yield normal

        def _yaml_metabolites(**kwargs):
            yield from mets

        from omnipath_metabo.datasets.cosmos.resources.gem import gem_interactions
        with (
            patch('pypath.inputs.metatlas._gem.metatlas_gem_network',
                  side_effect=_network),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_reactions',
                  side_effect=_yaml_reactions),
            patch('pypath.inputs.metatlas._gem.metatlas_gem_yaml_metabolites',
                  side_effect=_yaml_metabolites),
        ):
            records = list(gem_interactions(gem='TestGEM'))

        assert all(not r.attrs.get('orphan') for r in records)
