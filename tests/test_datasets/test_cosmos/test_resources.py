#!/usr/bin/env python

"""Tests for individual COSMOS resource processors.

All tests in this module are marked slow because they download data
from external databases.  Run with ``pytest -m slow`` to include them.
"""

import pytest

from omnipath_metabo.datasets.cosmos._record import Interaction

pytestmark = pytest.mark.slow


def _validate_interactions(interactions, resource_name, min_count=1):
    """Shared assertions for any list of Interaction records."""

    assert len(interactions) >= min_count, (
        f'{resource_name}: expected >= {min_count} records, '
        f'got {len(interactions)}'
    )

    for rec in interactions[:10]:
        assert isinstance(rec, Interaction)
        assert rec.resource == resource_name
        assert rec.source_type in ('small_molecule', 'protein')
        assert rec.target_type in ('small_molecule', 'protein')
        assert rec.mor in (1, -1, 0)
        assert isinstance(rec.locations, tuple)


class TestStitch:

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            stitch_interactions,
        )

        records = list(stitch_interactions(
            organism=9606,
            score_threshold=900,
        ))

        _validate_interactions(records, 'STITCH')

    def test_id_types(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            stitch_interactions,
        )

        rec = next(stitch_interactions(organism=9606, score_threshold=900))
        assert rec.id_type_a == 'pubchem'
        assert rec.id_type_b == 'ensp'

    def test_mode_filter(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            stitch_interactions,
        )

        records = list(stitch_interactions(
            organism=9606,
            score_threshold=900,
            mode='activation',
        ))

        assert all(r.mor == 1 for r in records)


class TestSlc:

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import slc_interactions

        records = list(slc_interactions())
        _validate_interactions(records, 'SLC')

    def test_has_locations(self):
        from omnipath_metabo.datasets.cosmos.resources import slc_interactions

        records = list(slc_interactions())
        with_locs = [r for r in records if r.locations]
        assert len(with_locs) > 0

    def test_locations_are_tuples(self):
        from omnipath_metabo.datasets.cosmos.resources import slc_interactions

        rec = next(slc_interactions())
        assert isinstance(rec.locations, tuple)
        assert all(isinstance(c, str) for c in rec.locations)

    def test_mor_is_one(self):
        from omnipath_metabo.datasets.cosmos.resources import slc_interactions

        records = list(slc_interactions())
        assert all(r.mor == 1 for r in records)

    def test_non_human_yields_nothing(self):
        from omnipath_metabo.datasets.cosmos.resources import slc_interactions

        records = list(slc_interactions(organism=10090))
        assert records == []


class TestBrenda:

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            brenda_regulations,
        )

        records = list(brenda_regulations())
        _validate_interactions(records, 'BRENDA', min_count=100)

    def test_has_activators_and_inhibitors(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            brenda_regulations,
        )

        records = list(brenda_regulations())
        mors = {r.mor for r in records}
        assert 1 in mors
        assert -1 in mors

    def test_id_types(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            brenda_regulations,
        )

        rec = next(brenda_regulations())
        assert rec.id_type_a == 'synonym'
        assert rec.id_type_b in ('uniprot', 'genesymbol')


class TestTcdb:

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            tcdb_interactions,
        )

        records = list(tcdb_interactions())
        _validate_interactions(records, 'TCDB')

    def test_has_locations(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            tcdb_interactions,
        )

        records = list(tcdb_interactions())
        with_locs = [r for r in records if r.locations]
        assert len(with_locs) > 0


class TestMrclinksdb:

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            mrclinksdb_interactions,
        )

        records = list(mrclinksdb_interactions())
        _validate_interactions(records, 'MRCLinksDB')

    def test_has_locations(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            mrclinksdb_interactions,
        )

        records = list(mrclinksdb_interactions())
        with_locs = [r for r in records if r.locations]
        assert len(with_locs) > 0


class TestStitchNewFeatures:
    """Integration tests for STITCH protein-role classification and filters."""

    def test_interaction_type_values(self):
        from omnipath_metabo.datasets.cosmos.resources import stitch_interactions

        records = list(stitch_interactions(organism=9606, score_threshold=900))
        valid_types = {'receptor', 'transporter', 'other'}
        for rec in records[:50]:
            assert rec.interaction_type in valid_types, (
                f'Unexpected interaction_type: {rec.interaction_type!r}'
            )

    def test_stitch_mode_in_attrs(self):
        from omnipath_metabo.datasets.cosmos.resources import stitch_interactions

        rec = next(stitch_interactions(organism=9606, score_threshold=900))
        assert 'stitch_mode' in rec.attrs
        assert rec.attrs['stitch_mode'] in (
            'activation', 'inhibition', 'binding',
            'pred_bind', 'expression', 'reaction', 'catalysis',
        )

    def test_a_is_acting_filter_reduces_records(self):
        from omnipath_metabo.datasets.cosmos.resources import stitch_interactions

        # Collect first 200 interactions with and without the filter.
        with_filter, without_filter = [], []

        for i, rec in enumerate(
            stitch_interactions(organism=9606, score_threshold=900, a_is_acting=True)
        ):
            with_filter.append(rec)
            if i >= 199:
                break

        for i, rec in enumerate(
            stitch_interactions(organism=9606, score_threshold=900, a_is_acting=False)
        ):
            without_filter.append(rec)
            if i >= 199:
                break

        # The directed filter should change the mode distribution.
        from collections import Counter
        modes_with = Counter(r.attrs['stitch_mode'] for r in with_filter)
        modes_without = Counter(r.attrs['stitch_mode'] for r in without_filter)
        # Without filter, binding dominates; with filter, it is reduced.
        assert modes_without.get('binding', 0) > modes_with.get('binding', 0)

    def test_binding_mode_included_in_default(self):
        from omnipath_metabo.datasets.cosmos.resources import stitch_interactions

        modes = set()
        for i, rec in enumerate(
            stitch_interactions(organism=9606, score_threshold=900, a_is_acting=False)
        ):
            modes.add(rec.attrs['stitch_mode'])
            if i >= 499:
                break

        assert 'binding' in modes


class TestGem:
    """Integration tests for the GEM processor (downloads Human-GEM)."""

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        records = list(gem_interactions(gem='Human-GEM'))
        assert len(records) > 1000

    def test_has_normal_and_orphan_edges(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        records = list(gem_interactions(gem='Human-GEM'))
        normal = [r for r in records if not r.attrs.get('orphan')]
        orphans = [r for r in records if r.attrs.get('orphan')]
        assert len(normal) > 0
        assert len(orphans) > 0

    def test_orphan_id_types(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        records = list(gem_interactions(gem='Human-GEM'))
        for rec in records:
            if rec.attrs.get('orphan'):
                if rec.source_type == 'small_molecule':
                    assert rec.id_type_b == 'reaction_id'
                else:
                    assert rec.id_type_a == 'reaction_id'

    def test_normal_edge_id_types(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        records = list(gem_interactions(gem='Human-GEM'))
        for rec in records[:50]:
            if not rec.attrs.get('orphan'):
                if rec.source_type == 'small_molecule':
                    assert rec.id_type_a == 'metatlas'
                    assert rec.id_type_b == 'ensembl'
                else:
                    assert rec.id_type_a == 'ensembl'
                    assert rec.id_type_b == 'metatlas'

    def test_transport_resource_label(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        records = list(gem_interactions(gem='Human-GEM'))
        resources = {r.resource for r in records}
        assert 'GEM:Human-GEM' in resources
        assert 'GEM_transporter:Human-GEM' in resources

    def test_include_orphans_false(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        records = list(gem_interactions(gem='Human-GEM', include_orphans=False))
        assert all(not r.attrs.get('orphan') for r in records)

    def test_reverse_edges_present(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        records = list(gem_interactions(gem='Human-GEM', include_reverse=True))
        assert any(r.attrs.get('reverse') for r in records)

    def test_gems_provenance_single(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        records = list(gem_interactions(gem='Human-GEM'))
        for rec in records[:20]:
            assert 'gems' in rec.attrs
            assert rec.attrs['gems'] == ['Human-GEM']


@pytest.mark.slow
class TestStitchMouse:
    """Integration tests verifying STITCH yields data for mouse (organism=10090)."""

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import stitch_interactions

        records = list(stitch_interactions(organism=10090, score_threshold=900))
        assert len(records) > 0

    def test_resource_label(self):
        from omnipath_metabo.datasets.cosmos.resources import stitch_interactions

        rec = next(stitch_interactions(organism=10090, score_threshold=900))
        assert rec.resource == 'STITCH'


@pytest.mark.slow
class TestMrclinksdbMouse:
    """Integration tests verifying MRCLinksDB yields data for mouse (organism=10090)."""

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import mrclinksdb_interactions

        records = list(mrclinksdb_interactions(organism=10090))
        assert len(records) > 0


@pytest.mark.slow
class TestGemMouse:
    """Integration tests verifying GEM yields data for Mouse-GEM."""

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        records = list(gem_interactions(gem='Mouse-GEM'))
        assert len(records) > 1000

    def test_gems_provenance(self):
        from omnipath_metabo.datasets.cosmos.resources import gem_interactions

        rec = next(gem_interactions(gem='Mouse-GEM'))
        assert rec.attrs['gems'] == ['Mouse-GEM']


@pytest.mark.slow
class TestBuildMouse:
    """Integration test verifying build() accepts organism=10090."""

    def test_build_mouse_organism(self):
        import pandas as pd

        from omnipath_metabo.datasets.cosmos import build

        df = build(
            organism=10090,
            slc=False,
            recon3d=False,
            stitch={'score_threshold': 900},
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0


class TestRecon3d:
    """Integration tests for the Recon3D transporter processor."""

    def test_yields_interactions(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            recon3d_transporter_interactions,
        )

        records = list(recon3d_transporter_interactions())
        assert len(records) > 1000

    def test_has_normal_and_orphan_edges(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            recon3d_transporter_interactions,
        )

        records = list(recon3d_transporter_interactions())
        normal = [r for r in records if not r.attrs.get('orphan')]
        orphans = [r for r in records if r.attrs.get('orphan')]
        assert len(normal) > 0
        assert len(orphans) > 0

    def test_interaction_type_is_transport(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            recon3d_transporter_interactions,
        )

        records = list(recon3d_transporter_interactions())
        assert all(r.interaction_type == 'transport' for r in records)

    def test_normal_edge_id_types(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            recon3d_transporter_interactions,
        )

        records = list(recon3d_transporter_interactions())
        for rec in records[:50]:
            if not rec.attrs.get('orphan'):
                if rec.source_type == 'small_molecule':
                    assert rec.id_type_a == 'bigg'
                    assert rec.id_type_b == 'entrez'
                else:
                    assert rec.id_type_a == 'entrez'
                    assert rec.id_type_b == 'bigg'

    def test_include_orphans_false(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            recon3d_transporter_interactions,
        )

        records = list(recon3d_transporter_interactions(include_orphans=False))
        assert all(not r.attrs.get('orphan') for r in records)

    def test_transport_attrs_present(self):
        from omnipath_metabo.datasets.cosmos.resources import (
            recon3d_transporter_interactions,
        )

        rec = next(recon3d_transporter_interactions())
        assert 'transport_from' in rec.attrs
        assert 'transport_to' in rec.attrs
        assert 'reaction_id' in rec.attrs
