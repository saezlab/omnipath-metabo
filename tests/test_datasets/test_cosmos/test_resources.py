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
