#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos._build module."""

import copy
from unittest.mock import patch

import pandas as pd
import pytest

from omnipath_metabo.datasets.cosmos._build import PROCESSORS, build
from omnipath_metabo.datasets.cosmos._record import Interaction


def _fake_resource(**kwargs):
    """A minimal fake resource generator for testing."""

    yield Interaction(
        source='CID001',
        target='P00001',
        source_type='small_molecule',
        target_type='protein',
        id_type_a='pubchem',
        id_type_b='uniprot',
        interaction_type='test',
        resource='FakeDB',
        mor=1,
    )
    yield Interaction(
        source='CID002',
        target='P00002',
        source_type='small_molecule',
        target_type='protein',
        id_type_a='pubchem',
        id_type_b='uniprot',
        interaction_type='test',
        resource='FakeDB',
        mor=-1,
    )


FAKE_PROCESSORS = {
    'fake_a': _fake_resource,
    'fake_b': _fake_resource,
}

FAKE_DEFAULT_CONFIG = {
    'organism': 9606,
    'translate_ids': False,  # fake IDs can't be translated; skip for unit tests
    'resources': {
        'fake_a': {},
        'fake_b': {},
    },
}


@pytest.fixture
def _mock_build():
    """Patch PROCESSORS and default_config for isolated build tests."""

    with (
        patch.dict(
            'omnipath_metabo.datasets.cosmos._build.PROCESSORS',
            FAKE_PROCESSORS,
            clear=True,
        ),
        patch(
            'omnipath_metabo.datasets.cosmos._config.default_config',
            side_effect=lambda: copy.deepcopy(FAKE_DEFAULT_CONFIG),
        ),
        patch(
            'omnipath_metabo.datasets.cosmos._build.config',
            wraps=None,
        ) as mock_config,
    ):
        from omnipath_metabo.datasets.cosmos._config import config
        mock_config.side_effect = config
        yield


@pytest.fixture
def _mock_processors():
    """Patch only PROCESSORS with a capture function."""

    captured = {}

    def capture_resource(**kwargs):
        captured.update(kwargs)
        yield from _fake_resource()

    processors = {'capture': capture_resource}
    default_cfg = {
        'organism': 9606,
        'resources': {'capture': {}},
    }

    with (
        patch.dict(
            'omnipath_metabo.datasets.cosmos._build.PROCESSORS',
            processors,
            clear=True,
        ),
        patch(
            'omnipath_metabo.datasets.cosmos._config.default_config',
            side_effect=lambda: copy.deepcopy(default_cfg),
        ),
        patch(
            'omnipath_metabo.datasets.cosmos._build.config',
            wraps=None,
        ) as mock_config,
    ):
        from omnipath_metabo.datasets.cosmos._config import config
        mock_config.side_effect = config
        yield captured


class TestBuildWithMocks:
    """Tests for build() using fake processors."""

    def test_returns_dataframe(self, _mock_build):
        df = build()
        assert isinstance(df, pd.DataFrame)

    def test_columns_match_interaction_fields(self, _mock_build):
        df = build()
        assert list(df.columns) == list(Interaction._fields)

    def test_default_runs_all(self, _mock_build):
        df = build()
        assert len(df) == 4  # 2 records x 2 fake resources

    def test_single_resource(self, _mock_build):
        df = build(fake_b=False)
        assert len(df) == 2

    def test_disable_resource(self, _mock_build):
        df = build(fake_a=False, fake_b=False)
        assert len(df) == 0

    def test_unknown_resource_raises(self, _mock_build):
        with pytest.raises(ValueError, match='Unknown resource'):
            build(resources={'nonexistent': {}})

    def test_organism_injected(self, _mock_processors):
        build(organism=10090)
        assert _mock_processors.get('organism') == 10090

    def test_per_resource_organism_overrides(self, _mock_processors):
        build(
            organism=9606,
            resources={'capture': {'organism': 10090}},
        )
        assert _mock_processors.get('organism') == 10090

    def test_organism_injected_arbitrary_value(self, _mock_processors):
        """organism= is forwarded for any NCBI taxon ID, not just 9606/10090."""

        build(organism=10116)  # rat
        assert _mock_processors.get('organism') == 10116

    def test_resource_params_passed(self, _mock_processors):
        build(resources={'capture': {'score_threshold': 500}})
        assert _mock_processors.get('score_threshold') == 500

    def test_empty_resources_override_is_noop(self, _mock_build):
        """Empty resources override doesn't change defaults."""

        df = build(resources={})
        assert len(df) == 4  # defaults still active


class TestProcessorsRegistry:
    """Tests for the PROCESSORS constant."""

    def test_all_resources_registered(self):
        expected = {'stitch', 'tcdb', 'slc', 'brenda', 'mrclinksdb', 'gem', 'recon3d'}
        assert set(PROCESSORS) == expected

    def test_all_callable(self):
        for name, func in PROCESSORS.items():
            assert callable(func), f'{name} is not callable'


@pytest.mark.slow
class TestBuildIntegration:
    """Integration tests that call real resource processors.

    These tests download data from external databases and may take
    several minutes.  Run with ``pytest -m slow`` to include them.
    """

    def test_build_stitch_only(self):
        df = build(
            stitch={'score_threshold': 900},
            tcdb=False, slc=False, brenda=False, mrclinksdb=False,
        )

        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert (df['resource'] == 'STITCH').all()
        assert set(df['source_type']) == {'small_molecule'}
        assert set(df['target_type']) == {'protein'}

    def test_build_slc_only(self):
        df = build(
            slc={},
            stitch=False, tcdb=False, brenda=False, mrclinksdb=False,
        )

        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert (df['resource'] == 'SLC').all()

    def test_build_brenda_only(self):
        df = build(
            brenda={},
            stitch=False, tcdb=False, slc=False, mrclinksdb=False,
        )

        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert (df['resource'] == 'BRENDA').all()
