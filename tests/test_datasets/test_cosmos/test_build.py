#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos._build module."""

import copy
from unittest.mock import patch

import pandas as pd
import pytest

from omnipath_metabo.datasets.cosmos._build import (
    PROCESSORS,
    build,
    build_allosteric,
    build_enzyme_metabolite,
    build_receptors,
    build_transporters,
)
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


# ---------------------------------------------------------------------------
# Shared fake rows for subset filter tests
# ---------------------------------------------------------------------------

def _row(interaction_type: str, resource: str, mor: int = 1) -> Interaction:
    return Interaction(
        source='CHEBI:1',
        target='ENSG00000001',
        source_type='small_molecule',
        target_type='protein',
        id_type_a='chebi',
        id_type_b='ensembl',
        interaction_type=interaction_type,
        resource=resource,
        mor=mor,
    )


# transport rows
_TRANSPORT_TCDB = _row('transport', 'TCDB')
_TRANSPORT_SLC = _row('transport', 'SLC')
_TRANSPORT_GEM = _row('transport', 'GEM_transporter:Human-GEM')
# caught only by resource prefix, not by interaction_type
_TRANSPORT_GEM_BY_RESOURCE = _row('gem_reaction', 'GEM_transporter:Human-GEM')
_TRANSPORT_STITCH = _row('transporter', 'STITCH')

# receptor rows
_RECEPTOR_MRC = _row('ligand_receptor', 'MRCLinksDB')
_RECEPTOR_STITCH = _row('receptor', 'STITCH')

# allosteric rows
_ALLOSTERIC_BRENDA = _row('allosteric_regulation', 'BRENDA')
_ALLOSTERIC_STITCH = _row('other', 'STITCH')

# enzyme-metabolite rows
_ENZYME_GEM = _row('reaction', 'GEM:Human-GEM')

_ALL_ROWS = [
    _TRANSPORT_TCDB, _TRANSPORT_SLC,
    _TRANSPORT_GEM, _TRANSPORT_GEM_BY_RESOURCE, _TRANSPORT_STITCH,
    _RECEPTOR_MRC, _RECEPTOR_STITCH,
    _ALLOSTERIC_BRENDA, _ALLOSTERIC_STITCH,
    _ENZYME_GEM,
]


@pytest.fixture
def _all_categories_df():
    return pd.DataFrame(_ALL_ROWS, columns=Interaction._fields)


@pytest.fixture
def _mock_build_fn(_all_categories_df):
    """Patch build() to return a controlled DataFrame for subset filter tests."""
    with patch(
        'omnipath_metabo.datasets.cosmos._build.build',
        side_effect=lambda *a, **kw: _all_categories_df.copy(),
    ) as mock:
        yield mock


# ---------------------------------------------------------------------------
# TestBuildTransporters
# ---------------------------------------------------------------------------

class TestBuildTransporters:
    """Tests for build_transporters() filtering and default-disabled resources."""

    def test_returns_dataframe(self, _mock_build_fn):
        assert isinstance(build_transporters(), pd.DataFrame)

    def test_index_reset(self, _mock_build_fn):
        df = build_transporters()
        assert list(df.index) == list(range(len(df)))

    def test_total_row_count(self, _mock_build_fn):
        # TCDB, SLC, GEM_transporter (×2: by type + by resource), STITCH transporter
        assert len(build_transporters()) == 5

    def test_keeps_transport_type(self, _mock_build_fn):
        df = build_transporters()
        assert 'transport' in df['interaction_type'].values

    def test_keeps_gem_transporter_by_resource_prefix(self, _mock_build_fn):
        df = build_transporters()
        # both GEM_transporter rows: one matched via type, one only via resource prefix
        gem_rows = df[df['resource'].str.startswith('GEM_transporter')]
        assert len(gem_rows) == 2

    def test_keeps_stitch_transporter(self, _mock_build_fn):
        df = build_transporters()
        stitch_rows = df[df['resource'].eq('STITCH')]
        assert len(stitch_rows) == 1
        assert stitch_rows.iloc[0]['interaction_type'] == 'transporter'

    def test_excludes_receptor_rows(self, _mock_build_fn):
        df = build_transporters()
        assert 'ligand_receptor' not in df['interaction_type'].values
        assert 'receptor' not in df['interaction_type'].values

    def test_excludes_allosteric_rows(self, _mock_build_fn):
        df = build_transporters()
        assert 'allosteric_regulation' not in df['interaction_type'].values

    def test_excludes_gem_metabolic(self, _mock_build_fn):
        df = build_transporters()
        # GEM: prefix (metabolic) must not appear; only GEM_transporter: is kept
        assert len(df[df['resource'].str.startswith('GEM:')]) == 0

    def test_disables_brenda_by_default(self, _mock_build_fn):
        build_transporters()
        assert _mock_build_fn.call_args.kwargs.get('brenda') is False

    def test_disables_mrclinksdb_by_default(self, _mock_build_fn):
        build_transporters()
        assert _mock_build_fn.call_args.kwargs.get('mrclinksdb') is False

    def test_can_reenable_brenda(self, _mock_build_fn):
        build_transporters(brenda={})
        assert _mock_build_fn.call_args.kwargs.get('brenda') == {}


# ---------------------------------------------------------------------------
# TestBuildReceptors
# ---------------------------------------------------------------------------

class TestBuildReceptors:
    """Tests for build_receptors() filtering and default-disabled resources."""

    def test_returns_dataframe(self, _mock_build_fn):
        assert isinstance(build_receptors(), pd.DataFrame)

    def test_index_reset(self, _mock_build_fn):
        df = build_receptors()
        assert list(df.index) == list(range(len(df)))

    def test_total_row_count(self, _mock_build_fn):
        # MRCLinksDB ligand_receptor + STITCH receptor
        assert len(build_receptors()) == 2

    def test_keeps_ligand_receptor(self, _mock_build_fn):
        df = build_receptors()
        assert 'ligand_receptor' in df['interaction_type'].values

    def test_keeps_stitch_receptor(self, _mock_build_fn):
        df = build_receptors()
        stitch_rows = df[df['resource'].eq('STITCH')]
        assert len(stitch_rows) == 1
        assert stitch_rows.iloc[0]['interaction_type'] == 'receptor'

    def test_excludes_transport_rows(self, _mock_build_fn):
        df = build_receptors()
        assert 'transport' not in df['interaction_type'].values
        assert 'transporter' not in df['interaction_type'].values

    def test_excludes_allosteric_rows(self, _mock_build_fn):
        df = build_receptors()
        assert 'allosteric_regulation' not in df['interaction_type'].values
        assert 'other' not in df['interaction_type'].values

    def test_excludes_gem_rows(self, _mock_build_fn):
        df = build_receptors()
        assert len(df[df['resource'].str.startswith('GEM')]) == 0

    def test_disables_tcdb_by_default(self, _mock_build_fn):
        build_receptors()
        assert _mock_build_fn.call_args.kwargs.get('tcdb') is False

    def test_disables_slc_by_default(self, _mock_build_fn):
        build_receptors()
        assert _mock_build_fn.call_args.kwargs.get('slc') is False

    def test_disables_brenda_by_default(self, _mock_build_fn):
        build_receptors()
        assert _mock_build_fn.call_args.kwargs.get('brenda') is False

    def test_disables_gem_by_default(self, _mock_build_fn):
        build_receptors()
        assert _mock_build_fn.call_args.kwargs.get('gem') is False

    def test_disables_recon3d_by_default(self, _mock_build_fn):
        build_receptors()
        assert _mock_build_fn.call_args.kwargs.get('recon3d') is False

    def test_can_reenable_tcdb(self, _mock_build_fn):
        build_receptors(tcdb={})
        assert _mock_build_fn.call_args.kwargs.get('tcdb') == {}


# ---------------------------------------------------------------------------
# TestBuildAllosteric
# ---------------------------------------------------------------------------

class TestBuildAllosteric:
    """Tests for build_allosteric() filtering and default-disabled resources."""

    def test_returns_dataframe(self, _mock_build_fn):
        assert isinstance(build_allosteric(), pd.DataFrame)

    def test_index_reset(self, _mock_build_fn):
        df = build_allosteric()
        assert list(df.index) == list(range(len(df)))

    def test_total_row_count(self, _mock_build_fn):
        # BRENDA allosteric_regulation + STITCH other
        assert len(build_allosteric()) == 2

    def test_keeps_allosteric_regulation(self, _mock_build_fn):
        df = build_allosteric()
        assert 'allosteric_regulation' in df['interaction_type'].values

    def test_keeps_stitch_other(self, _mock_build_fn):
        df = build_allosteric()
        stitch_rows = df[df['resource'].eq('STITCH')]
        assert len(stitch_rows) == 1
        assert stitch_rows.iloc[0]['interaction_type'] == 'other'

    def test_excludes_transport_rows(self, _mock_build_fn):
        df = build_allosteric()
        assert 'transport' not in df['interaction_type'].values
        assert 'transporter' not in df['interaction_type'].values

    def test_excludes_receptor_rows(self, _mock_build_fn):
        df = build_allosteric()
        assert 'ligand_receptor' not in df['interaction_type'].values
        assert 'receptor' not in df['interaction_type'].values

    def test_excludes_gem_rows(self, _mock_build_fn):
        df = build_allosteric()
        assert len(df[df['resource'].str.startswith('GEM')]) == 0

    def test_disables_tcdb_by_default(self, _mock_build_fn):
        build_allosteric()
        assert _mock_build_fn.call_args.kwargs.get('tcdb') is False

    def test_disables_slc_by_default(self, _mock_build_fn):
        build_allosteric()
        assert _mock_build_fn.call_args.kwargs.get('slc') is False

    def test_disables_mrclinksdb_by_default(self, _mock_build_fn):
        build_allosteric()
        assert _mock_build_fn.call_args.kwargs.get('mrclinksdb') is False

    def test_disables_gem_by_default(self, _mock_build_fn):
        build_allosteric()
        assert _mock_build_fn.call_args.kwargs.get('gem') is False

    def test_disables_recon3d_by_default(self, _mock_build_fn):
        build_allosteric()
        assert _mock_build_fn.call_args.kwargs.get('recon3d') is False

    def test_can_reenable_gem(self, _mock_build_fn):
        build_allosteric(gem={})
        assert _mock_build_fn.call_args.kwargs.get('gem') == {}


# ---------------------------------------------------------------------------
# TestBuildEnzymeMetabolite
# ---------------------------------------------------------------------------

class TestBuildEnzymeMetabolite:
    """Tests for build_enzyme_metabolite() filtering and default-disabled resources.

    Only stoichiometric GEM metabolic reactions (resource starts with 'GEM:')
    are kept.  Transporter GEM edges ('GEM_transporter:') and allosteric
    resources (BRENDA, STITCH) are excluded.
    """

    def test_returns_dataframe(self, _mock_build_fn):
        assert isinstance(build_enzyme_metabolite(), pd.DataFrame)

    def test_index_reset(self, _mock_build_fn):
        df = build_enzyme_metabolite()
        assert list(df.index) == list(range(len(df)))

    def test_total_row_count(self, _mock_build_fn):
        # only GEM:Human-GEM metabolic row
        assert len(build_enzyme_metabolite()) == 1

    def test_keeps_gem_metabolic(self, _mock_build_fn):
        df = build_enzyme_metabolite()
        assert df['resource'].str.startswith('GEM:').all()

    def test_excludes_gem_transporter(self, _mock_build_fn):
        df = build_enzyme_metabolite()
        # 'GEM_transporter:...' does NOT start with 'GEM:' → excluded
        assert len(df[df['resource'].str.startswith('GEM_transporter')]) == 0

    def test_excludes_allosteric_rows(self, _mock_build_fn):
        df = build_enzyme_metabolite()
        assert 'allosteric_regulation' not in df['interaction_type'].values

    def test_excludes_stitch(self, _mock_build_fn):
        df = build_enzyme_metabolite()
        assert 'STITCH' not in df['resource'].values

    def test_excludes_transport_rows(self, _mock_build_fn):
        df = build_enzyme_metabolite()
        assert 'transport' not in df['interaction_type'].values

    def test_excludes_receptor_rows(self, _mock_build_fn):
        df = build_enzyme_metabolite()
        assert 'ligand_receptor' not in df['interaction_type'].values

    def test_disables_brenda_by_default(self, _mock_build_fn):
        build_enzyme_metabolite()
        assert _mock_build_fn.call_args.kwargs.get('brenda') is False

    def test_disables_stitch_by_default(self, _mock_build_fn):
        build_enzyme_metabolite()
        assert _mock_build_fn.call_args.kwargs.get('stitch') is False

    def test_disables_tcdb_by_default(self, _mock_build_fn):
        build_enzyme_metabolite()
        assert _mock_build_fn.call_args.kwargs.get('tcdb') is False

    def test_disables_slc_by_default(self, _mock_build_fn):
        build_enzyme_metabolite()
        assert _mock_build_fn.call_args.kwargs.get('slc') is False

    def test_disables_mrclinksdb_by_default(self, _mock_build_fn):
        build_enzyme_metabolite()
        assert _mock_build_fn.call_args.kwargs.get('mrclinksdb') is False

    def test_disables_recon3d_by_default(self, _mock_build_fn):
        build_enzyme_metabolite()
        assert _mock_build_fn.call_args.kwargs.get('recon3d') is False

    def test_can_reenable_brenda(self, _mock_build_fn):
        build_enzyme_metabolite(brenda={})
        assert _mock_build_fn.call_args.kwargs.get('brenda') == {}
