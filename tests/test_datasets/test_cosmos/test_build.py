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
from omnipath_metabo.datasets.cosmos._bundle import CosmosBundle
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
        # Disable all real resources so _auto_select_gem doesn't inject 'gem'.
        'gem': False,
        'recon3d': False,
        'tcdb': False,
        'slc': False,
        'brenda': False,
        'mrclinksdb': False,
        'mrclinksdb_transporter': False,
        'stitch': False,
        'ppi': False,
        'grn': False,
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
        'translate_ids': False,  # fake IDs can't be translated; skip for unit tests
        'resources': {
            'capture': {},
            'gem': False,
        },
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

    def test_returns_bundle(self, _mock_build):
        assert isinstance(build(), CosmosBundle)

    def test_network_is_list_of_interactions(self, _mock_build):
        bundle = build()
        assert isinstance(bundle.network, list)
        for row in bundle.network:
            assert isinstance(row, Interaction)

    def test_default_runs_all(self, _mock_build):
        assert len(build().network) == 4  # 2 records x 2 fake resources

    def test_single_resource(self, _mock_build):
        assert len(build(fake_b=False).network) == 2

    def test_disable_resource(self, _mock_build):
        assert len(build(fake_a=False, fake_b=False).network) == 0

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

        assert len(build(resources={}).network) == 4  # defaults still active

    def test_bundle_provenance_empty_without_translation(self, _mock_build):
        """With translate_ids=False, provenance lists are empty."""
        bundle = build()
        assert bundle.metabolites == []
        assert bundle.proteins == []
        assert bundle.reactions == []


class TestProcessorsRegistry:
    """Tests for the PROCESSORS constant."""

    def test_all_resources_registered(self):
        expected = {
            'stitch', 'tcdb', 'slc', 'brenda',
            'mrclinksdb', 'mrclinksdb_transporter',
            'gem', 'recon3d',
            'ppi', 'grn',
        }
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
        bundle = build(
            stitch={'score_threshold': 900},
            tcdb=False, slc=False, brenda=False, mrclinksdb=False,
            gem=False, recon3d=False,
        )

        assert isinstance(bundle, CosmosBundle)
        assert len(bundle.network) > 0
        df = pd.DataFrame(bundle.network)
        assert (df['resource'] == 'STITCH').all()
        assert set(df['source_type']) == {'small_molecule'}
        assert set(df['target_type']) == {'protein'}

    def test_build_slc_only(self):
        bundle = build(
            slc={},
            stitch=False, tcdb=False, brenda=False, mrclinksdb=False,
            gem=False, recon3d=False,
        )

        assert isinstance(bundle, CosmosBundle)
        assert len(bundle.network) > 0
        df = pd.DataFrame(bundle.network)
        assert (df['resource'] == 'SLC').all()

    def test_build_brenda_only(self):
        bundle = build(
            brenda={},
            stitch=False, tcdb=False, slc=False, mrclinksdb=False,
            gem=False, recon3d=False,
        )

        assert isinstance(bundle, CosmosBundle)
        assert len(bundle.network) > 0
        df = pd.DataFrame(bundle.network)
        assert (df['resource'] == 'BRENDA').all()


# ---------------------------------------------------------------------------
# Shared fake rows for subset filter tests
# ---------------------------------------------------------------------------

def _row(interaction_type: str, resource: str, mor: int = 1) -> Interaction:
    return Interaction(
        source='CHEBI:1',
        target='P00001',
        source_type='small_molecule',
        target_type='protein',
        id_type_a='chebi',
        id_type_b='uniprot',
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
# MRCLinksDB uses 'ligand_receptor'; STITCH uses 'receptor' (set by
# _classify_protein — the protein's role, not the interaction type).
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
def _all_categories_bundle():
    return CosmosBundle(network=list(_ALL_ROWS))


@pytest.fixture
def _mock_build_fn(_all_categories_bundle):
    """Patch build() to return a controlled CosmosBundle for subset filter tests."""

    def _fake_build(*a, row_filter=None, **kw):
        rows = list(_ALL_ROWS)
        if row_filter is not None:
            rows = [r for r in rows if row_filter(r)]
        return CosmosBundle(network=rows)

    with patch(
        'omnipath_metabo.datasets.cosmos._build.build',
        side_effect=_fake_build,
    ) as mock:
        yield mock


def _net(bundle: CosmosBundle) -> pd.DataFrame:
    """Convert bundle.network to a DataFrame for column-based assertions."""
    return pd.DataFrame(bundle.network)


# ---------------------------------------------------------------------------
# TestBuildTransporters
# ---------------------------------------------------------------------------

class TestBuildTransporters:
    """Tests for build_transporters() filtering and default-disabled resources."""

    def test_returns_bundle(self, _mock_build_fn):
        assert isinstance(build_transporters(), CosmosBundle)

    def test_total_row_count(self, _mock_build_fn):
        # TCDB, SLC, GEM_transporter (×2: by type + by resource), STITCH transporter
        assert len(build_transporters().network) == 5

    def test_keeps_transport_type(self, _mock_build_fn):
        df = _net(build_transporters())
        assert 'transport' in df['interaction_type'].values

    def test_keeps_gem_transporter_by_resource_prefix(self, _mock_build_fn):
        df = _net(build_transporters())
        # both GEM_transporter rows: one matched via type, one only via resource prefix
        gem_rows = df[df['resource'].str.startswith('GEM_transporter')]
        assert len(gem_rows) == 2

    def test_keeps_stitch_transporter(self, _mock_build_fn):
        df = _net(build_transporters())
        stitch_rows = df[df['resource'].eq('STITCH')]
        assert len(stitch_rows) == 1
        assert stitch_rows.iloc[0]['interaction_type'] == 'transporter'

    def test_excludes_receptor_rows(self, _mock_build_fn):
        df = _net(build_transporters())
        assert 'ligand_receptor' not in df['interaction_type'].values
        assert 'receptor' not in df['interaction_type'].values

    def test_excludes_allosteric_rows(self, _mock_build_fn):
        df = _net(build_transporters())
        assert 'allosteric_regulation' not in df['interaction_type'].values

    def test_excludes_gem_metabolic(self, _mock_build_fn):
        df = _net(build_transporters())
        # GEM: prefix (metabolic) must not appear; only GEM_transporter: is kept
        assert len(df[df['resource'].str.startswith('GEM:')]) == 0

    def test_disables_brenda_by_default(self, _mock_build_fn):
        build_transporters()
        assert _mock_build_fn.call_args.kwargs.get('brenda') is False

    def test_enables_mrclinksdb_by_default(self, _mock_build_fn):
        # MRCLinksDB is intentionally included: transport-classified records
        # (interaction_type='transport') contribute to the transporter subset.
        build_transporters()
        assert _mock_build_fn.call_args.kwargs.get('mrclinksdb') is not False

    def test_can_reenable_brenda(self, _mock_build_fn):
        build_transporters(brenda={})
        assert _mock_build_fn.call_args.kwargs.get('brenda') == {}


# ---------------------------------------------------------------------------
# TestBuildReceptors
# ---------------------------------------------------------------------------

class TestBuildReceptors:
    """Tests for build_receptors() filtering and default-disabled resources."""

    def test_returns_bundle(self, _mock_build_fn):
        assert isinstance(build_receptors(), CosmosBundle)

    def test_total_row_count(self, _mock_build_fn):
        # MRCLinksDB ligand_receptor + STITCH receptor
        assert len(build_receptors().network) == 2

    def test_keeps_ligand_receptor(self, _mock_build_fn):
        df = _net(build_receptors())
        assert 'ligand_receptor' in df['interaction_type'].values

    def test_keeps_stitch_receptor(self, _mock_build_fn):
        df = _net(build_receptors())
        stitch_rows = df[df['resource'].eq('STITCH')]
        assert len(stitch_rows) == 1
        assert stitch_rows.iloc[0]['interaction_type'] == 'receptor'

    def test_excludes_transport_rows(self, _mock_build_fn):
        df = _net(build_receptors())
        assert 'transport' not in df['interaction_type'].values
        assert 'transporter' not in df['interaction_type'].values

    def test_excludes_allosteric_rows(self, _mock_build_fn):
        df = _net(build_receptors())
        assert 'allosteric_regulation' not in df['interaction_type'].values
        assert 'other' not in df['interaction_type'].values

    def test_excludes_gem_rows(self, _mock_build_fn):
        df = _net(build_receptors())
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
# Helpers for cell_surface_only tests
# ---------------------------------------------------------------------------

def _receptor_row(locations: tuple) -> Interaction:
    """ligand_receptor row with given locations tuple."""
    return Interaction(
        source='CHEBI:1',
        target='P00001',
        source_type='small_molecule',
        target_type='protein',
        id_type_a='chebi',
        id_type_b='uniprot',
        interaction_type='ligand_receptor',
        resource='MRCLinksDB',
        mor=1,
        locations=locations,
    )


@pytest.fixture
def _mock_build_receptors_fn():
    """Patch build() for cell_surface_only receptor tests with location-bearing rows."""

    rows = [
        _receptor_row(('c', 'e')),   # surface receptor → kept with cell_surface_only
        _receptor_row(('c', 'm')),   # mitochondrial receptor → dropped
        _receptor_row(()),           # no location → dropped
        _receptor_row(('e',)),       # surface only → kept
        Interaction(                  # transport row with 'e' → not a receptor, dropped
            source='CHEBI:2',
            target='P00002',
            source_type='small_molecule',
            target_type='protein',
            id_type_a='chebi',
            id_type_b='uniprot',
            interaction_type='transport',
            resource='MRCLinksDB',
            mor=1,
            locations=('e',),
        ),
    ]

    def _fake_build(*a, row_filter=None, **kw):
        result = list(rows)
        if row_filter is not None:
            result = [r for r in result if row_filter(r)]
        return CosmosBundle(network=result)

    with patch(
        'omnipath_metabo.datasets.cosmos._build.build',
        side_effect=_fake_build,
    ):
        yield


class TestBuildReceptorsCellSurfaceOnly:
    """Tests for build_receptors(cell_surface_only=True/False)."""

    def test_default_keeps_all_ligand_receptor(self, _mock_build_receptors_fn):
        # Without cell_surface_only, all ligand_receptor rows survive regardless of location.
        bundle = build_receptors()
        df = _net(bundle)
        assert len(df[df['interaction_type'] == 'ligand_receptor']) == 4

    def test_cell_surface_only_keeps_e_location(self, _mock_build_receptors_fn):
        bundle = build_receptors(cell_surface_only=True)
        df = _net(bundle)
        surviving = [row for row in bundle.network if 'e' in row.locations]
        assert len(surviving) == 2

    def test_cell_surface_only_drops_no_e_location(self, _mock_build_receptors_fn):
        bundle = build_receptors(cell_surface_only=True)
        assert all('e' in row.locations for row in bundle.network)

    def test_cell_surface_only_drops_empty_location(self, _mock_build_receptors_fn):
        bundle = build_receptors(cell_surface_only=True)
        assert all(row.locations != () for row in bundle.network)

    def test_cell_surface_only_still_requires_ligand_receptor(self, _mock_build_receptors_fn):
        # transport row with 'e' must NOT survive — interaction_type check comes first
        bundle = build_receptors(cell_surface_only=True)
        df = _net(bundle)
        assert 'transport' not in df['interaction_type'].values

    def test_cell_surface_only_false_equivalent_to_default(self, _mock_build_receptors_fn):
        assert (
            len(build_receptors(cell_surface_only=False).network)
            == len(build_receptors().network)
        )


# ---------------------------------------------------------------------------
# TestEnrichStitchLocations
# ---------------------------------------------------------------------------

class TestEnrichStitchLocations:
    """Tests for _enrich_stitch_locations frozenset handling."""

    def _make_df(self, target_val) -> pd.DataFrame:
        """Build a minimal one-row DataFrame with STITCH resource."""
        return pd.DataFrame([{
            'source': 'CHEBI:1',
            'target': target_val,
            'source_type': 'small_molecule',
            'target_type': 'protein',
            'id_type_a': 'chebi',
            'id_type_b': 'uniprot',
            'interaction_type': 'ligand_receptor',
            'resource': 'STITCH',
            'mor': 1,
            'locations': (),
            'attrs': {},
        }])

    # Patch path: functions are imported locally inside _enrich_stitch_locations
    # via `from .location import ...`, so we patch them at their source module.
    _LOC_MOD = 'omnipath_metabo.datasets.cosmos.location'

    def test_string_target_resolved(self):
        """Plain string UniProt AC gets locations resolved."""
        from omnipath_metabo.datasets.cosmos._build import _enrich_stitch_locations

        with (
            patch(f'{self._LOC_MOD}.uniprot_locations', return_value={}),
            patch(f'{self._LOC_MOD}.tcdb_locations', return_value={'Cell membrane': 'e'}),
            patch(
                f'{self._LOC_MOD}.resolve_protein_locations',
                side_effect=lambda uid, *_: {'e'} if uid == 'P00001' else None,
            ),
        ):
            df = self._make_df('P00001')
            result = _enrich_stitch_locations(df, organism=9606)
            assert result.iloc[0]['locations'] == ('e',)

    def test_frozenset_target_resolved(self):
        """Frozenset target (post-translate_pkn) gets locations resolved."""
        from omnipath_metabo.datasets.cosmos._build import _enrich_stitch_locations

        with (
            patch(f'{self._LOC_MOD}.uniprot_locations', return_value={}),
            patch(f'{self._LOC_MOD}.tcdb_locations', return_value={'Cell membrane': 'e'}),
            patch(
                f'{self._LOC_MOD}.resolve_protein_locations',
                side_effect=lambda uid, *_: {'e'} if uid == 'P00001' else None,
            ),
        ):
            df = self._make_df(frozenset({'P00001'}))
            result = _enrich_stitch_locations(df, organism=9606)
            assert result.iloc[0]['locations'] == ('e',)

    def test_frozenset_multi_ac_union(self):
        """Multi-element frozenset: locations are union of all ACs."""
        from omnipath_metabo.datasets.cosmos._build import _enrich_stitch_locations

        def _fake_resolve(uid, *_):
            return {'e', 'c'} if uid == 'P00001' else {'n'} if uid == 'P00002' else None

        with (
            patch(f'{self._LOC_MOD}.uniprot_locations', return_value={}),
            patch(f'{self._LOC_MOD}.tcdb_locations', return_value={}),
            patch(
                f'{self._LOC_MOD}.resolve_protein_locations',
                side_effect=_fake_resolve,
            ),
        ):
            df = self._make_df(frozenset({'P00001', 'P00002'}))
            result = _enrich_stitch_locations(df, organism=9606)
            assert set(result.iloc[0]['locations']) == {'c', 'e', 'n'}

    def test_frozenset_no_location_returns_empty_tuple(self):
        """Frozenset where no AC has location data returns ()."""
        from omnipath_metabo.datasets.cosmos._build import _enrich_stitch_locations

        with (
            patch(f'{self._LOC_MOD}.uniprot_locations', return_value={}),
            patch(f'{self._LOC_MOD}.tcdb_locations', return_value={}),
            patch(f'{self._LOC_MOD}.resolve_protein_locations', return_value=None),
        ):
            df = self._make_df(frozenset({'P99999'}))
            result = _enrich_stitch_locations(df, organism=9606)
            assert result.iloc[0]['locations'] == ()


# ---------------------------------------------------------------------------
# TestBuildAllosteric
# ---------------------------------------------------------------------------

class TestBuildAllosteric:
    """Tests for build_allosteric() filtering and default-disabled resources."""

    def test_returns_bundle(self, _mock_build_fn):
        assert isinstance(build_allosteric(), CosmosBundle)

    def test_total_row_count(self, _mock_build_fn):
        # BRENDA allosteric_regulation + STITCH other
        assert len(build_allosteric().network) == 2

    def test_keeps_allosteric_regulation(self, _mock_build_fn):
        df = _net(build_allosteric())
        assert 'allosteric_regulation' in df['interaction_type'].values

    def test_keeps_stitch_other(self, _mock_build_fn):
        df = _net(build_allosteric())
        stitch_rows = df[df['resource'].eq('STITCH')]
        assert len(stitch_rows) == 1
        assert stitch_rows.iloc[0]['interaction_type'] == 'other'

    def test_excludes_transport_rows(self, _mock_build_fn):
        df = _net(build_allosteric())
        assert 'transport' not in df['interaction_type'].values
        assert 'transporter' not in df['interaction_type'].values

    def test_excludes_receptor_rows(self, _mock_build_fn):
        df = _net(build_allosteric())
        assert 'ligand_receptor' not in df['interaction_type'].values
        assert 'receptor' not in df['interaction_type'].values

    def test_excludes_gem_rows(self, _mock_build_fn):
        df = _net(build_allosteric())
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

    def test_returns_bundle(self, _mock_build_fn):
        assert isinstance(build_enzyme_metabolite(), CosmosBundle)

    def test_total_row_count(self, _mock_build_fn):
        # only GEM:Human-GEM metabolic row
        assert len(build_enzyme_metabolite().network) == 1

    def test_keeps_gem_metabolic(self, _mock_build_fn):
        df = _net(build_enzyme_metabolite())
        assert df['resource'].str.startswith('GEM:').all()

    def test_excludes_gem_transporter(self, _mock_build_fn):
        df = _net(build_enzyme_metabolite())
        # 'GEM_transporter:...' does NOT start with 'GEM:' → excluded
        assert len(df[df['resource'].str.startswith('GEM_transporter')]) == 0

    def test_excludes_allosteric_rows(self, _mock_build_fn):
        df = _net(build_enzyme_metabolite())
        assert 'allosteric_regulation' not in df['interaction_type'].values

    def test_excludes_stitch(self, _mock_build_fn):
        df = _net(build_enzyme_metabolite())
        assert 'STITCH' not in df['resource'].values

    def test_excludes_transport_rows(self, _mock_build_fn):
        df = _net(build_enzyme_metabolite())
        assert 'transport' not in df['interaction_type'].values

    def test_excludes_receptor_rows(self, _mock_build_fn):
        df = _net(build_enzyme_metabolite())
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
