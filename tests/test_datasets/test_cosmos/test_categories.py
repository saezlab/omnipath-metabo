#!/usr/bin/env python

"""Tests for category-specific build functions.

All tests use mocked processors so no external network calls are made.
"""

import copy
from unittest.mock import patch

import pandas as pd
import pytest

from omnipath_metabo.datasets.cosmos._record import Interaction


# ---------------------------------------------------------------------------
# Synthetic interactions for each category
# ---------------------------------------------------------------------------

def _make_interaction(
    interaction_type,
    resource,
    source_type='small_molecule',
    target_type='protein',
):
    return Interaction(
        source='MET001',
        target='PROT001',
        source_type=source_type,
        target_type=target_type,
        id_type_a='pubchem',
        id_type_b='uniprot',
        interaction_type=interaction_type,
        resource=resource,
        mor=0,
    )


# One record per category, to verify that each wrapper keeps and excludes correctly
_TRANSPORT_TCDB = _make_interaction('transport', 'TCDB')
_TRANSPORT_SLC = _make_interaction('transport', 'SLC')
_TRANSPORT_GEM_T = _make_interaction('catalysis', 'GEM_transporter:Human-GEM')
_TRANSPORT_RECON = _make_interaction('transport', 'Recon3D')
_TRANSPORT_STITCH = _make_interaction('transporter', 'STITCH')

_RECEPTOR_MRC = _make_interaction('ligand_receptor', 'MRCLinksDB')
_RECEPTOR_STITCH = _make_interaction('receptor', 'STITCH')

_ENZYME_BRENDA = _make_interaction('allosteric_regulation', 'BRENDA')
_ENZYME_GEM = _make_interaction('catalysis', 'GEM:Human-GEM')
_ENZYME_STITCH = _make_interaction('other', 'STITCH')

_ALL_INTERACTIONS = [
    _TRANSPORT_TCDB,
    _TRANSPORT_SLC,
    _TRANSPORT_GEM_T,
    _TRANSPORT_RECON,
    _TRANSPORT_STITCH,
    _RECEPTOR_MRC,
    _RECEPTOR_STITCH,
    _ENZYME_BRENDA,
    _ENZYME_GEM,
    _ENZYME_STITCH,
]


FAKE_DEFAULT_CONFIG = {
    'organism': 9606,
    'translate_ids': False,
    'apply_blacklist': False,
    'resources': {
        'fake': {},
    },
}


def _make_mock_build(interactions):
    """Return a context manager that makes build() return the given interactions."""

    def _fake_resource(**kwargs):
        yield from interactions

    fake_processors = {'fake': _fake_resource}
    cfg = copy.deepcopy(FAKE_DEFAULT_CONFIG)

    return (
        patch.dict(
            'omnipath_metabo.datasets.cosmos._build.PROCESSORS',
            fake_processors,
            clear=True,
        ),
        patch(
            'omnipath_metabo.datasets.cosmos._config.default_config',
            side_effect=lambda: copy.deepcopy(cfg),
        ),
    )


# ---------------------------------------------------------------------------
# Tests: build_transporters
# ---------------------------------------------------------------------------

class TestBuildTransporters:

    def _run(self, interactions=None):
        from omnipath_metabo.datasets.cosmos import build_transporters

        ints = interactions if interactions is not None else _ALL_INTERACTIONS
        p0, p1 = _make_mock_build(ints)
        with p0, p1:
            return build_transporters()

    def test_returns_dataframe(self):
        df = self._run()
        assert isinstance(df, pd.DataFrame)

    def test_transport_interaction_type_kept(self):
        df = self._run([_TRANSPORT_TCDB, _TRANSPORT_SLC])
        assert len(df) == 2
        assert (df['interaction_type'] == 'transport').all()

    def test_gem_transporter_resource_kept(self):
        df = self._run([_TRANSPORT_GEM_T])
        assert len(df) == 1
        assert df.iloc[0]['resource'] == 'GEM_transporter:Human-GEM'

    def test_stitch_transporter_kept(self):
        df = self._run([_TRANSPORT_STITCH])
        assert len(df) == 1
        assert df.iloc[0]['interaction_type'] == 'transporter'

    def test_no_receptor_rows(self):
        df = self._run(_ALL_INTERACTIONS)
        assert not df['interaction_type'].isin(['ligand_receptor', 'receptor']).any()

    def test_no_enzyme_rows(self):
        df = self._run(_ALL_INTERACTIONS)
        # allosteric_regulation and plain GEM: should be absent
        assert not df['interaction_type'].eq('allosteric_regulation').any()
        # GEM: (non-transporter) resource should not be present
        gem_metabolic = df['resource'].eq('GEM:Human-GEM')
        assert not gem_metabolic.any()

    def test_index_reset(self):
        df = self._run(_ALL_INTERACTIONS)
        assert list(df.index) == list(range(len(df)))


# ---------------------------------------------------------------------------
# Tests: build_receptors
# ---------------------------------------------------------------------------

class TestBuildReceptors:

    def _run(self, interactions=None):
        from omnipath_metabo.datasets.cosmos import build_receptors

        ints = interactions if interactions is not None else _ALL_INTERACTIONS
        p0, p1 = _make_mock_build(ints)
        with p0, p1:
            return build_receptors()

    def test_returns_dataframe(self):
        df = self._run()
        assert isinstance(df, pd.DataFrame)

    def test_ligand_receptor_kept(self):
        df = self._run([_RECEPTOR_MRC])
        assert len(df) == 1
        assert df.iloc[0]['interaction_type'] == 'ligand_receptor'

    def test_stitch_receptor_kept(self):
        df = self._run([_RECEPTOR_STITCH])
        assert len(df) == 1
        assert df.iloc[0]['interaction_type'] == 'receptor'

    def test_no_transport_rows(self):
        df = self._run(_ALL_INTERACTIONS)
        transport_mask = (
            df['interaction_type'].eq('transport') |
            df['resource'].str.startswith('GEM_transporter')
        )
        assert not transport_mask.any()

    def test_no_enzyme_rows(self):
        df = self._run(_ALL_INTERACTIONS)
        assert not df['interaction_type'].eq('allosteric_regulation').any()

    def test_index_reset(self):
        df = self._run(_ALL_INTERACTIONS)
        assert list(df.index) == list(range(len(df)))


# ---------------------------------------------------------------------------
# Tests: build_enzyme_metabolite
# ---------------------------------------------------------------------------

class TestBuildEnzymeMetabolite:

    def _run(self, interactions=None):
        from omnipath_metabo.datasets.cosmos import build_enzyme_metabolite

        ints = interactions if interactions is not None else _ALL_INTERACTIONS
        p0, p1 = _make_mock_build(ints)
        with p0, p1:
            return build_enzyme_metabolite()

    def test_returns_dataframe(self):
        df = self._run()
        assert isinstance(df, pd.DataFrame)

    def test_allosteric_kept(self):
        df = self._run([_ENZYME_BRENDA])
        assert len(df) == 1
        assert df.iloc[0]['interaction_type'] == 'allosteric_regulation'

    def test_gem_metabolic_kept(self):
        df = self._run([_ENZYME_GEM])
        assert len(df) == 1
        assert df.iloc[0]['resource'] == 'GEM:Human-GEM'

    def test_gem_transporter_excluded(self):
        # GEM_transporter does NOT start with 'GEM:' so must not appear
        df = self._run([_TRANSPORT_GEM_T, _ENZYME_GEM])
        resources = set(df['resource'])
        assert 'GEM_transporter:Human-GEM' not in resources
        assert 'GEM:Human-GEM' in resources

    def test_stitch_other_kept(self):
        df = self._run([_ENZYME_STITCH])
        assert len(df) == 1
        assert df.iloc[0]['interaction_type'] == 'other'

    def test_no_transport_rows(self):
        df = self._run(_ALL_INTERACTIONS)
        transport_mask = (
            df['interaction_type'].eq('transport') |
            df['resource'].str.startswith('GEM_transporter')
        )
        assert not transport_mask.any()

    def test_no_receptor_rows(self):
        df = self._run(_ALL_INTERACTIONS)
        receptor_mask = (
            df['interaction_type'].eq('ligand_receptor') |
            (df['resource'].eq('STITCH') & df['interaction_type'].eq('receptor'))
        )
        assert not receptor_mask.any()

    def test_index_reset(self):
        df = self._run(_ALL_INTERACTIONS)
        assert list(df.index) == list(range(len(df)))
