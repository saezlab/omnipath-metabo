#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos._cache module."""

import shutil
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from omnipath_metabo.datasets.cosmos._cache import (
    _save_bundle,
    list_cached,
    load_cached,
)
from omnipath_metabo.datasets.cosmos._record import Interaction


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_cache(tmp_path):
    """Temporary cache directory."""
    cache_dir = tmp_path / 'cosmos_cache'
    cache_dir.mkdir()
    return cache_dir


def _make_bundle(network):
    """Build a minimal mock bundle with a .network attribute."""
    bundle = MagicMock()
    bundle.network = network
    return bundle


def _sample_interactions():
    """Return a list of Interaction records for testing."""
    return [
        Interaction(
            source='CHEBI:15422',
            target='P00533',
            source_type='small_molecule',
            target_type='protein',
            id_type_a='chebi',
            id_type_b='uniprot',
            interaction_type='transport',
            resource='SLC',
            mor=1,
            locations=('e', 'c'),
            attrs={'key': 'val'},
        ),
        Interaction(
            source=frozenset({'CHEBI:1', 'CHEBI:2'}),
            target='P00534',
            source_type='small_molecule',
            target_type='protein',
            id_type_a='chebi',
            id_type_b='uniprot',
            interaction_type='transport',
            resource='SLC',
            mor=-1,
            locations=('e',),
            attrs={},
        ),
        Interaction(
            source='P00535',
            target=frozenset({'P00536', 'P00537'}),
            source_type='protein',
            target_type='protein',
            id_type_a='uniprot',
            id_type_b='uniprot',
            interaction_type='signaling',
            resource='OmniPath:omnipath',
            mor=1,
            locations=(),
            attrs={},
        ),
    ]


# ---------------------------------------------------------------------------
# _save_bundle
# ---------------------------------------------------------------------------

class TestSaveBundle:
    """Tests for _save_bundle() Parquet serialisation."""

    def test_creates_parquet_file(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions())
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        path = tmp_cache / 'transporters_9606.parquet'
        assert path.exists()

    def test_frozenset_source_converted_to_string(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions())
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        df = pd.read_parquet(tmp_cache / 'transporters_9606.parquet')
        # Row with frozenset source should be semicolon-joined sorted string
        row = df[df['mor'] == -1].iloc[0]
        assert isinstance(row['source'], str)
        assert row['source'] == 'CHEBI:1;CHEBI:2'

    def test_frozenset_target_converted_to_string(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions())
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        df = pd.read_parquet(tmp_cache / 'transporters_9606.parquet')
        # Row with frozenset target
        row = df[df['resource'] == 'OmniPath:omnipath'].iloc[0]
        assert isinstance(row['target'], str)
        assert row['target'] == 'P00536;P00537'

    def test_plain_string_preserved(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions())
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        df = pd.read_parquet(tmp_cache / 'transporters_9606.parquet')
        row = df[df['source'] == 'CHEBI:15422'].iloc[0]
        assert row['source'] == 'CHEBI:15422'
        assert row['target'] == 'P00533'

    def test_file_naming_convention(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions()[:1])
        _save_bundle(bundle, 'ppi', 10090, tmp_cache)

        assert (tmp_cache / 'ppi_10090.parquet').exists()

    def test_row_count_matches(self, tmp_cache):
        interactions = _sample_interactions()
        bundle = _make_bundle(interactions)
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        df = pd.read_parquet(tmp_cache / 'transporters_9606.parquet')
        assert len(df) == len(interactions)


# ---------------------------------------------------------------------------
# load_cached
# ---------------------------------------------------------------------------

class TestLoadCached:
    """Tests for load_cached() Parquet reading."""

    def test_returns_dataframe(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions())
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        result = load_cached('transporters', 9606, tmp_cache)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 3

    def test_missing_file_returns_none(self, tmp_cache):
        result = load_cached('nonexistent', 9606, tmp_cache)
        assert result is None

    def test_columns_present(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions())
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        result = load_cached('transporters', 9606, tmp_cache)
        assert 'source' in result.columns
        assert 'target' in result.columns
        assert 'resource' in result.columns
        assert 'mor' in result.columns

    def test_round_trip_preserves_values(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions())
        _save_bundle(bundle, 'receptors', 9606, tmp_cache)

        df = load_cached('receptors', 9606, tmp_cache)
        first_row = df.iloc[0]
        assert first_row['source'] == 'CHEBI:15422'
        assert first_row['target'] == 'P00533'
        assert first_row['mor'] == 1


# ---------------------------------------------------------------------------
# list_cached
# ---------------------------------------------------------------------------

class TestListCached:
    """Tests for list_cached() cache enumeration."""

    def test_empty_cache(self, tmp_cache):
        result = list_cached(tmp_cache)
        assert result == []

    def test_nonexistent_dir(self, tmp_path):
        result = list_cached(tmp_path / 'does_not_exist')
        assert result == []

    def test_single_file(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions()[:1])
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        result = list_cached(tmp_cache)
        assert len(result) == 1
        assert result[0]['category'] == 'transporters'
        assert result[0]['organism'] == 9606

    def test_multiple_files(self, tmp_cache):
        for cat, org in [('transporters', 9606), ('receptors', 9606), ('ppi', 10090)]:
            bundle = _make_bundle(_sample_interactions()[:1])
            _save_bundle(bundle, cat, org, tmp_cache)

        result = list_cached(tmp_cache)
        assert len(result) == 3

        categories = {r['category'] for r in result}
        organisms = {r['organism'] for r in result}
        assert categories == {'transporters', 'receptors', 'ppi'}
        assert organisms == {9606, 10090}

    def test_entry_has_expected_keys(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions()[:1])
        _save_bundle(bundle, 'allosteric', 9606, tmp_cache)

        result = list_cached(tmp_cache)
        entry = result[0]
        assert 'category' in entry
        assert 'organism' in entry
        assert 'size_mb' in entry
        assert 'modified' in entry

    def test_size_is_positive(self, tmp_cache):
        bundle = _make_bundle(_sample_interactions())
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        result = list_cached(tmp_cache)
        assert result[0]['size_mb'] > 0

    def test_ignores_non_parquet_files(self, tmp_cache):
        (tmp_cache / 'readme.txt').write_text('not a parquet file')
        bundle = _make_bundle(_sample_interactions()[:1])
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        result = list_cached(tmp_cache)
        assert len(result) == 1

    def test_ignores_malformed_filenames(self, tmp_cache):
        """Parquet files without the category_organism pattern are skipped."""
        (tmp_cache / 'badname.parquet').write_bytes(b'fake')
        bundle = _make_bundle(_sample_interactions()[:1])
        _save_bundle(bundle, 'transporters', 9606, tmp_cache)

        result = list_cached(tmp_cache)
        # badname.parquet has no underscore-separated organism -> skipped
        # (or if rsplit yields non-integer, also skipped)
        valid = [r for r in result if r['category'] == 'transporters']
        assert len(valid) == 1
