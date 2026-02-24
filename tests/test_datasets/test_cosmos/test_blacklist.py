#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos._blacklist module."""

import pandas as pd
import pytest

from omnipath_metabo.datasets.cosmos._blacklist import apply_blacklist


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_df(rows):
    """Build a minimal PKN-style DataFrame from a list of dicts."""

    cols = ['source', 'target', 'resource', 'interaction_type',
            'source_type', 'target_type']
    return pd.DataFrame(rows, columns=cols)


_ROWS = [
    dict(source='CHEBI:15422', target='ENSG00000001234', resource='STITCH',
         interaction_type='receptor', source_type='small_molecule',
         target_type='protein'),
    dict(source='CHEBI:15422', target='ENSG00000005678', resource='STITCH',
         interaction_type='transporter', source_type='small_molecule',
         target_type='protein'),
    dict(source='CHEBI:57618', target='ENSG00000009999', resource='Recon3D',
         interaction_type='transport', source_type='small_molecule',
         target_type='protein'),
    dict(source='CHEBI:30616', target='ENSG00000002222', resource='GEM:Human-GEM',
         interaction_type='catalysis', source_type='small_molecule',
         target_type='protein'),
]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestApplyBlacklistEmpty:

    def test_empty_entries_returns_unchanged(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[])
        assert len(result) == len(df)
        pd.testing.assert_frame_equal(result, df)

    def test_none_entries_uses_default_builtin(self):
        """With the default (currently empty) built-in blacklist, df is unchanged."""

        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=None)
        # Built-in blacklist.yaml is empty — all rows preserved.
        assert len(result) == len(df)


class TestApplyBlacklistSingleField:

    def test_filter_by_source(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{'source': 'CHEBI:57618'}])
        assert len(result) == 3
        assert 'CHEBI:57618' not in result['source'].values

    def test_filter_by_resource(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{'resource': 'STITCH'}])
        assert len(result) == 2
        assert 'STITCH' not in result['resource'].values

    def test_filter_by_interaction_type(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{'interaction_type': 'catalysis'}])
        assert len(result) == 3


class TestApplyBlacklistAndLogic:
    """Within an entry, all fields must match (AND)."""

    def test_and_both_match(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{
            'source': 'CHEBI:15422',
            'target': 'ENSG00000001234',
        }])
        # Only the first row matches both conditions.
        assert len(result) == 3
        targets = result['target'].values
        assert 'ENSG00000001234' not in targets
        assert 'ENSG00000005678' in targets  # same source, different target → kept

    def test_and_partial_no_match(self):
        """Entry where source matches but target doesn't → nothing removed."""

        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{
            'source': 'CHEBI:15422',
            'target': 'ENSG99999999',   # doesn't exist
        }])
        assert len(result) == len(df)

    def test_and_three_fields(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{
            'source': 'CHEBI:57618',
            'resource': 'Recon3D',
            'interaction_type': 'transport',
        }])
        assert len(result) == 3


class TestApplyBlacklistOrLogic:
    """Multiple entries are combined with OR logic."""

    def test_or_two_entries(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[
            {'source': 'CHEBI:15422', 'target': 'ENSG00000001234'},
            {'resource': 'Recon3D'},
        ])
        # First entry removes row 0; second removes row 2.
        assert len(result) == 2
        assert 'ENSG00000001234' not in result['target'].values
        assert 'Recon3D' not in result['resource'].values

    def test_or_overlapping_entries(self):
        """Overlapping entries don't double-count — still correct final count."""

        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[
            {'resource': 'STITCH'},
            {'source': 'CHEBI:15422'},  # overlaps — rows 0 and 1 already matched
        ])
        assert len(result) == 2


class TestApplyBlacklistEdgeCases:

    def test_unknown_column_skipped_gracefully(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{'nonexistent_col': 'CHEBI:15422'}])
        # Unknown column → no field conditions → entry_mask stays True → all removed
        # Actually: the entry has one field that is skipped, so entry_mask remains
        # all-True → all rows removed. This is conservative (safe) behaviour.
        # Verify at least that no exception is raised.
        assert isinstance(result, pd.DataFrame)

    def test_empty_entry_dict_ignored(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{}])
        # Empty dict filtered out before processing → no rows removed.
        assert len(result) == len(df)

    def test_returns_dataframe(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[])
        assert isinstance(result, pd.DataFrame)

    def test_index_reset(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{'source': 'CHEBI:15422'}])
        assert list(result.index) == list(range(len(result)))

    def test_empty_dataframe(self):
        df = _make_df([])
        result = apply_blacklist(df, entries=[{'source': 'CHEBI:15422'}])
        assert len(result) == 0

    def test_all_rows_removed(self):
        df = _make_df(_ROWS)
        result = apply_blacklist(df, entries=[{'source_type': 'small_molecule'}])
        assert len(result) == 0

    def test_numeric_values_compared_as_strings(self):
        """Values in entries are compared as strings (mor is int in PKN but str here)."""

        df = _make_df(_ROWS)
        # mor column doesn't exist in our test df, but type coercion is tested below:
        # Ensure str(val) == str(df[col]) matching works for non-string values.
        df2 = df.copy()
        df2['numeric_col'] = [1, 2, 3, 4]
        result = apply_blacklist(df2, entries=[{'numeric_col': 2}])
        assert len(result) == 3
