#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos._orthology module."""

from unittest.mock import patch

import pandas as pd
import pytest

from omnipath_metabo.datasets.cosmos._orthology import (
    _translate_id,
    translate_bundle_by_orthology,
)


# ---------------------------------------------------------------------------
# Fixtures: minimal DataFrames and orthology tables
# ---------------------------------------------------------------------------

def _make_df(rows: list[dict]) -> pd.DataFrame:
    """Build a minimal PKN-like DataFrame."""
    return pd.DataFrame(rows)


def _human_row(source='P00533', target='P01133', resource='SLC'):
    """One protein-protein row from a human-only resource."""
    return {
        'source': source,
        'target': target,
        'source_type': 'protein',
        'target_type': 'protein',
        'resource': resource,
        'mor': 1,
    }


def _metabolite_row(source='CHEBI:15422', target='P00533', resource='SLC'):
    """One metabolite-protein row."""
    return {
        'source': source,
        'target': target,
        'source_type': 'small_molecule',
        'target_type': 'protein',
        'resource': resource,
        'mor': 1,
    }


def _other_resource_row(source='Q99999', target='P88888', resource='TCDB'):
    """Row from a non-human-only resource (should not be translated)."""
    return {
        'source': source,
        'target': target,
        'source_type': 'protein',
        'target_type': 'protein',
        'resource': resource,
        'mor': 1,
    }


# Simple orthology table: human -> mouse
_ORTH_TABLE = {
    'P00533': {'Q01279'},      # EGFR -> mouse EGFR
    'P01133': {'P01132'},      # EGF -> mouse EGF
}

# Orthology table with 1:many mapping
_ORTH_TABLE_MULTI = {
    'P00533': {'Q01279', 'Q99999'},
}


# ---------------------------------------------------------------------------
# _translate_id
# ---------------------------------------------------------------------------

class TestTranslateId:
    """Tests for the single-identifier translator."""

    def test_protein_found(self):
        result = _translate_id('P00533', 'protein', _ORTH_TABLE)
        assert result == 'Q01279'

    def test_protein_not_found(self):
        result = _translate_id('P99999', 'protein', _ORTH_TABLE)
        assert result is None

    def test_metabolite_passes_through(self):
        result = _translate_id('CHEBI:15422', 'small_molecule', _ORTH_TABLE)
        assert result == 'CHEBI:15422'

    def test_frozenset_translated(self):
        result = _translate_id(frozenset({'P00533'}), 'protein', _ORTH_TABLE)
        assert result == frozenset({'Q01279'})

    def test_frozenset_empty_returns_none(self):
        result = _translate_id(frozenset({'P99999'}), 'protein', _ORTH_TABLE)
        assert result is None

    def test_frozenset_multi_input(self):
        result = _translate_id(
            frozenset({'P00533', 'P01133'}),
            'protein',
            _ORTH_TABLE,
        )
        assert result == frozenset({'Q01279', 'P01132'})

    def test_one_to_many_returns_frozenset(self):
        result = _translate_id('P00533', 'protein', _ORTH_TABLE_MULTI)
        assert isinstance(result, frozenset)
        assert result == frozenset({'Q01279', 'Q99999'})


# ---------------------------------------------------------------------------
# translate_bundle_by_orthology
# ---------------------------------------------------------------------------

class TestTranslateBundleByOrthology:
    """Tests for the bundle-level orthology translator."""

    @patch(
        'omnipath_metabo.datasets.cosmos._orthology._get_orthology_table',
    )
    def test_same_organism_returns_unchanged(self, mock_orth):
        df = _make_df([_human_row()])
        result = translate_bundle_by_orthology(
            df,
            source_organism=9606,
            target_organism=9606,
        )
        mock_orth.assert_not_called()
        pd.testing.assert_frame_equal(result, df)

    @patch(
        'omnipath_metabo.datasets.cosmos._orthology._get_orthology_table',
        return_value=_ORTH_TABLE,
    )
    def test_translates_human_only_resources(self, mock_orth):
        df = _make_df([_human_row(source='P00533', target='P01133')])
        result = translate_bundle_by_orthology(df, target_organism=10090)

        assert len(result) == 1
        assert result.iloc[0]['source'] == 'Q01279'
        assert result.iloc[0]['target'] == 'P01132'

    @patch(
        'omnipath_metabo.datasets.cosmos._orthology._get_orthology_table',
        return_value=_ORTH_TABLE,
    )
    def test_drops_rows_without_ortholog(self, mock_orth):
        df = _make_df([
            _human_row(source='P00533', target='P01133'),  # both have orthologs
            _human_row(source='P00533', target='PXXXXX'),  # target has no ortholog
        ])
        result = translate_bundle_by_orthology(df, target_organism=10090)

        assert len(result) == 1
        assert result.iloc[0]['source'] == 'Q01279'

    @patch(
        'omnipath_metabo.datasets.cosmos._orthology._get_orthology_table',
        return_value=_ORTH_TABLE,
    )
    def test_non_target_resources_unchanged(self, mock_orth):
        df = _make_df([
            _human_row(source='P00533', target='P01133'),
            _other_resource_row(),
        ])
        result = translate_bundle_by_orthology(df, target_organism=10090)

        # The TCDB row should be unchanged
        tcdb_rows = result[result['resource'] == 'TCDB']
        assert len(tcdb_rows) == 1
        assert tcdb_rows.iloc[0]['source'] == 'Q99999'

    @patch(
        'omnipath_metabo.datasets.cosmos._orthology._get_orthology_table',
        return_value=_ORTH_TABLE,
    )
    def test_metabolite_ids_not_translated(self, mock_orth):
        df = _make_df([
            _metabolite_row(source='CHEBI:15422', target='P00533'),
        ])
        result = translate_bundle_by_orthology(df, target_organism=10090)

        assert len(result) == 1
        assert result.iloc[0]['source'] == 'CHEBI:15422'  # metabolite unchanged
        assert result.iloc[0]['target'] == 'Q01279'       # protein translated

    @patch(
        'omnipath_metabo.datasets.cosmos._orthology._get_orthology_table',
        return_value={},
    )
    def test_empty_orthology_table_drops_human_rows(self, mock_orth):
        df = _make_df([
            _human_row(),
            _other_resource_row(),
        ])
        result = translate_bundle_by_orthology(df, target_organism=10090)

        # Human-only rows dropped, non-human-only rows kept
        assert len(result) == 1
        assert result.iloc[0]['resource'] == 'TCDB'

    @patch(
        'omnipath_metabo.datasets.cosmos._orthology._get_orthology_table',
        return_value=_ORTH_TABLE,
    )
    def test_no_human_only_rows_returns_unchanged(self, mock_orth):
        df = _make_df([_other_resource_row()])
        result = translate_bundle_by_orthology(df, target_organism=10090)

        assert len(result) == 1
        assert result.iloc[0]['resource'] == 'TCDB'

    @patch(
        'omnipath_metabo.datasets.cosmos._orthology._get_orthology_table',
        return_value=_ORTH_TABLE,
    )
    def test_custom_resources_set(self, mock_orth):
        """Can specify a custom set of resources to translate."""
        df = _make_df([
            _other_resource_row(source='P00533', target='P01133', resource='TCDB'),
        ])
        result = translate_bundle_by_orthology(
            df,
            target_organism=10090,
            resources={'TCDB'},
        )

        assert len(result) == 1
        assert result.iloc[0]['source'] == 'Q01279'

    @patch(
        'omnipath_metabo.datasets.cosmos._orthology._get_orthology_table',
        return_value=_ORTH_TABLE,
    )
    def test_recon3d_rows_translated(self, mock_orth):
        """Recon3D is in the default human-only resources set."""
        df = _make_df([
            _human_row(source='P00533', target='P01133', resource='Recon3D'),
        ])
        result = translate_bundle_by_orthology(df, target_organism=10090)

        assert len(result) == 1
        assert result.iloc[0]['source'] == 'Q01279'
