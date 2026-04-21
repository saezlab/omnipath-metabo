#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos.resources.omnipath module.

Tests use mocks to avoid real HTTP calls to omnipathdb.org.
"""

from unittest.mock import MagicMock, patch

import pytest

from omnipath_metabo.datasets.cosmos._record import Interaction
from omnipath_metabo.datasets.cosmos.resources.omnipath import (
    _query_omnipath,
    _sign_to_mor,
    grn_interactions,
    ppi_interactions,
)


# ---------------------------------------------------------------------------
# Fake OmniPath TSV responses
# ---------------------------------------------------------------------------

_PPI_TSV = (
    'source\ttarget\tsource_genesymbol\ttarget_genesymbol\t'
    'is_directed\tis_stimulation\tis_inhibition\tsources\treferences\n'
    'P00533\tP01133\tEGFR\tEGF\t1\t1\t0\tSIGNOR\t1234567\n'
    'P04049\tP15056\tRAF1\tBRAF\t1\t0\t1\tSIGNOR\t2345678\n'
    'P06241\tP12931\tFYN\tSRC\t1\t0\t0\tSIGNOR\t3456789\n'
)

_GRN_TSV = (
    'source\ttarget\tsource_genesymbol\ttarget_genesymbol\t'
    'is_directed\tis_stimulation\tis_inhibition\tsources\treferences\n'
    'P04637\tP21397\tTP53\tMAOA\t1\t1\t0\tCollecTRI\t9876543\n'
    'Q01094\tP38936\tE2F1\tCDKN1A\t1\t0\t1\tCollecTRI\t8765432\n'
)


def _mock_urlopen(tsv_text):
    """Create a mock context manager for urllib.request.urlopen."""

    resp = MagicMock()
    resp.read.return_value = tsv_text.encode('utf-8')
    resp.__enter__ = MagicMock(return_value=resp)
    resp.__exit__ = MagicMock(return_value=False)
    return resp


# ---------------------------------------------------------------------------
# _sign_to_mor
# ---------------------------------------------------------------------------

class TestSignToMor:
    """Tests for the stimulation/inhibition to MOR converter."""

    def test_stimulation_only(self):
        assert _sign_to_mor({'is_stimulation': '1', 'is_inhibition': '0'}) == 1

    def test_inhibition_only(self):
        assert _sign_to_mor({'is_stimulation': '0', 'is_inhibition': '1'}) == -1

    def test_both_returns_zero(self):
        assert _sign_to_mor({'is_stimulation': '1', 'is_inhibition': '1'}) == 0

    def test_neither_returns_zero(self):
        assert _sign_to_mor({'is_stimulation': '0', 'is_inhibition': '0'}) == 0

    def test_missing_keys_returns_zero(self):
        assert _sign_to_mor({}) == 0

    def test_stimulation_missing_inhibition(self):
        assert _sign_to_mor({'is_stimulation': '1'}) == 1

    def test_inhibition_missing_stimulation(self):
        assert _sign_to_mor({'is_inhibition': '1'}) == -1


# ---------------------------------------------------------------------------
# _query_omnipath
# ---------------------------------------------------------------------------

class TestQueryOmnipath:
    """Tests for OmniPath web API query function."""

    @patch('urllib.request.urlopen')
    def test_returns_list_of_dicts(self, mock_urlopen):
        mock_urlopen.return_value = _mock_urlopen(_PPI_TSV)
        rows = _query_omnipath(datasets='omnipath')
        assert isinstance(rows, list)
        assert len(rows) == 3
        assert isinstance(rows[0], dict)

    @patch('urllib.request.urlopen')
    def test_url_contains_endpoint_and_params(self, mock_urlopen):
        mock_urlopen.return_value = _mock_urlopen(_PPI_TSV)
        _query_omnipath(
            endpoint='interactions',
            datasets='omnipath,ligrecextra',
            organism=9606,
        )

        call_args = mock_urlopen.call_args
        url = call_args[0][0]
        assert 'omnipathdb.org/interactions?' in url
        assert 'datasets=omnipath,ligrecextra' in url
        assert 'organisms=9606' in url
        assert 'genesymbols=yes' in url

    @patch('urllib.request.urlopen')
    def test_tsv_fields_parsed(self, mock_urlopen):
        mock_urlopen.return_value = _mock_urlopen(_PPI_TSV)
        rows = _query_omnipath(datasets='omnipath')
        assert rows[0]['source_genesymbol'] == 'EGFR'
        assert rows[0]['target_genesymbol'] == 'EGF'
        assert rows[0]['is_stimulation'] == '1'

    @patch('urllib.request.urlopen')
    def test_organism_in_url(self, mock_urlopen):
        mock_urlopen.return_value = _mock_urlopen(_PPI_TSV)
        _query_omnipath(organism=10090)

        url = mock_urlopen.call_args[0][0]
        assert 'organisms=10090' in url

    @patch('urllib.request.urlopen')
    def test_resources_filter_in_url(self, mock_urlopen):
        mock_urlopen.return_value = _mock_urlopen(_PPI_TSV)
        _query_omnipath(resources='SIGNOR')

        url = mock_urlopen.call_args[0][0]
        assert 'resources=SIGNOR' in url

    @patch('urllib.request.urlopen')
    def test_no_datasets_omits_param(self, mock_urlopen):
        mock_urlopen.return_value = _mock_urlopen(_PPI_TSV)
        _query_omnipath(datasets=None)

        url = mock_urlopen.call_args[0][0]
        assert 'datasets=' not in url


# ---------------------------------------------------------------------------
# ppi_interactions
# ---------------------------------------------------------------------------

class TestPpiInteractions:
    """Tests for PPI interaction generator."""

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_yields_interaction_records(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': 'EGFR',
                'target_genesymbol': 'EGF',
                'is_stimulation': '1',
                'is_inhibition': '0',
                'sources': 'SIGNOR',
                'references': '1234567',
                'is_directed': '1',
            },
        ]

        records = list(ppi_interactions())
        assert len(records) == 1
        assert isinstance(records[0], Interaction)

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_correct_fields(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': 'EGFR',
                'target_genesymbol': 'EGF',
                'is_stimulation': '1',
                'is_inhibition': '0',
                'sources': 'SIGNOR',
                'references': '1234567',
                'is_directed': '1',
            },
        ]

        rec = list(ppi_interactions())[0]
        assert rec.source == 'EGFR'
        assert rec.target == 'EGF'
        assert rec.source_type == 'protein'
        assert rec.target_type == 'protein'
        assert rec.id_type_a == 'genesymbol'
        assert rec.id_type_b == 'genesymbol'
        assert rec.interaction_type == 'signaling'
        assert rec.mor == 1
        assert rec.locations == ()

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_resource_includes_datasets(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': 'A',
                'target_genesymbol': 'B',
                'is_stimulation': '0',
                'is_inhibition': '0',
                'sources': '',
                'references': '',
                'is_directed': '0',
            },
        ]

        rec = list(ppi_interactions(datasets='omnipath,ligrecextra'))[0]
        assert rec.resource == 'OmniPath:omnipath,ligrecextra'

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_inhibition_mor(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': 'RAF1',
                'target_genesymbol': 'BRAF',
                'is_stimulation': '0',
                'is_inhibition': '1',
                'sources': '',
                'references': '',
                'is_directed': '1',
            },
        ]

        rec = list(ppi_interactions())[0]
        assert rec.mor == -1

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_skips_empty_genesymbols(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': '',
                'target_genesymbol': 'B',
                'is_stimulation': '1',
                'is_inhibition': '0',
                'sources': '',
                'references': '',
                'is_directed': '1',
            },
            {
                'source_genesymbol': 'A',
                'target_genesymbol': '',
                'is_stimulation': '1',
                'is_inhibition': '0',
                'sources': '',
                'references': '',
                'is_directed': '1',
            },
        ]

        records = list(ppi_interactions())
        assert len(records) == 0

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_attrs_contain_sources_and_references(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': 'EGFR',
                'target_genesymbol': 'EGF',
                'is_stimulation': '1',
                'is_inhibition': '0',
                'sources': 'SIGNOR;Reactome',
                'references': '123;456',
                'is_directed': '1',
            },
        ]

        rec = list(ppi_interactions())[0]
        assert rec.attrs['sources'] == 'SIGNOR;Reactome'
        assert rec.attrs['references'] == '123;456'
        assert rec.attrs['directed'] is True

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_organism_forwarded(self, mock_query):
        mock_query.return_value = []
        list(ppi_interactions(organism=10090))
        mock_query.assert_called_once()
        assert mock_query.call_args.kwargs.get('organism') == 10090


# ---------------------------------------------------------------------------
# grn_interactions
# ---------------------------------------------------------------------------

class TestGrnInteractions:
    """Tests for GRN interaction generator."""

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_yields_interaction_records(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': 'TP53',
                'target_genesymbol': 'MAOA',
                'is_stimulation': '1',
                'is_inhibition': '0',
                'sources': 'CollecTRI',
                'references': '9876543',
                'is_directed': '1',
            },
        ]

        records = list(grn_interactions())
        assert len(records) == 1
        assert isinstance(records[0], Interaction)

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_interaction_type_is_gene_regulation(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': 'TP53',
                'target_genesymbol': 'MAOA',
                'is_stimulation': '1',
                'is_inhibition': '0',
                'sources': 'CollecTRI',
                'references': '9876543',
                'is_directed': '1',
            },
        ]

        rec = list(grn_interactions())[0]
        assert rec.interaction_type == 'gene_regulation'

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_default_datasets_is_collectri(self, mock_query):
        mock_query.return_value = []
        list(grn_interactions())
        assert mock_query.call_args.kwargs.get('datasets') == 'collectri'

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_resource_label_includes_datasets(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': 'A',
                'target_genesymbol': 'B',
                'is_stimulation': '0',
                'is_inhibition': '0',
                'sources': '',
                'references': '',
                'is_directed': '1',
            },
        ]

        rec = list(grn_interactions(datasets='collectri,dorothea'))[0]
        assert rec.resource == 'OmniPath:collectri,dorothea'

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_dorothea_levels_default(self, mock_query):
        """When dorothea is in datasets, default levels A,B,C are sent."""
        mock_query.return_value = []
        list(grn_interactions(datasets='collectri,dorothea'))
        assert mock_query.call_args.kwargs.get('dorothea_levels') == 'A,B,C'

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_dorothea_levels_custom(self, mock_query):
        mock_query.return_value = []
        list(grn_interactions(datasets='dorothea', dorothea_levels='A,B'))
        assert mock_query.call_args.kwargs.get('dorothea_levels') == 'A,B'

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_no_dorothea_no_levels(self, mock_query):
        """When dorothea is not in datasets, no levels param is sent."""
        mock_query.return_value = []
        list(grn_interactions(datasets='collectri'))
        assert 'dorothea_levels' not in mock_query.call_args.kwargs

    @patch(
        'omnipath_metabo.datasets.cosmos.resources.omnipath._query_omnipath',
    )
    def test_protein_types(self, mock_query):
        mock_query.return_value = [
            {
                'source_genesymbol': 'TP53',
                'target_genesymbol': 'CDKN1A',
                'is_stimulation': '1',
                'is_inhibition': '0',
                'sources': '',
                'references': '',
                'is_directed': '1',
            },
        ]

        rec = list(grn_interactions())[0]
        assert rec.source_type == 'protein'
        assert rec.target_type == 'protein'
        assert rec.id_type_a == 'genesymbol'
        assert rec.id_type_b == 'genesymbol'
