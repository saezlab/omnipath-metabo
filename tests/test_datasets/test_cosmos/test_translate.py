#!/usr/bin/env python

"""Unit tests for COSMOS PKN ID translation.

Fast tests cover the pass-through cases (chebi, ensembl, reaction_id)
that require no external data.  Tests that exercise pypath BioMart
mappings are marked slow.
"""

import pandas as pd
import pytest

from unittest.mock import patch

from omnipath_metabo.datasets.cosmos._record import Interaction
from omnipath_metabo.datasets.cosmos._translate import (
    _normalise_hmdb,
    _to_chebi,
    _to_ensg,
    translate_pkn,
)


def _make_df(rows):
    """Build a PKN DataFrame from a list of Interaction namedtuples."""

    return pd.DataFrame(rows, columns=list(Interaction._fields))


# ---------------------------------------------------------------------------
# _to_chebi dispatch
# ---------------------------------------------------------------------------

class TestToChebi:

    def test_chebi_passthrough(self):
        assert _to_chebi('CHEBI:12345', 'chebi') == 'CHEBI:12345'

    def test_unknown_id_type_returns_none(self):
        result = _to_chebi('X', 'unknown_type')
        assert result is None


# ---------------------------------------------------------------------------
# _to_ensg dispatch
# ---------------------------------------------------------------------------

class TestToEnsg:

    def test_ensembl_passthrough(self):
        assert _to_ensg('ENSG00000001234', 'ensembl', 9606) == 'ENSG00000001234'

    def test_reaction_id_passthrough(self):
        # Orphan pseudo-enzyme IDs must pass through unchanged.
        assert _to_ensg('MAR04831', 'reaction_id', 9606) == 'MAR04831'

    def test_reaction_id_with_underscores(self):
        assert _to_ensg('ORPHAN_NA_TRANS', 'reaction_id', 9606) == 'ORPHAN_NA_TRANS'

    def test_unknown_id_type_returns_none(self):
        result = _to_ensg('X', 'unknown_type', 9606)
        assert result is None


# ---------------------------------------------------------------------------
# translate_pkn — fast cases using chebi + ensembl/reaction_id id_types
# ---------------------------------------------------------------------------

class TestTranslatePkn:

    def _make_met_protein_row(self, source='CHEBI:1', target='ENSG001',
                               id_type_a='chebi', id_type_b='ensembl'):
        return Interaction(
            source=source,
            target=target,
            source_type='small_molecule',
            target_type='protein',
            id_type_a=id_type_a,
            id_type_b=id_type_b,
            interaction_type='catalysis',
            resource='GEM:Human-GEM',
            mor=0,
        )

    def _make_protein_met_row(self, source='ENSG001', target='CHEBI:1',
                               id_type_a='ensembl', id_type_b='chebi'):
        return Interaction(
            source=source,
            target=target,
            source_type='protein',
            target_type='small_molecule',
            id_type_a=id_type_a,
            id_type_b=id_type_b,
            interaction_type='catalysis',
            resource='GEM:Human-GEM',
            mor=0,
        )

    def test_returns_dataframe(self):
        df = _make_df([self._make_met_protein_row()])
        result = translate_pkn(df)
        assert isinstance(result, pd.DataFrame)

    def test_chebi_source_preserved(self):
        df = _make_df([self._make_met_protein_row(source='CHEBI:30616')])
        result = translate_pkn(df)
        assert result.iloc[0]['source'] == 'CHEBI:30616'

    def test_ensembl_target_preserved(self):
        df = _make_df([self._make_met_protein_row(target='ENSG00000001234')])
        result = translate_pkn(df)
        assert result.iloc[0]['target'] == 'ENSG00000001234'

    def test_id_type_a_updated_to_chebi(self):
        df = _make_df([self._make_met_protein_row()])
        result = translate_pkn(df)
        assert result.iloc[0]['id_type_a'] == 'chebi'

    def test_id_type_b_updated_to_ensg(self):
        df = _make_df([self._make_met_protein_row()])
        result = translate_pkn(df)
        assert result.iloc[0]['id_type_b'] == 'ensg'

    def test_drops_row_when_source_untranslatable(self):
        # pubchem id_type with no mapping available → source becomes None → dropped
        row = self._make_met_protein_row(
            source='NOTACID',
            id_type_a='pubchem',
        )
        df = _make_df([row])
        result = translate_pkn(df)
        assert len(result) == 0

    def test_drops_row_when_target_untranslatable(self):
        row = self._make_met_protein_row(
            target='NOTANID',
            id_type_b='uniprot',
        )
        df = _make_df([row])
        result = translate_pkn(df)
        assert len(result) == 0

    def test_reaction_id_preserved_as_target(self):
        # Orphan edge: small_molecule → reaction_id
        row = Interaction(
            source='CHEBI:30616',
            target='MAR04831',
            source_type='small_molecule',
            target_type='protein',
            id_type_a='chebi',
            id_type_b='reaction_id',
            interaction_type='catalysis',
            resource='GEM:Human-GEM',
            mor=0,
            attrs={'orphan': True, 'reaction_id': 'MAR04831'},
        )
        df = _make_df([row])
        result = translate_pkn(df)
        assert len(result) == 1
        assert result.iloc[0]['target'] == 'MAR04831'
        assert result.iloc[0]['id_type_b'] == 'reaction_id'

    def test_reaction_id_preserved_as_source(self):
        # Orphan edge: reaction_id → small_molecule
        row = Interaction(
            source='MAR04831',
            target='CHEBI:30616',
            source_type='protein',
            target_type='small_molecule',
            id_type_a='reaction_id',
            id_type_b='chebi',
            interaction_type='catalysis',
            resource='GEM:Human-GEM',
            mor=0,
            attrs={'orphan': True, 'reaction_id': 'MAR04831'},
        )
        df = _make_df([row])
        result = translate_pkn(df)
        assert len(result) == 1
        assert result.iloc[0]['source'] == 'MAR04831'
        assert result.iloc[0]['id_type_a'] == 'reaction_id'

    def test_mixed_normal_and_orphan_rows(self):
        normal = self._make_met_protein_row(
            source='CHEBI:30616', target='ENSG00000001234',
        )
        orphan = Interaction(
            source='CHEBI:30616',
            target='MAR04831',
            source_type='small_molecule',
            target_type='protein',
            id_type_a='chebi',
            id_type_b='reaction_id',
            interaction_type='catalysis',
            resource='GEM:Human-GEM',
            mor=0,
            attrs={'orphan': True},
        )
        df = _make_df([normal, orphan])
        result = translate_pkn(df)
        assert len(result) == 2
        assert result[result['id_type_b'] == 'ensg'].iloc[0]['target'] == 'ENSG00000001234'
        assert result[result['id_type_b'] == 'reaction_id'].iloc[0]['target'] == 'MAR04831'

    def test_direction_aware_gem_protein_source(self):
        # GEM produces protein-source edges (enzyme → metabolite); ensure
        # translate_pkn handles both directions.
        row = self._make_protein_met_row(
            source='ENSG00000001234', target='CHEBI:30616',
        )
        df = _make_df([row])
        result = translate_pkn(df)
        assert len(result) == 1
        assert result.iloc[0]['id_type_a'] == 'ensg'
        assert result.iloc[0]['id_type_b'] == 'chebi'

    def test_index_reset(self):
        rows = [self._make_met_protein_row() for _ in range(3)]
        df = _make_df(rows)
        result = translate_pkn(df)
        assert list(result.index) == list(range(len(result)))


# ---------------------------------------------------------------------------
# HMDB normalisation
# ---------------------------------------------------------------------------

class TestNormaliseHmdb:

    def test_old_format_5digit(self):
        assert _normalise_hmdb('HMDB00001') == 'HMDB0000001'

    def test_old_format_no_leading_zeros(self):
        assert _normalise_hmdb('HMDB00190') == 'HMDB0000190'

    def test_new_format_already_7digit(self):
        assert _normalise_hmdb('HMDB0000001') == 'HMDB0000001'

    def test_idempotent_on_new_format(self):
        val = 'HMDB0000190'
        assert _normalise_hmdb(val) == val

    def test_large_id(self):
        # IDs with 7 digits already (e.g. HMDB0251697)
        assert _normalise_hmdb('HMDB0251697') == 'HMDB0251697'

    def test_no_prefix_passthrough(self):
        assert _normalise_hmdb('CHEBI:12345') == 'CHEBI:12345'

    def test_empty_string_passthrough(self):
        assert _normalise_hmdb('') == ''


# ---------------------------------------------------------------------------
# _to_chebi with hmdb id_type (mocked UniChem)
# ---------------------------------------------------------------------------

class TestToChebiHmdb:

    _FAKE_MAPPING = {
        'HMDB0000001': 'CHEBI:16015',
        'HMDB0000190': 'CHEBI:17289',
    }

    def test_new_format_translates(self):
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
            return_value=self._FAKE_MAPPING,
        ):
            result = _to_chebi('HMDB0000001', 'hmdb')
        assert result == 'CHEBI:16015'

    def test_old_format_normalised_and_translates(self):
        """Old 5-digit HMDB ID is normalised before lookup."""
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
            return_value=self._FAKE_MAPPING,
        ):
            result = _to_chebi('HMDB00190', 'hmdb')
        assert result == 'CHEBI:17289'

    def test_unknown_hmdb_returns_none(self):
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
            return_value=self._FAKE_MAPPING,
        ):
            result = _to_chebi('HMDB9999999', 'hmdb')
        assert result is None

    def test_old_and_new_same_result(self):
        """Old-format and new-format of the same ID produce the same ChEBI."""
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
            return_value=self._FAKE_MAPPING,
        ):
            old = _to_chebi('HMDB00001', 'hmdb')
            new = _to_chebi('HMDB0000001', 'hmdb')
        assert old == new == 'CHEBI:16015'


# ---------------------------------------------------------------------------
# Vectorized translate_pkn — output equivalence on a synthetic DataFrame
# ---------------------------------------------------------------------------

class TestTranslatePknVectorized:
    """Verify that the vectorized translate_pkn produces identical results
    to the expected output on a small synthetic DataFrame.

    Uses only pass-through id_types (chebi, ensembl, reaction_id) so no
    external calls are needed.
    """

    def _make_row(self, source, target, id_type_a, id_type_b,
                  source_type='small_molecule', target_type='protein',
                  resource='GEM:Human-GEM', interaction_type='catalysis',
                  **attrs_kw):
        return Interaction(
            source=source, target=target,
            source_type=source_type, target_type=target_type,
            id_type_a=id_type_a, id_type_b=id_type_b,
            interaction_type=interaction_type, resource=resource,
            mor=0, attrs=dict(attrs_kw),
        )

    def test_chebi_ensembl_passthrough_unchanged(self):
        """chebi + ensembl rows are returned with values identical to input."""
        rows = [
            self._make_row('CHEBI:30616', 'ENSG00000001234', 'chebi', 'ensembl'),
            self._make_row('CHEBI:15422', 'ENSG00000141510', 'chebi', 'ensembl'),
        ]
        df = _make_df(rows)
        result = translate_pkn(df)
        assert list(result['source']) == ['CHEBI:30616', 'CHEBI:15422']
        assert list(result['target']) == ['ENSG00000001234', 'ENSG00000141510']

    def test_multiple_id_types_in_one_call(self):
        """DataFrame with mixed id_types handled correctly in a single call."""
        rows = [
            # chebi + ensembl: pass-through
            self._make_row('CHEBI:30616', 'ENSG00000001234', 'chebi', 'ensembl'),
            # chebi + reaction_id: orphan, pass-through
            self._make_row('CHEBI:30616', 'MAR04831', 'chebi', 'reaction_id',
                           target_type='protein', orphan=True),
            # protein-source GEM edge: ensembl + chebi
            self._make_row('ENSG00000001234', 'CHEBI:30616', 'ensembl', 'chebi',
                           source_type='protein', target_type='small_molecule'),
        ]
        df = _make_df(rows)
        result = translate_pkn(df)
        assert len(result) == 3
        # chebi+ensembl row
        row0 = result[result['id_type_a'] == 'chebi'].iloc[0]
        assert row0['source'] == 'CHEBI:30616'
        assert row0['id_type_b'] == 'ensg'
        # orphan row — reaction_id preserved
        orphan_rows = result[result['id_type_b'] == 'reaction_id']
        assert len(orphan_rows) == 1
        assert orphan_rows.iloc[0]['target'] == 'MAR04831'

    def test_pubchem_no_mapping_drops_row(self):
        """A pubchem ID with no mapping in the dict → row is dropped."""
        rows = [
            self._make_row('CHEBI:30616', 'ENSG001', 'chebi', 'ensembl'),
            self._make_row('9999999999', 'ENSG001', 'pubchem', 'ensembl'),
        ]
        df = _make_df(rows)
        result = translate_pkn(df)
        # Only the chebi row survives; the pubchem one is dropped.
        assert len(result) == 1
        assert result.iloc[0]['source'] == 'CHEBI:30616'

    def test_result_index_is_contiguous(self):
        """Index is always reset to 0..n-1 after translation."""
        rows = [
            self._make_row('CHEBI:30616', 'ENSG001', 'chebi', 'ensembl'),
            self._make_row('CHEBI:15422', 'ENSG002', 'chebi', 'ensembl'),
            self._make_row('CHEBI:57540', 'ENSG003', 'chebi', 'ensembl'),
        ]
        df = _make_df(rows)
        result = translate_pkn(df)
        assert list(result.index) == list(range(len(result)))
