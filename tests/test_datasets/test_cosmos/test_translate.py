#!/usr/bin/env python

"""Unit tests for COSMOS PKN ID translation.

Fast tests cover the pass-through cases (chebi, ensembl, reaction_id)
that require no external data.  Tests that exercise pypath BioMart
mappings are marked slow.
"""

import pandas as pd
import pytest

from omnipath_metabo.datasets.cosmos._record import Interaction
from omnipath_metabo.datasets.cosmos._translate import (
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
