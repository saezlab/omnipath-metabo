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
    _ensg_to_uniprot_rest,
    _ensp_to_uniprot_rest,
    _entrez_to_uniprot_rest,
    _lipidmaps_to_chebi,
    _metatlas_to_chebi,
    _normalise_hmdb,
    _to_chebi,
    _to_uniprot,
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
# _to_uniprot dispatch
# ---------------------------------------------------------------------------

class TestToUniprot:

    def test_uniprot_passthrough(self):
        assert _to_uniprot('P00533', 'uniprot', 9606) == 'P00533'

    def test_reaction_id_passthrough(self):
        # Orphan pseudo-enzyme IDs must pass through unchanged.
        assert _to_uniprot('MAR04831', 'reaction_id', 9606) == 'MAR04831'

    def test_reaction_id_with_underscores(self):
        assert _to_uniprot('ORPHAN_NA_TRANS', 'reaction_id', 9606) == 'ORPHAN_NA_TRANS'

    def test_unknown_id_type_returns_none(self):
        result = _to_uniprot('X', 'unknown_type', 9606)
        assert result is None


# ---------------------------------------------------------------------------
# translate_pkn — fast cases using chebi + ensembl/reaction_id id_types
# ---------------------------------------------------------------------------

class TestTranslatePkn:

    def _make_met_protein_row(self, source='CHEBI:1', target='P12345',
                               id_type_a='chebi', id_type_b='uniprot'):
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

    def _make_protein_met_row(self, source='P12345', target='CHEBI:1',
                               id_type_a='uniprot', id_type_b='chebi'):
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

    def test_uniprot_target_preserved(self):
        df = _make_df([self._make_met_protein_row(target='P00533')])
        result = translate_pkn(df)
        assert result.iloc[0]['target'] == 'P00533'

    def test_id_type_a_updated_to_chebi(self):
        df = _make_df([self._make_met_protein_row()])
        result = translate_pkn(df)
        assert result.iloc[0]['id_type_a'] == 'chebi'

    def test_id_type_b_updated_to_uniprot(self):
        df = _make_df([self._make_met_protein_row()])
        result = translate_pkn(df)
        assert result.iloc[0]['id_type_b'] == 'uniprot'

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
        # ensp with a fake ID: pypath returns empty set with unloaded tables → dropped
        row = self._make_met_protein_row(
            target='ENSP999FAKE',
            id_type_b='ensp',
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
            source='CHEBI:30616', target='P00533',
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
        assert result[result['id_type_b'] == 'uniprot'].iloc[0]['target'] == 'P00533'
        assert result[result['id_type_b'] == 'reaction_id'].iloc[0]['target'] == 'MAR04831'

    def test_direction_aware_gem_protein_source(self):
        # GEM produces protein-source edges (enzyme → metabolite); ensure
        # translate_pkn handles both directions.
        row = self._make_protein_met_row(
            source='P00533', target='CHEBI:30616',
        )
        df = _make_df([row])
        result = translate_pkn(df)
        assert len(result) == 1
        assert result.iloc[0]['id_type_a'] == 'uniprot'
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

    def test_chebi_uniprot_passthrough_unchanged(self):
        """chebi + uniprot rows are returned with values identical to input."""
        rows = [
            self._make_row('CHEBI:30616', 'P00533', 'chebi', 'uniprot'),
            self._make_row('CHEBI:15422', 'P04637', 'chebi', 'uniprot'),
        ]
        df = _make_df(rows)
        result = translate_pkn(df)
        assert list(result['source']) == ['CHEBI:30616', 'CHEBI:15422']
        assert list(result['target']) == ['P00533', 'P04637']

    def test_multiple_id_types_in_one_call(self):
        """DataFrame with mixed id_types handled correctly in a single call."""
        rows = [
            # chebi + uniprot: pass-through
            self._make_row('CHEBI:30616', 'P00533', 'chebi', 'uniprot'),
            # chebi + reaction_id: orphan, pass-through
            self._make_row('CHEBI:30616', 'MAR04831', 'chebi', 'reaction_id',
                           target_type='protein', orphan=True),
            # protein-source edge: uniprot + chebi
            self._make_row('P00533', 'CHEBI:30616', 'uniprot', 'chebi',
                           source_type='protein', target_type='small_molecule'),
        ]
        df = _make_df(rows)
        result = translate_pkn(df)
        assert len(result) == 3
        # chebi+uniprot row
        row0 = result[result['id_type_a'] == 'chebi'].iloc[0]
        assert row0['source'] == 'CHEBI:30616'
        assert row0['id_type_b'] == 'uniprot'
        # orphan row — reaction_id preserved
        orphan_rows = result[result['id_type_b'] == 'reaction_id']
        assert len(orphan_rows) == 1
        assert orphan_rows.iloc[0]['target'] == 'MAR04831'

    def test_pubchem_no_mapping_drops_row(self):
        """A pubchem ID with no mapping in the dict → row is dropped."""
        rows = [
            self._make_row('CHEBI:30616', 'P00533', 'chebi', 'uniprot'),
            self._make_row('9999999999', 'P00533', 'pubchem', 'uniprot'),
        ]
        df = _make_df(rows)
        result = translate_pkn(df)
        # Only the chebi row survives; the pubchem one is dropped.
        assert len(result) == 1
        assert result.iloc[0]['source'] == 'CHEBI:30616'

    def test_result_index_is_contiguous(self):
        """Index is always reset to 0..n-1 after translation."""
        rows = [
            self._make_row('CHEBI:30616', 'P00533', 'chebi', 'uniprot'),
            self._make_row('CHEBI:15422', 'P04637', 'chebi', 'uniprot'),
            self._make_row('CHEBI:57540', 'P12345', 'chebi', 'uniprot'),
        ]
        df = _make_df(rows)
        result = translate_pkn(df)
        assert list(result.index) == list(range(len(result)))


# ---------------------------------------------------------------------------
# translate_pkn vectorized path — id_type='hmdb' (mocked UniChem)
# ---------------------------------------------------------------------------

class TestTranslatePknHmdb:
    """Verify translate_pkn normalises old HMDB IDs before the dict lookup."""

    _FAKE_MAPPING = {
        'HMDB0000001': 'CHEBI:16015',
        'HMDB0000190': 'CHEBI:17289',
    }

    def _make_row(self, source, target='P00533', id_type_a='hmdb',
                  id_type_b='uniprot'):
        return Interaction(
            source=source,
            target=target,
            source_type='small_molecule',
            target_type='protein',
            id_type_a=id_type_a,
            id_type_b=id_type_b,
            interaction_type='catalysis',
            resource='STITCH',
            mor=1,
        )

    def test_new_format_hmdb_translates(self):
        """7-digit HMDB ID is looked up directly."""
        rows = [self._make_row('HMDB0000001')]
        df = _make_df(rows)
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
            return_value=self._FAKE_MAPPING,
        ):
            result = translate_pkn(df)
        assert len(result) == 1
        assert result.iloc[0]['source'] == 'CHEBI:16015'

    def test_old_format_hmdb_normalised_then_translates(self):
        """5-digit HMDB ID is normalised to 7-digit before lookup."""
        rows = [self._make_row('HMDB00190')]
        df = _make_df(rows)
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
            return_value=self._FAKE_MAPPING,
        ):
            result = translate_pkn(df)
        assert len(result) == 1
        assert result.iloc[0]['source'] == 'CHEBI:17289'

    def test_old_and_new_same_id_same_output(self):
        """Old-format and new-format of the same HMDB ID produce the same ChEBI."""
        rows = [self._make_row('HMDB00001'), self._make_row('HMDB0000001')]
        df = _make_df(rows)
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
            return_value=self._FAKE_MAPPING,
        ):
            result = translate_pkn(df)
        assert len(result) == 2
        assert (result['source'] == 'CHEBI:16015').all()

    def test_unmapped_hmdb_drops_row(self):
        """HMDB ID with no ChEBI mapping causes the row to be dropped."""
        rows = [
            self._make_row('HMDB9999999'),   # not in mapping → dropped
            self._make_row('HMDB0000001'),   # maps → kept
        ]
        df = _make_df(rows)
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
            return_value=self._FAKE_MAPPING,
        ):
            result = translate_pkn(df)
        assert len(result) == 1
        assert result.iloc[0]['source'] == 'CHEBI:16015'

    def test_id_type_a_updated_to_chebi(self):
        """After translation, id_type_a is updated from 'hmdb' to 'chebi'."""
        rows = [self._make_row('HMDB0000001')]
        df = _make_df(rows)
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
            return_value=self._FAKE_MAPPING,
        ):
            result = translate_pkn(df)
        assert result.iloc[0]['id_type_a'] == 'chebi'


# ---------------------------------------------------------------------------
# _lipidmaps_to_chebi (mocked UniChem)
# ---------------------------------------------------------------------------

class TestLipidmapsToChebi:

    _FAKE_MAPPING = {
        'LMFA01010001': {'CHEBI:15756'},
        'LMFA01010002': {'CHEBI:15904'},
    }

    def test_returns_dict(self):
        with patch(
            'pypath.inputs.unichem.unichem_mapping',
            return_value=self._FAKE_MAPPING,
        ):
            _lipidmaps_to_chebi.cache_clear()
            result = _lipidmaps_to_chebi()
        assert isinstance(result, dict)

    def test_known_id_resolves(self):
        with patch(
            'pypath.inputs.unichem.unichem_mapping',
            return_value=self._FAKE_MAPPING,
        ):
            _lipidmaps_to_chebi.cache_clear()
            result = _lipidmaps_to_chebi()
        assert result['LMFA01010001'] == 'CHEBI:15756'

    def test_all_entries_present(self):
        with patch(
            'pypath.inputs.unichem.unichem_mapping',
            return_value=self._FAKE_MAPPING,
        ):
            _lipidmaps_to_chebi.cache_clear()
            result = _lipidmaps_to_chebi()
        assert len(result) == 2

    def test_unknown_id_absent(self):
        with patch(
            'pypath.inputs.unichem.unichem_mapping',
            return_value=self._FAKE_MAPPING,
        ):
            _lipidmaps_to_chebi.cache_clear()
            result = _lipidmaps_to_chebi()
        assert 'LMFA_NOTEXIST' not in result


# ---------------------------------------------------------------------------
# _metatlas_to_chebi fallback chain (mocked external sources)
# ---------------------------------------------------------------------------

class TestMetatlasToChebiFallbacks:
    """Verify the four-step fallback chain in _metatlas_to_chebi.

    Uses mocked external sources so no network calls are made.
    Each test covers a specific fallback step being the one that resolves
    an ID that has no direct metChEBIID.
    """

    # Fake GEM rows: one with direct ChEBI, one per fallback step
    _FAKE_ROWS = [
        # direct ChEBI
        {'metsNoComp': 'MAM00001', 'metChEBIID': 'CHEBI:111',
         'metMetaNetXID': '', 'metLipidMapsID': '', 'metPubChemID': '', 'metKEGGID': '', 'metHMDBID': ''},
        # resolved via MetaNetX
        {'metsNoComp': 'MAM00002', 'metChEBIID': '',
         'metMetaNetXID': 'MNXM999', 'metLipidMapsID': '', 'metPubChemID': '', 'metKEGGID': '', 'metHMDBID': ''},
        # resolved via LipidMaps
        {'metsNoComp': 'MAM00003', 'metChEBIID': '',
         'metMetaNetXID': '', 'metLipidMapsID': 'LMFA01010001', 'metPubChemID': '', 'metKEGGID': '', 'metHMDBID': ''},
        # resolved via PubChem
        {'metsNoComp': 'MAM00004', 'metChEBIID': '',
         'metMetaNetXID': '', 'metLipidMapsID': '', 'metPubChemID': '5793', 'metKEGGID': '', 'metHMDBID': ''},
        # resolved via KEGG
        {'metsNoComp': 'MAM00005', 'metChEBIID': '',
         'metMetaNetXID': '', 'metLipidMapsID': '', 'metPubChemID': '', 'metKEGGID': 'C00001', 'metHMDBID': ''},
        # resolved via HMDB
        {'metsNoComp': 'MAM00006', 'metChEBIID': '',
         'metMetaNetXID': '', 'metLipidMapsID': '', 'metPubChemID': '', 'metKEGGID': '', 'metHMDBID': 'HMDB0000001'},
        # unresolvable — no IDs
        {'metsNoComp': 'MAM00007', 'metChEBIID': '',
         'metMetaNetXID': '', 'metLipidMapsID': '', 'metPubChemID': '', 'metKEGGID': '', 'metHMDBID': ''},
    ]

    def _run(self):
        """Run _metatlas_to_chebi with all external sources mocked."""
        _metatlas_to_chebi.cache_clear()
        with (
            patch(
                'pypath.inputs.metatlas._gem.metatlas_gem_metabolites',
                return_value=self._FAKE_ROWS,
            ),
            patch(
                'pypath.inputs.metanetx.metanetx_metabolite_chebi',
                return_value={'MNXM999': 'CHEBI:222'},
            ),
            patch(
                'omnipath_metabo.datasets.cosmos._translate._lipidmaps_to_chebi',
                return_value={'LMFA01010001': 'CHEBI:333'},
            ),
            patch(
                'omnipath_metabo.datasets.cosmos._translate._pubchem_to_chebi',
                return_value={'5793': 'CHEBI:444'},
            ),
            patch(
                'omnipath_metabo.datasets.cosmos._translate._pubchem_to_chebi_ramp',
                return_value={},
            ),
            patch(
                'pypath.inputs.kegg.kegg_compound_chebi',
                return_value={'C00001': 'CHEBI:15377'},
            ),
            patch(
                'omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi',
                return_value={'HMDB0000001': 'CHEBI:555'},
            ),
        ):
            return _metatlas_to_chebi('FakeGEM')

    def test_direct_chebi_present(self):
        result = self._run()
        assert result['MAM00001'] == 'CHEBI:111'

    def test_metanetx_fallback(self):
        result = self._run()
        assert result['MAM00002'] == 'CHEBI:222'

    def test_lipidmaps_fallback(self):
        result = self._run()
        assert result['MAM00003'] == 'CHEBI:333'

    def test_pubchem_fallback(self):
        result = self._run()
        assert result['MAM00004'] == 'CHEBI:444'

    def test_kegg_fallback(self):
        result = self._run()
        assert result['MAM00005'] == 'CHEBI:15377'

    def test_hmdb_fallback(self):
        result = self._run()
        assert result['MAM00006'] == 'CHEBI:555'

    def test_unresolvable_absent(self):
        result = self._run()
        assert 'MAM00007' not in result

    def test_total_resolved_count(self):
        result = self._run()
        assert len(result) == 6

    def test_semicolon_separated_metanetx(self):
        """If metMetaNetXID has multiple values (';'-separated), all are tried."""
        rows = [
            {'metsNoComp': 'MAM00010', 'metChEBIID': '',
             'metMetaNetXID': 'MNXM_MISS;MNXM999',
             'metLipidMapsID': '', 'metPubChemID': '', 'metKEGGID': '', 'metHMDBID': ''},
        ]
        _metatlas_to_chebi.cache_clear()
        with (
            patch(
                'pypath.inputs.metatlas._gem.metatlas_gem_metabolites',
                return_value=rows,
            ),
            patch(
                'pypath.inputs.metanetx.metanetx_metabolite_chebi',
                return_value={'MNXM999': 'CHEBI:222'},
            ),
            patch('omnipath_metabo.datasets.cosmos._translate._lipidmaps_to_chebi', return_value={}),
            patch('omnipath_metabo.datasets.cosmos._translate._pubchem_to_chebi', return_value={}),
            patch('omnipath_metabo.datasets.cosmos._translate._pubchem_to_chebi_ramp', return_value={}),
            patch('pypath.inputs.kegg.kegg_compound_chebi', return_value={}),
            patch('omnipath_metabo.datasets.cosmos._translate._hmdb_to_chebi', return_value={}),
        ):
            result = _metatlas_to_chebi('FakeGEM')
        assert result['MAM00010'] == 'CHEBI:222'


# ---------------------------------------------------------------------------
# _ensp_to_uniprot_rest — mocked HTTP
# ---------------------------------------------------------------------------

class TestEnspToUniprotRest:
    """Unit tests for _ensp_to_uniprot_rest using mocked HTTP responses."""

    def _mock_responses(self, results_payload, monkeypatch):
        """
        Wire up requests.post (submit) and requests.get (status + results).
        """
        import requests as _requests

        class _FakeSubmitResp:
            status_code = 200
            def raise_for_status(self): pass
            def json(self): return {'jobId': 'TESTJOB'}

        class _FakeStatusResp:
            status_code = 200
            def json(self): return {'jobStatus': 'FINISHED'}

        class _FakeResultsResp:
            status_code = 200
            def raise_for_status(self): pass
            def json(self): return results_payload

        call_count = {'get': 0}

        def _fake_post(url, data=None, timeout=None):
            return _FakeSubmitResp()

        def _fake_get(url, timeout=None, allow_redirects=True):
            call_count['get'] += 1
            if 'status' in url:
                return _FakeStatusResp()
            return _FakeResultsResp()

        monkeypatch.setattr(_requests, 'post', _fake_post)
        monkeypatch.setattr(_requests, 'get', _fake_get)

    def test_returns_dict(self, monkeypatch):
        self._mock_responses({'results': [], 'failedIds': ['ENSP0001']}, monkeypatch)
        result = _ensp_to_uniprot_rest(['ENSP0001'])
        assert isinstance(result, dict)

    def test_resolved_ensp_present(self, monkeypatch):
        self._mock_responses({
            'results': [{'from': 'ENSP0001', 'to': 'P12345'}],
            'failedIds': [],
        }, monkeypatch)
        result = _ensp_to_uniprot_rest(['ENSP0001'])
        assert result['ENSP0001'] == 'P12345'

    def test_failed_ensp_absent(self, monkeypatch):
        self._mock_responses({
            'results': [],
            'failedIds': ['ENSP0002'],
        }, monkeypatch)
        result = _ensp_to_uniprot_rest(['ENSP0002'])
        assert 'ENSP0002' not in result

    def test_mixed_resolved_and_failed(self, monkeypatch):
        self._mock_responses({
            'results': [{'from': 'ENSP0001', 'to': 'P12345'}],
            'failedIds': ['ENSP0002'],
        }, monkeypatch)
        result = _ensp_to_uniprot_rest(['ENSP0001', 'ENSP0002'])
        assert result['ENSP0001'] == 'P12345'
        assert 'ENSP0002' not in result

    def test_to_as_dict_with_primary_accession(self, monkeypatch):
        """UniProt sometimes returns 'to' as a dict with primaryAccession."""
        self._mock_responses({
            'results': [{'from': 'ENSP0003', 'to': {'primaryAccession': 'Q99999'}}],
            'failedIds': [],
        }, monkeypatch)
        result = _ensp_to_uniprot_rest(['ENSP0003'])
        assert result['ENSP0003'] == 'Q99999'

    def test_empty_input_returns_empty(self, monkeypatch):
        # Should not make any HTTP calls
        result = _ensp_to_uniprot_rest([])
        assert result == {}

    def test_submit_failure_returns_empty(self, monkeypatch):
        import requests as _requests

        def _bad_post(url, data=None, timeout=None):
            raise ConnectionError('network down')

        monkeypatch.setattr(_requests, 'post', _bad_post)
        result = _ensp_to_uniprot_rest(['ENSP0001'])
        assert result == {}


# ---------------------------------------------------------------------------
# _ensg_to_uniprot_rest — delegates to _uniprot_idmap_batch
# ---------------------------------------------------------------------------

class TestEnsgToUniprotRest:
    """_ensg_to_uniprot_rest is a thin wrapper — just verify delegation."""

    def test_returns_dict(self, monkeypatch):
        import requests as _requests

        class _FakeSubmit:
            status_code = 200
            def raise_for_status(self): pass
            def json(self): return {'jobId': 'JOB2'}

        class _FakeStatus:
            status_code = 200
            def json(self): return {'jobStatus': 'FINISHED'}

        class _FakeResults:
            status_code = 200
            def raise_for_status(self): pass
            def json(self): return {
                'results': [{'from': 'ENSG0001', 'to': 'Q11111'}],
                'failedIds': [],
            }

        def _fake_post(url, data=None, timeout=None):
            assert data.get('from') == 'Ensembl', 'wrong from_db'
            return _FakeSubmit()

        def _fake_get(url, timeout=None, allow_redirects=True):
            if 'status' in url:
                return _FakeStatus()
            return _FakeResults()

        monkeypatch.setattr(_requests, 'post', _fake_post)
        monkeypatch.setattr(_requests, 'get', _fake_get)

        result = _ensg_to_uniprot_rest(['ENSG0001'])
        assert result == {'ENSG0001': 'Q11111'}

    def test_empty_input_returns_empty(self, monkeypatch):
        result = _ensg_to_uniprot_rest([])
        assert result == {}


# ---------------------------------------------------------------------------
# _entrez_to_uniprot_rest — delegates to _uniprot_idmap_batch with 'GeneID'
# ---------------------------------------------------------------------------

class TestEntrezToUniprotRest:
    """_entrez_to_uniprot_rest is a thin wrapper — verify delegation and basic flow."""

    def _mock_responses(self, results_payload, monkeypatch):
        import requests as _requests

        class _FakeSubmit:
            status_code = 200
            def raise_for_status(self): pass
            def json(self): return {'jobId': 'ENTREZJOB'}

        class _FakeStatus:
            status_code = 200
            def json(self): return {'jobStatus': 'FINISHED'}

        class _FakeResults:
            status_code = 200
            def raise_for_status(self): pass
            def json(self): return results_payload

        def _fake_post(url, data=None, timeout=None):
            assert data.get('from') == 'GeneID', 'wrong from_db'
            return _FakeSubmit()

        def _fake_get(url, timeout=None, allow_redirects=True):
            if 'status' in url:
                return _FakeStatus()
            return _FakeResults()

        monkeypatch.setattr(_requests, 'post', _fake_post)
        monkeypatch.setattr(_requests, 'get', _fake_get)

    def test_uses_gene_id_from_db(self, monkeypatch):
        """Verifies 'from=GeneID' is passed to the UniProt API."""
        self._mock_responses({'results': [], 'failedIds': ['1591']}, monkeypatch)
        result = _entrez_to_uniprot_rest(['1591'])
        assert isinstance(result, dict)

    def test_resolved_entrez_present(self, monkeypatch):
        self._mock_responses({
            'results': [{'from': '1591', 'to': 'P05108'}],
            'failedIds': [],
        }, monkeypatch)
        result = _entrez_to_uniprot_rest(['1591'])
        assert result['1591'] == 'P05108'

    def test_failed_entrez_absent(self, monkeypatch):
        self._mock_responses({'results': [], 'failedIds': ['9999999']}, monkeypatch)
        result = _entrez_to_uniprot_rest(['9999999'])
        assert '9999999' not in result

    def test_empty_input_returns_empty(self, monkeypatch):
        result = _entrez_to_uniprot_rest([])
        assert result == {}


# ---------------------------------------------------------------------------
# ensembl compound-ID split logic in translate_pkn
# ---------------------------------------------------------------------------

class TestEnsemblCompoundIdSplit:
    """Verify that _-joined compound ENSG IDs are split and resolved."""

    def _make_gem_df(self, ensg_id):
        """One GEM metabolic edge with a compound ENSG protein target."""
        from omnipath_metabo.datasets.cosmos._record import Interaction
        row = Interaction(
            source='MAM00001',
            target=ensg_id,
            source_type='small_molecule',
            target_type='protein',
            id_type_a='metatlas',
            id_type_b='ensembl',
            interaction_type='reaction',
            resource='GEM:Human-GEM',
            mor=1,
            locations=['c'],
            attrs={},
        )
        return pd.DataFrame([row], columns=list(Interaction._fields))

    def test_compound_ensg_split_resolves(self, monkeypatch):
        """A _-joined ENSG that fails as a whole resolves via first component."""
        import pypath.utils.mapping as mapping_mod

        compound = 'ENSG00000001_ENSG00000002'

        def _fake_map(uid, src, tgt, ncbi_tax_id=9606):
            # Whole compound ID fails; first component resolves
            if uid == 'ENSG00000001':
                return {'P99999'}
            return set()

        monkeypatch.setattr(mapping_mod, 'map_name', _fake_map)

        # Also mock metatlas→ChEBI so the metabolite side passes through
        import omnipath_metabo.datasets.cosmos._translate as tr
        monkeypatch.setattr(tr, '_metatlas_to_chebi', lambda gem: {'MAM00001': 'CHEBI:1'})

        df = self._make_gem_df(compound)
        result = translate_pkn(df, organism=9606)
        assert len(result) == 1
        assert result.iloc[0]['target'] == frozenset({'P99999'})

    def test_single_ensg_no_split(self, monkeypatch):
        """A single ENSG that resolves directly is not split."""
        import pypath.utils.mapping as mapping_mod

        def _fake_map(uid, src, tgt, ncbi_tax_id=9606):
            if uid == 'ENSG00000001':
                return {'P11111'}
            return set()

        monkeypatch.setattr(mapping_mod, 'map_name', _fake_map)

        import omnipath_metabo.datasets.cosmos._translate as tr
        monkeypatch.setattr(tr, '_metatlas_to_chebi', lambda gem: {'MAM00001': 'CHEBI:1'})

        df = self._make_gem_df('ENSG00000001')
        result = translate_pkn(df, organism=9606)
        assert len(result) == 1
        assert result.iloc[0]['target'] == frozenset({'P11111'})

    def test_unresolvable_compound_drops_row(self, monkeypatch):
        """A compound ENSG where all components fail → row dropped."""
        import pypath.utils.mapping as mapping_mod
        import omnipath_metabo.datasets.cosmos._translate as tr

        monkeypatch.setattr(mapping_mod, 'map_name', lambda *a, **kw: set())
        monkeypatch.setattr(tr, '_ensg_to_uniprot_rest', lambda ids: {})
        monkeypatch.setattr(tr, '_metatlas_to_chebi', lambda gem: {'MAM00001': 'CHEBI:1'})

        df = self._make_gem_df('ENSG00000999_ENSG00000998')
        result = translate_pkn(df, organism=9606)
        assert len(result) == 0
