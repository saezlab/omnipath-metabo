#!/usr/bin/env python

"""Unit tests for COSMOS PKN ID translation.

Fast tests cover the pass-through cases (chebi, ensembl, reaction_id)
that require no external data.  Tests that exercise the omnipath-utils
adapter mappings are mocked.
"""

import pandas as pd
import pytest

from unittest.mock import MagicMock, patch

from omnipath_metabo.datasets.cosmos._record import Interaction
from omnipath_metabo.datasets.cosmos._translate import (
    _build_metab_mapping,
    _build_protein_mapping,
    _hmdb_synonyms_chebi,
    _looks_like_chemical_name,
    _metatlas_to_chebi,
    _name_to_chebi,
    _ramp_synonyms_chebi,
    translate_pkn,
)


def _make_df(rows):
    """Build a PKN DataFrame from a list of Interaction namedtuples."""

    return pd.DataFrame(rows, columns=list(Interaction._fields))


# ---------------------------------------------------------------------------
# _looks_like_chemical_name
# ---------------------------------------------------------------------------

class TestLooksLikeChemicalName:

    def test_simple_name(self):
        assert _looks_like_chemical_name('ATP') is True

    def test_percentage_rejected(self):
        assert _looks_like_chemical_name('10.2% residual activity') is False

    def test_concentration_rejected(self):
        assert _looks_like_chemical_name('100 nM compound') is False

    def test_assay_vocabulary_rejected(self):
        assert _looks_like_chemical_name('protein signaling pathway') is False

    def test_long_name_rejected(self):
        assert _looks_like_chemical_name('x' * 200) is False


# ---------------------------------------------------------------------------
# Adapter tests (_mapping.py)
# ---------------------------------------------------------------------------

class TestMappingAdapter:
    """Verify the dual-mode adapter imports and exports the expected symbols."""

    def test_mapping_translate_callable(self):
        from omnipath_metabo.datasets.cosmos._mapping import mapping_translate
        assert callable(mapping_translate)

    def test_mapping_table_callable(self):
        from omnipath_metabo.datasets.cosmos._mapping import mapping_table
        assert callable(mapping_table)

    def test_mapping_mode_is_string(self):
        from omnipath_metabo.datasets.cosmos._mapping import mapping_mode
        assert isinstance(mapping_mode, str)
        assert mapping_mode in ('database', 'http')


# ---------------------------------------------------------------------------
# Metabolite mapping table functions (mocked adapter)
# ---------------------------------------------------------------------------

class TestPubchemChebiTable:

    def test_returns_dict(self):
        from omnipath_metabo.datasets.cosmos._translate import _pubchem_chebi_table
        _pubchem_chebi_table.cache_clear()
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_table',
            return_value={'5793': {'CHEBI:18367'}, '2244': {'CHEBI:15365'}},
        ):
            result = _pubchem_chebi_table()
        assert isinstance(result, dict)
        assert result['5793'] == 'CHEBI:18367'
        assert result['2244'] == 'CHEBI:15365'


class TestLipidmapsChebiTable:

    def test_returns_dict(self):
        from omnipath_metabo.datasets.cosmos._translate import _lipidmaps_chebi_table
        _lipidmaps_chebi_table.cache_clear()
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_table',
            return_value={
                'LMFA01010001': {'CHEBI:15756'},
                'LMFA01010002': {'CHEBI:15904'},
            },
        ):
            result = _lipidmaps_chebi_table()
        assert isinstance(result, dict)
        assert result['LMFA01010001'] == 'CHEBI:15756'
        assert len(result) == 2

    def test_unknown_id_absent(self):
        from omnipath_metabo.datasets.cosmos._translate import _lipidmaps_chebi_table
        _lipidmaps_chebi_table.cache_clear()
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_table',
            return_value={'LMFA01010001': {'CHEBI:15756'}},
        ):
            result = _lipidmaps_chebi_table()
        assert 'LMFA_NOTEXIST' not in result


class TestHmdbChebiTable:

    def test_returns_dict(self):
        from omnipath_metabo.datasets.cosmos._translate import _hmdb_chebi_table
        _hmdb_chebi_table.cache_clear()
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_table',
            return_value={
                'HMDB0000001': {'CHEBI:16015'},
                'HMDB0000190': {'CHEBI:17289'},
            },
        ):
            result = _hmdb_chebi_table()
        assert isinstance(result, dict)
        assert result['HMDB0000001'] == 'CHEBI:16015'


class TestBiggChebiTable:

    def test_returns_dict(self):
        from omnipath_metabo.datasets.cosmos._translate import _bigg_chebi_table
        _bigg_chebi_table.cache_clear()
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_table',
            return_value={'glc__D': {'CHEBI:17634'}},
        ):
            result = _bigg_chebi_table()
        assert isinstance(result, dict)
        assert result['glc__D'] == 'CHEBI:17634'


class TestKeggChebiTable:

    def test_returns_dict(self):
        from omnipath_metabo.datasets.cosmos._translate import _kegg_chebi_table
        _kegg_chebi_table.cache_clear()
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_table',
            return_value={'C00001': {'CHEBI:15377'}},
        ):
            result = _kegg_chebi_table()
        assert isinstance(result, dict)
        assert result['C00001'] == 'CHEBI:15377'


class TestMetanetxChebiTable:

    def test_returns_dict(self):
        from omnipath_metabo.datasets.cosmos._translate import _metanetx_chebi_table
        _metanetx_chebi_table.cache_clear()
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_table',
            return_value={'MNXM999': {'CHEBI:222'}},
        ):
            result = _metanetx_chebi_table()
        assert isinstance(result, dict)
        assert result['MNXM999'] == 'CHEBI:222'


# ---------------------------------------------------------------------------
# _name_to_chebi, _hmdb_synonyms_chebi, _ramp_synonyms_chebi
# ---------------------------------------------------------------------------

class TestNameToChebi:

    def test_hmdb_hit(self):
        _hmdb_synonyms_chebi.cache_clear()
        _ramp_synonyms_chebi.cache_clear()
        with (
            patch(
                'pypath.inputs.hmdb.metabolites.synonyms_chebi',
                return_value={'atp': 'CHEBI:30616'},
            ),
            patch(
                'pypath.inputs.ramp._mapping.ramp_synonyms_chebi',
                return_value={},
            ),
        ):
            result = _name_to_chebi('ATP')
        assert result == 'CHEBI:30616'

    def test_ramp_hit(self):
        _hmdb_synonyms_chebi.cache_clear()
        _ramp_synonyms_chebi.cache_clear()
        with (
            patch(
                'pypath.inputs.hmdb.metabolites.synonyms_chebi',
                return_value={},
            ),
            patch(
                'pypath.inputs.ramp._mapping.ramp_synonyms_chebi',
                return_value={'nadh': 'CHEBI:57945'},
            ),
        ):
            result = _name_to_chebi('NADH')
        assert result == 'CHEBI:57945'

    def test_returns_none_for_unknown(self):
        _hmdb_synonyms_chebi.cache_clear()
        _ramp_synonyms_chebi.cache_clear()
        with (
            patch(
                'pypath.inputs.hmdb.metabolites.synonyms_chebi',
                return_value={},
            ),
            patch(
                'pypath.inputs.ramp._mapping.ramp_synonyms_chebi',
                return_value={},
            ),
            patch(
                'pypath.inputs.pubchem.pubchem_name_cids',
                return_value=[],
            ),
        ):
            result = _name_to_chebi('totally_unknown_compound_xyz_999')
        assert result is None


class TestHmdbSynonymsChebi:

    def test_returns_dict(self):
        _hmdb_synonyms_chebi.cache_clear()
        with patch(
            'pypath.inputs.hmdb.metabolites.synonyms_chebi',
            return_value={'glucose': 'CHEBI:17234'},
        ):
            result = _hmdb_synonyms_chebi()
        assert isinstance(result, dict)
        assert result['glucose'] == 'CHEBI:17234'

    def test_returns_empty_on_failure(self):
        _hmdb_synonyms_chebi.cache_clear()
        with patch(
            'pypath.inputs.hmdb.metabolites.synonyms_chebi',
            side_effect=ConnectionError('cloudflare block'),
        ):
            result = _hmdb_synonyms_chebi()
        assert result == {}


class TestRampSynonymsChebi:

    def test_returns_dict(self):
        _ramp_synonyms_chebi.cache_clear()
        with patch(
            'pypath.inputs.ramp._mapping.ramp_synonyms_chebi',
            return_value={'lactate': 'CHEBI:24996'},
        ):
            result = _ramp_synonyms_chebi()
        assert isinstance(result, dict)
        assert result['lactate'] == 'CHEBI:24996'

    def test_returns_empty_on_failure(self):
        _ramp_synonyms_chebi.cache_clear()
        with patch(
            'pypath.inputs.ramp._mapping.ramp_synonyms_chebi',
            side_effect=RuntimeError('download failed'),
        ):
            result = _ramp_synonyms_chebi()
        assert result == {}


# ---------------------------------------------------------------------------
# translate_pkn -- fast cases using chebi + ensembl/reaction_id id_types
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
        # pubchem id_type with no mapping available -> source becomes None -> dropped
        row = self._make_met_protein_row(
            source='NOTACID',
            id_type_a='pubchem',
        )
        df = _make_df([row])
        result = translate_pkn(df)
        assert len(result) == 0

    def test_drops_row_when_target_untranslatable(self):
        # ensp with a fake ID: mapping_translate returns empty -> dropped
        row = self._make_met_protein_row(
            target='ENSP999FAKE',
            id_type_b='ensp',
        )
        df = _make_df([row])
        result = translate_pkn(df)
        assert len(result) == 0

    def test_reaction_id_preserved_as_target(self):
        # Orphan edge: small_molecule -> reaction_id
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
        # Orphan edge: reaction_id -> small_molecule
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
        # GEM produces protein-source edges (enzyme -> metabolite); ensure
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
# Vectorized translate_pkn -- output equivalence on a synthetic DataFrame
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
        # orphan row -- reaction_id preserved
        orphan_rows = result[result['id_type_b'] == 'reaction_id']
        assert len(orphan_rows) == 1
        assert orphan_rows.iloc[0]['target'] == 'MAR04831'

    def test_pubchem_no_mapping_drops_row(self):
        """A pubchem ID with no mapping in the dict -> row is dropped."""
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
# translate_pkn vectorized path -- id_type='hmdb' (mocked table)
# ---------------------------------------------------------------------------

class TestTranslatePknHmdb:
    """Verify translate_pkn handles HMDB IDs via the bulk table lookup."""

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

    def test_hmdb_translates(self):
        """HMDB ID is looked up in the _hmdb_chebi_table dict."""
        rows = [self._make_row('HMDB0000001')]
        df = _make_df(rows)
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_chebi_table',
            return_value=self._FAKE_MAPPING,
        ):
            result = translate_pkn(df)
        assert len(result) == 1
        assert result.iloc[0]['source'] == 'CHEBI:16015'

    def test_unmapped_hmdb_drops_row(self):
        """HMDB ID with no ChEBI mapping causes the row to be dropped."""
        rows = [
            self._make_row('HMDB9999999'),   # not in mapping -> dropped
            self._make_row('HMDB0000001'),   # maps -> kept
        ]
        df = _make_df(rows)
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_chebi_table',
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
            'omnipath_metabo.datasets.cosmos._translate._hmdb_chebi_table',
            return_value=self._FAKE_MAPPING,
        ):
            result = translate_pkn(df)
        assert result.iloc[0]['id_type_a'] == 'chebi'


# ---------------------------------------------------------------------------
# _metatlas_to_chebi fallback chain (mocked external sources)
# ---------------------------------------------------------------------------

class TestMetatlasToChebiFallbacks:
    """Verify the six-step fallback chain in _metatlas_to_chebi.

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
        # unresolvable -- no IDs
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
                'omnipath_metabo.datasets.cosmos._translate._metanetx_chebi_table',
                return_value={'MNXM999': 'CHEBI:222'},
            ),
            patch(
                'omnipath_metabo.datasets.cosmos._translate._lipidmaps_chebi_table',
                return_value={'LMFA01010001': 'CHEBI:333'},
            ),
            patch(
                'omnipath_metabo.datasets.cosmos._translate._pubchem_chebi_table',
                return_value={'5793': 'CHEBI:444'},
            ),
            patch(
                'omnipath_metabo.datasets.cosmos._translate._kegg_chebi_table',
                return_value={'C00001': 'CHEBI:15377'},
            ),
            patch(
                'omnipath_metabo.datasets.cosmos._translate._hmdb_chebi_table',
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
                'omnipath_metabo.datasets.cosmos._translate._metanetx_chebi_table',
                return_value={'MNXM999': 'CHEBI:222'},
            ),
            patch('omnipath_metabo.datasets.cosmos._translate._lipidmaps_chebi_table', return_value={}),
            patch('omnipath_metabo.datasets.cosmos._translate._pubchem_chebi_table', return_value={}),
            patch('omnipath_metabo.datasets.cosmos._translate._kegg_chebi_table', return_value={}),
            patch('omnipath_metabo.datasets.cosmos._translate._hmdb_chebi_table', return_value={}),
        ):
            result = _metatlas_to_chebi('FakeGEM')
        assert result['MAM00010'] == 'CHEBI:222'


# ---------------------------------------------------------------------------
# _build_protein_mapping -- mocked mapping_translate
# ---------------------------------------------------------------------------

class TestBuildProteinMapping:
    """Tests for the batch protein ID translation dispatcher."""

    def test_uniprot_passthrough(self):
        ids = pd.Series(['P00533', 'P04637'])
        result = _build_protein_mapping('uniprot', ids, 9606)
        assert result['P00533'] == frozenset({'P00533'})
        assert result['P04637'] == frozenset({'P04637'})

    def test_reaction_id_passthrough(self):
        ids = pd.Series(['MAR04831', 'ORPHAN_NA_TRANS'])
        result = _build_protein_mapping('reaction_id', ids, 9606)
        assert result['MAR04831'] == frozenset({'MAR04831'})
        assert result['ORPHAN_NA_TRANS'] == frozenset({'ORPHAN_NA_TRANS'})

    def test_ensp_uses_mapping_translate(self):
        ids = pd.Series(['ENSP00000269305', 'ENSP00000350283'])
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_translate',
            return_value={
                'ENSP00000269305': {'P04637'},
                'ENSP00000350283': {'P00533'},
            },
        ):
            result = _build_protein_mapping('ensp', ids, 9606)
        assert result['ENSP00000269305'] == frozenset({'P04637'})
        assert result['ENSP00000350283'] == frozenset({'P00533'})

    def test_ensp_missing_returns_none(self):
        ids = pd.Series(['ENSP_FAKE'])
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_translate',
            return_value={},
        ):
            result = _build_protein_mapping('ensp', ids, 9606)
        assert result['ENSP_FAKE'] is None

    def test_ensembl_simple_ids(self):
        ids = pd.Series(['ENSG00000141510'])
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_translate',
            return_value={'ENSG00000141510': {'P04637'}},
        ):
            result = _build_protein_mapping('ensembl', ids, 9606)
        assert result['ENSG00000141510'] == frozenset({'P04637'})

    def test_ensembl_compound_ids(self):
        """Compound ENSG IDs (A_B) are split and results unioned."""
        ids = pd.Series(['ENSG0001_ENSG0002'])
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_translate',
            return_value={
                'ENSG0001': {'P11111'},
                'ENSG0002': {'P22222'},
            },
        ):
            result = _build_protein_mapping('ensembl', ids, 9606)
        assert result['ENSG0001_ENSG0002'] == frozenset({'P11111', 'P22222'})

    def test_genesymbol_uses_mapping_translate(self):
        ids = pd.Series(['TP53', 'EGFR'])
        with patch(
            'omnipath_metabo.datasets.cosmos._translate.mapping_translate',
            return_value={
                'TP53': {'P04637'},
                'EGFR': {'P00533'},
            },
        ):
            result = _build_protein_mapping('genesymbol', ids, 9606)
        assert result['TP53'] == frozenset({'P04637'})
        assert result['EGFR'] == frozenset({'P00533'})

    def test_entrez_with_bigg_fallback(self):
        """Entrez IDs try BiGG gene symbol -> UniProt first, then batch."""
        ids = pd.Series(['1591', '999999'])
        with (
            patch(
                'omnipath_metabo.datasets.cosmos._translate._entrez_to_uniprot_bigg',
                return_value={'1591': 'P05108'},
            ),
            patch(
                'omnipath_metabo.datasets.cosmos._translate.mapping_translate',
                return_value={'999999': {'Q99999'}},
            ),
        ):
            result = _build_protein_mapping('entrez', ids, 9606)
        assert result['1591'] == frozenset({'P05108'})
        assert result['999999'] == frozenset({'Q99999'})

    def test_unknown_type_returns_none(self):
        ids = pd.Series(['X123'])
        result = _build_protein_mapping('unknown_type', ids, 9606)
        assert result['X123'] is None


# ---------------------------------------------------------------------------
# _build_metab_mapping -- mocked table functions
# ---------------------------------------------------------------------------

class TestBuildMetabMapping:
    """Tests for the batch metabolite ID translation dispatcher."""

    def test_chebi_passthrough(self):
        ids = pd.Series(['CHEBI:30616', 'CHEBI:15422'])
        resource = pd.Series(['STITCH', 'STITCH'])
        result = _build_metab_mapping('chebi', ids, resource)
        assert result['CHEBI:30616'] == 'CHEBI:30616'
        assert result['CHEBI:15422'] == 'CHEBI:15422'

    def test_pubchem_uses_table(self):
        ids = pd.Series(['5793', '2244'])
        resource = pd.Series(['STITCH', 'STITCH'])
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._pubchem_chebi_table',
            return_value={'5793': 'CHEBI:18367', '2244': 'CHEBI:15365'},
        ):
            result = _build_metab_mapping('pubchem', ids, resource)
        assert result['5793'] == 'CHEBI:18367'
        assert result['2244'] == 'CHEBI:15365'

    def test_hmdb_uses_table(self):
        ids = pd.Series(['HMDB0000001'])
        resource = pd.Series(['STITCH'])
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._hmdb_chebi_table',
            return_value={'HMDB0000001': 'CHEBI:16015'},
        ):
            result = _build_metab_mapping('hmdb', ids, resource)
        assert result['HMDB0000001'] == 'CHEBI:16015'

    def test_bigg_uses_table(self):
        ids = pd.Series(['glc__D'])
        resource = pd.Series(['BiGG'])
        with patch(
            'omnipath_metabo.datasets.cosmos._translate._bigg_chebi_table',
            return_value={'glc__D': 'CHEBI:17634'},
        ):
            result = _build_metab_mapping('bigg', ids, resource)
        assert result['glc__D'] == 'CHEBI:17634'

    def test_unknown_type_returns_none(self):
        ids = pd.Series(['X123'])
        resource = pd.Series(['UnknownDB'])
        result = _build_metab_mapping('unknown_type', ids, resource)
        assert result['X123'] is None

    def test_synonym_hmdb_hit(self):
        """Synonym lookup tries HMDB synonyms first."""
        ids = pd.Series(['ATP'])
        resource = pd.Series(['BRENDA'])
        _hmdb_synonyms_chebi.cache_clear()
        _ramp_synonyms_chebi.cache_clear()
        with (
            patch(
                'pypath.inputs.hmdb.metabolites.synonyms_chebi',
                return_value={'atp': 'CHEBI:30616'},
            ),
            patch(
                'pypath.inputs.ramp._mapping.ramp_synonyms_chebi',
                return_value={},
            ),
        ):
            result = _build_metab_mapping('synonym', ids, resource)
        assert result['ATP'] == 'CHEBI:30616'

    def test_synonym_ramp_fallback(self):
        """Synonym lookup falls back to RaMP when HMDB misses."""
        ids = pd.Series(['NADH'])
        resource = pd.Series(['BRENDA'])
        _hmdb_synonyms_chebi.cache_clear()
        _ramp_synonyms_chebi.cache_clear()
        with (
            patch(
                'pypath.inputs.hmdb.metabolites.synonyms_chebi',
                return_value={},
            ),
            patch(
                'pypath.inputs.ramp._mapping.ramp_synonyms_chebi',
                return_value={'nadh': 'CHEBI:57945'},
            ),
        ):
            result = _build_metab_mapping('synonym', ids, resource)
        assert result['NADH'] == 'CHEBI:57945'
