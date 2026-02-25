#!/usr/bin/env python

"""Fast unit tests for omnipath_metabo.datasets.cosmos._format."""

import pandas as pd
import pytest

from omnipath_metabo.datasets.cosmos._format import (
    _assign_n,
    _fmt_gene,
    _fmt_met,
    _is_pre_expanded,
    _other_comp,
    _row_category,
    format_pkn,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_INTERACTION_COLS = [
    'source', 'target', 'source_type', 'target_type',
    'id_type_a', 'id_type_b', 'interaction_type', 'resource',
    'mor', 'locations', 'attrs',
]


def _row(
    source='CHEBI:15422',
    target='ENSG00000141510',
    source_type='small_molecule',
    target_type='protein',
    id_type_a='chebi',
    id_type_b='ensg',
    interaction_type='transport',
    resource='TCDB',
    mor=1,
    locations=('e', 'c'),
    attrs=None,
):
    """Build a minimal post-translation PKN row dict."""
    return {
        'source': source,
        'target': target,
        'source_type': source_type,
        'target_type': target_type,
        'id_type_a': id_type_a,
        'id_type_b': id_type_b,
        'interaction_type': interaction_type,
        'resource': resource,
        'mor': mor,
        'locations': locations,
        'attrs': attrs if attrs is not None else {},
    }


def _df(*rows):
    """Build a PKN DataFrame from one or more row dicts."""
    return pd.DataFrame(list(rows), columns=_INTERACTION_COLS)


# ---------------------------------------------------------------------------
# Node ID formatting helpers
# ---------------------------------------------------------------------------

class TestFmtMet:
    def test_with_compartment(self):
        assert _fmt_met('CHEBI:15422', 'c') == 'Metab__CHEBI:15422_c'

    def test_other_compartment(self):
        assert _fmt_met('CHEBI:9251', 'l') == 'Metab__CHEBI:9251_l'

    def test_no_compartment(self):
        assert _fmt_met('CHEBI:15422', '') == 'Metab__CHEBI:15422'


class TestFmtGene:
    def test_forward(self):
        assert _fmt_gene('ENSG00000141510', 1) == 'Gene1__ENSG00000141510'

    def test_forward_large_n(self):
        assert _fmt_gene('ENSG00000141510', 42) == 'Gene42__ENSG00000141510'

    def test_reverse(self):
        assert _fmt_gene('ENSG00000141510', 1, reverse=True) == 'Gene1__ENSG00000141510_rev'

    def test_forward_explicit_false(self):
        assert _fmt_gene('ENSG00000141510', 5, reverse=False) == 'Gene5__ENSG00000141510'


class TestOtherComp:
    def test_plasma_membrane(self):
        assert _other_comp(('e', 'c')) == 'e'

    def test_mitochondrial(self):
        assert _other_comp(('m', 'c')) == 'm'

    def test_first_non_c(self):
        assert _other_comp(('e', 'm', 'c')) == 'e'

    def test_only_c(self):
        assert _other_comp(('c',)) == ''

    def test_empty(self):
        assert _other_comp(()) == ''


# ---------------------------------------------------------------------------
# Resource classification
# ---------------------------------------------------------------------------

class TestIsPreExpanded:
    def test_gem_metabolic(self):
        assert _is_pre_expanded('GEM:Human-GEM')

    def test_gem_transporter(self):
        assert _is_pre_expanded('GEM_transporter:Human-GEM')

    def test_recon3d(self):
        assert _is_pre_expanded('Recon3D')

    def test_tcdb(self):
        assert not _is_pre_expanded('TCDB')

    def test_slc(self):
        assert not _is_pre_expanded('SLC')

    def test_stitch(self):
        assert not _is_pre_expanded('STITCH')

    def test_brenda(self):
        assert not _is_pre_expanded('BRENDA')


class TestRowCategory:
    def test_tcdb_transport(self):
        assert _row_category('transport', 'TCDB') == 'transporter'

    def test_slc_transport(self):
        assert _row_category('transport', 'SLC') == 'transporter'

    def test_gem_transporter(self):
        assert _row_category('transport', 'GEM_transporter:Human-GEM') == 'transporter'

    def test_stitch_transporter(self):
        assert _row_category('transporter', 'STITCH') == 'transporter'

    def test_mrclinksdb_ligand_receptor(self):
        assert _row_category('ligand_receptor', 'MRCLinksDB') == 'receptor'

    def test_stitch_receptor(self):
        assert _row_category('receptor', 'STITCH') == 'receptor'

    def test_brenda_allosteric(self):
        assert _row_category('allosteric_regulation', 'BRENDA') == 'other'

    def test_gem_metabolic(self):
        assert _row_category('transport', 'GEM:Human-GEM') == 'transporter'

    def test_stitch_other(self):
        assert _row_category('other', 'STITCH') == 'other'


# ---------------------------------------------------------------------------
# N assignment
# ---------------------------------------------------------------------------

class TestAssignN:
    def test_non_gem_each_row_gets_own_n(self):
        df = _df(
            _row(interaction_type='transport', resource='TCDB'),
            _row(interaction_type='transport', resource='TCDB'),
        )
        df['_category'] = 'transporter'
        ns = _assign_n(df)
        assert list(ns) == [1, 2]

    def test_gem_rows_share_n_for_same_reaction(self):
        """All 4 rows of a GEM reaction share the same N."""
        rows = [
            _row(
                resource='GEM:Human-GEM',
                interaction_type='transport',
                attrs={'reaction_id': 'MAR001', 'reverse': False},
            )
            for _ in range(4)
        ]
        df = _df(*rows)
        df['_category'] = 'transporter'
        ns = _assign_n(df)
        assert list(ns) == [1, 1, 1, 1]

    def test_gem_different_reactions_get_different_n(self):
        rows = [
            _row(resource='GEM:Human-GEM', attrs={'reaction_id': 'MAR001', 'reverse': False}),
            _row(resource='GEM:Human-GEM', attrs={'reaction_id': 'MAR002', 'reverse': False}),
        ]
        df = _df(*rows)
        df['_category'] = 'other'
        ns = _assign_n(df)
        assert list(ns) == [1, 2]

    def test_n_resets_between_categories(self):
        rows = [
            _row(interaction_type='transport', resource='TCDB'),
            _row(interaction_type='ligand_receptor', resource='MRCLinksDB'),
        ]
        df = _df(*rows)
        df['_category'] = df.apply(
            lambda r: _row_category(r['interaction_type'], r['resource']), axis=1
        )
        ns = _assign_n(df)
        # transporter gets N=1, receptor also starts at N=1
        assert list(ns) == [1, 1]


# ---------------------------------------------------------------------------
# Transporter expansion
# ---------------------------------------------------------------------------

class TestFormatTransporter:
    """Non-pre-expanded transporter: 1 input row → 4 output rows."""

    def _run(self, locations=('e', 'c')):
        df = _df(_row(
            source='CHEBI:15422',
            target='ENSG00000141510',
            interaction_type='transport',
            resource='TCDB',
            locations=locations,
        ))
        return format_pkn(df, include_orphans=True)

    def test_four_main_rows_generated(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 4

    def test_row1_met_other_to_gene_fwd(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        r = main.iloc[0]
        assert r['source'] == 'Metab__CHEBI:15422_e'
        assert r['target'] == 'Gene1__ENSG00000141510'
        assert r['attrs']['reverse'] is False

    def test_row2_gene_fwd_to_met_c(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        r = main.iloc[1]
        assert r['source'] == 'Gene1__ENSG00000141510'
        assert r['target'] == 'Metab__CHEBI:15422_c'
        assert r['attrs']['reverse'] is False

    def test_row3_met_c_to_gene_rev(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        r = main.iloc[2]
        assert r['source'] == 'Metab__CHEBI:15422_c'
        assert r['target'] == 'Gene1__ENSG00000141510_rev'
        assert r['attrs']['reverse'] is True

    def test_row4_gene_rev_to_met_other(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        r = main.iloc[3]
        assert r['source'] == 'Gene1__ENSG00000141510_rev'
        assert r['target'] == 'Metab__CHEBI:15422_e'
        assert r['attrs']['reverse'] is True

    def test_stitch_no_compartment(self):
        """STITCH has no compartment; both metabolite nodes are the same bare ID."""
        result = self._run(locations=())
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 4
        sources = set(main['source'])
        targets = set(main['target'])
        # Both metabolite nodes should be the same bare Metab__ (no suffix)
        assert 'Metab__CHEBI:15422' in sources | targets
        # Gene nodes should differ only by _rev
        assert 'Gene1__ENSG00000141510' in sources | targets
        assert 'Gene1__ENSG00000141510_rev' in sources | targets

    def test_cosmos_formatted_flag_set(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        assert all(r['cosmos_formatted'] is True for r in main['attrs'])

    def test_source_target_type_correct(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        # Rows 0 and 2 are met→gene
        assert main.iloc[0]['source_type'] == 'small_molecule'
        assert main.iloc[0]['target_type'] == 'protein'
        # Rows 1 and 3 are gene→met
        assert main.iloc[1]['source_type'] == 'protein'
        assert main.iloc[1]['target_type'] == 'small_molecule'


# ---------------------------------------------------------------------------
# Connector edges
# ---------------------------------------------------------------------------

class TestConnectorEdges:
    def test_connector_edges_present(self):
        df = _df(_row(
            source='CHEBI:15422',
            target='ENSG00000141510',
            interaction_type='transport',
            resource='TCDB',
            locations=('e', 'c'),
        ))
        result = format_pkn(df)
        conn = result[result['interaction_type'] == 'connector']
        assert len(conn) > 0

    def test_gene_connector_present(self):
        df = _df(_row(
            source='CHEBI:15422',
            target='ENSG00000141510',
            interaction_type='transport',
            resource='TCDB',
            locations=('e', 'c'),
        ))
        result = format_pkn(df)
        conn = result[result['interaction_type'] == 'connector']
        sources = set(conn['source'])
        targets = set(conn['target'])
        # Bare ENSG → Gene1__ENSG (forward connector)
        assert 'ENSG00000141510' in sources
        assert 'Gene1__ENSG00000141510' in targets
        # Reverse connector
        assert 'Gene1__ENSG00000141510_rev' in targets

    def test_metabolite_connector_present(self):
        df = _df(_row(
            source='CHEBI:15422',
            target='ENSG00000141510',
            interaction_type='transport',
            resource='TCDB',
            locations=('e', 'c'),
        ))
        result = format_pkn(df)
        conn = result[result['interaction_type'] == 'connector']
        targets = set(conn['target'])
        assert 'Metab__CHEBI:15422_e' in targets
        assert 'Metab__CHEBI:15422_c' in targets

    def test_connectors_deduplicated(self):
        """Two rows with the same gene/metabolite → connectors not duplicated."""
        df = _df(
            _row(source='CHEBI:15422', target='ENSG00000141510',
                 interaction_type='transport', resource='TCDB',
                 locations=('e', 'c')),
            _row(source='CHEBI:15422', target='ENSG00000141510',
                 interaction_type='transport', resource='TCDB',
                 locations=('e', 'c')),
        )
        result = format_pkn(df)
        conn = result[result['interaction_type'] == 'connector']
        # No duplicate (source, target) pairs in connectors
        dupes = conn.duplicated(subset=['source', 'target'])
        assert not dupes.any()

    def test_connector_resource_label(self):
        df = _df(_row(interaction_type='transport', resource='TCDB',
                      locations=('e', 'c')))
        result = format_pkn(df)
        conn = result[result['interaction_type'] == 'connector']
        assert (conn['resource'] == 'COSMOS_formatter').all()


# ---------------------------------------------------------------------------
# Pre-expanded rows (GEM / Recon3D)
# ---------------------------------------------------------------------------

class TestFormatPreExpanded:
    def _gem_rows(self, resource='GEM:Human-GEM'):
        """Four synthetic GEM rows: met→gene and gene→met, forward and reverse."""
        return _df(
            _row(source='CHEBI:15422', target='ENSG00000141510',
                 source_type='small_molecule', target_type='protein',
                 interaction_type='transport', resource=resource,
                 locations=('c',),
                 attrs={'reaction_id': 'MAR001', 'reverse': False}),
            _row(source='ENSG00000141510', target='CHEBI:30031',
                 source_type='protein', target_type='small_molecule',
                 interaction_type='transport', resource=resource,
                 locations=('e',),
                 attrs={'reaction_id': 'MAR001', 'reverse': False}),
            _row(source='CHEBI:30031', target='ENSG00000141510',
                 source_type='small_molecule', target_type='protein',
                 interaction_type='transport', resource=resource,
                 locations=('e',),
                 attrs={'reaction_id': 'MAR001', 'reverse': True}),
            _row(source='ENSG00000141510', target='CHEBI:15422',
                 source_type='protein', target_type='small_molecule',
                 interaction_type='transport', resource=resource,
                 locations=('c',),
                 attrs={'reaction_id': 'MAR001', 'reverse': True}),
        )

    def test_gem_no_extra_rows_generated(self):
        """Pre-expanded rows should not be expanded further."""
        df = self._gem_rows()
        result = format_pkn(df)
        main = result[result['interaction_type'] != 'connector']
        # 4 in → 4 out (no new rows)
        assert len(main) == 4

    def test_gem_all_rows_share_same_n(self):
        df = self._gem_rows()
        result = format_pkn(df)
        main = result[result['interaction_type'] != 'connector']
        gene_nodes = (
            main['source'].str.startswith('Gene').fillna(False) |
            main['target'].str.startswith('Gene').fillna(False)
        )
        gene_cols = pd.concat([
            main.loc[main['source'].str.startswith('Gene', na=False), 'source'],
            main.loc[main['target'].str.startswith('Gene', na=False), 'target'],
        ])
        # All gene nodes should start with Gene1__
        assert all(n.startswith('Gene1__') for n in gene_cols)

    def test_gem_reverse_row_gets_rev_suffix(self):
        df = self._gem_rows()
        result = format_pkn(df)
        main = result[result['interaction_type'] != 'connector']
        rev_rows = main[main['attrs'].apply(lambda a: a.get('reverse', False))]
        for _, r in rev_rows.iterrows():
            gene_node = r['source'] if r['source_type'] == 'protein' else r['target']
            assert gene_node.endswith('_rev'), f'{gene_node!r} should end with _rev'

    def test_gem_forward_row_no_rev_suffix(self):
        df = self._gem_rows()
        result = format_pkn(df)
        main = result[result['interaction_type'] != 'connector']
        fwd_rows = main[~main['attrs'].apply(lambda a: a.get('reverse', False))]
        for _, r in fwd_rows.iterrows():
            gene_node = r['source'] if r['source_type'] == 'protein' else r['target']
            assert not gene_node.endswith('_rev'), f'{gene_node!r} should not end with _rev'

    def test_recon3d_treated_as_pre_expanded(self):
        df = self._gem_rows(resource='Recon3D')
        result = format_pkn(df)
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 4


# ---------------------------------------------------------------------------
# Receptor and other (simple) rows
# ---------------------------------------------------------------------------

class TestFormatSimple:
    def test_receptor_no_reverse_gene(self):
        df = _df(_row(
            source='CHEBI:9251',
            target='ENSG00000123456',
            interaction_type='ligand_receptor',
            resource='MRCLinksDB',
            locations=(),
        ))
        result = format_pkn(df)
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 1
        assert main.iloc[0]['source'] == 'Metab__CHEBI:9251'
        assert main.iloc[0]['target'] == 'Gene1__ENSG00000123456'

    def test_other_brenda_no_reverse(self):
        df = _df(_row(
            source='CHEBI:30031',
            target='ENSG00000099876',
            interaction_type='allosteric_regulation',
            resource='BRENDA',
            locations=('c',),
        ))
        result = format_pkn(df)
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 1
        assert main.iloc[0]['source'] == 'Metab__CHEBI:30031_c'
        assert main.iloc[0]['target'] == 'Gene1__ENSG00000099876'

    def test_receptor_no_rev_connector(self):
        """Receptors should not produce a *_rev gene connector."""
        df = _df(_row(
            source='CHEBI:9251',
            target='ENSG00000123456',
            interaction_type='ligand_receptor',
            resource='MRCLinksDB',
            locations=(),
        ))
        result = format_pkn(df)
        conn = result[result['interaction_type'] == 'connector']
        rev_conn = conn[conn['target'].str.endswith('_rev', na=False)]
        assert len(rev_conn) == 0


# ---------------------------------------------------------------------------
# N counter across categories
# ---------------------------------------------------------------------------

class TestNCounterCategories:
    def test_n_independent_per_category(self):
        """Transporter and receptor counters are independent."""
        df = _df(
            _row(interaction_type='transport', resource='TCDB',
                 locations=('e', 'c')),
            _row(interaction_type='transport', resource='TCDB',
                 locations=('m', 'c')),
            _row(interaction_type='ligand_receptor', resource='MRCLinksDB',
                 locations=()),
        )
        result = format_pkn(df)
        main = result[result['interaction_type'] != 'connector']
        gene_nodes = pd.concat([
            main.loc[main['source'].str.startswith('Gene', na=False), 'source'],
            main.loc[main['target'].str.startswith('Gene', na=False), 'target'],
        ]).unique()

        # Transporters: Gene1__... and Gene2__...
        # Receptor: also Gene1__... (counter resets)
        prefixes = {n.split('__')[0] for n in gene_nodes}
        assert 'Gene1' in prefixes
        assert 'Gene2' in prefixes


# ---------------------------------------------------------------------------
# Orphan handling
# ---------------------------------------------------------------------------

class TestOrphanHandling:
    def _orphan_df(self):
        return _df(
            _row(
                source='CHEBI:15422',
                target='MAR99999',
                target_type='protein',
                id_type_b='reaction_id',
                interaction_type='transport',
                resource='Recon3D',
                locations=('c',),
                attrs={'reaction_id': 'MAR99999', 'orphan': True, 'reverse': False},
            )
        )

    def test_orphan_kept_by_default(self):
        result = format_pkn(self._orphan_df())
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 1

    def test_orphan_dropped_when_excluded(self):
        result = format_pkn(self._orphan_df(), include_orphans=False)
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 0


# ---------------------------------------------------------------------------
# cosmos_formatted flag
# ---------------------------------------------------------------------------

class TestCosmosFormattedFlag:
    def test_flag_set_on_all_rows(self):
        df = _df(
            _row(interaction_type='transport', resource='TCDB', locations=('e', 'c')),
            _row(interaction_type='ligand_receptor', resource='MRCLinksDB', locations=()),
        )
        result = format_pkn(df)
        for _, r in result.iterrows():
            assert r['attrs']['cosmos_formatted'] is True

    def test_empty_input(self):
        df = pd.DataFrame(columns=_INTERACTION_COLS)
        result = format_pkn(df)
        assert result.empty
