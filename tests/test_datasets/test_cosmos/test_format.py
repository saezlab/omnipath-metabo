#!/usr/bin/env python

"""Fast unit tests for omnipath_metabo.datasets.cosmos._format."""

import pandas as pd
import pytest

from omnipath_metabo.datasets.cosmos._format import (
    _assign_n,
    _fmt_gene,
    _fmt_met,
    _fmt_rxn,
    _is_pre_expanded,
    _row_category,
    format_pkn,
    format_allosteric,
    format_enzyme_metabolite,
    format_receptors,
    format_transporters,
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


def _call_format(df, **kw) -> pd.DataFrame:
    """Call format_pkn and return all edges (incl. connectors) as a DataFrame."""
    from omnipath_metabo.datasets.cosmos._record import CosmosEdge
    bundle = format_pkn(df, **kw)
    if not bundle.network:
        return pd.DataFrame(columns=list(CosmosEdge._fields))
    return pd.DataFrame(bundle.network)


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

    def test_single_ac(self):
        assert _fmt_gene('P12345', 1) == 'Gene1__P12345'

    def test_single_ac_reverse(self):
        assert _fmt_gene('P12345', 2, reverse=True) == 'Gene2__P12345_rev'


class TestFmtRxn:
    def test_forward(self):
        assert _fmt_rxn('MAR99999', 1) == 'Gene1__orphanReacMAR99999'

    def test_forward_large_n(self):
        assert _fmt_rxn('MAR99999', 42) == 'Gene42__orphanReacMAR99999'

    def test_reverse(self):
        assert _fmt_rxn('MAR99999', 1, reverse=True) == 'Gene1__orphanReacMAR99999_rev'

    def test_forward_explicit_false(self):
        assert _fmt_rxn('MAR99999', 5, reverse=False) == 'Gene5__orphanReacMAR99999'

    def test_recon3d_reaction_id(self):
        assert _fmt_rxn('R_TKT', 3) == 'Gene3__orphanReacR_TKT'



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
        return _call_format(df, include_orphans=True)

    def test_eight_main_rows_generated(self):
        """locations=('e','c') → 2 comps × 4 rows each = 8 main rows."""
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 8

    def test_four_main_rows_single_location(self):
        """Single location → exactly 4 main rows."""
        result = self._run(locations=('e',))
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 4

    def test_row1_met_e_to_gene_fwd(self):
        """First group (src_comp='e'): met_e → Gene_fwd."""
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

    def test_row4_gene_rev_to_met_e(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        r = main.iloc[3]
        assert r['source'] == 'Gene1__ENSG00000141510_rev'
        assert r['target'] == 'Metab__CHEBI:15422_e'
        assert r['attrs']['reverse'] is True

    def test_second_group_uses_c_as_src(self):
        """Second group (src_comp='c'): met_c → Gene_fwd."""
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        r = main.iloc[4]
        assert r['source'] == 'Metab__CHEBI:15422_c'
        assert r['target'] == 'Gene1__ENSG00000141510'

    def test_cosmos_formatted_flag_set(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        assert all(r['cosmos_formatted'] is True for r in main['attrs'])

    def test_source_target_type_correct(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        assert main.iloc[0]['source_type'] == 'small_molecule'
        assert main.iloc[0]['target_type'] == 'protein'
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
        result = _call_format(df)
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
        result = _call_format(df)
        conn = result[result['interaction_type'] == 'connector']
        sources = set(conn['source'])
        targets = set(conn['target'])
        assert 'ENSG00000141510' in sources
        assert 'Gene1__ENSG00000141510' in targets
        assert 'Gene1__ENSG00000141510_rev' in targets

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
        result = _call_format(df)
        conn = result[result['interaction_type'] == 'connector']
        # No duplicate (source, target) pairs in connectors
        dupes = conn.duplicated(subset=['source', 'target'])
        assert not dupes.any()

    def test_connector_resource_label(self):
        df = _df(_row(interaction_type='transport', resource='TCDB',
                      locations=('e', 'c')))
        result = _call_format(df)
        conn = result[result['interaction_type'] == 'connector']
        assert (conn['resource'] == 'COSMOS_formatter').all()

    def test_frozenset_gene_expands_to_separate_rows(self):
        """frozenset target → one set of rows per AC, not joined with ';'."""
        df = _df(_row(
            target=frozenset({'Q99999', 'P12345'}),
            interaction_type='transport',
            resource='TCDB',
            locations=('e',),
        ))
        result = _call_format(df)
        main = result[result['interaction_type'] != 'connector']
        # 1 location × 2 ACs × 4 rows = 8 main rows
        assert len(main) == 8
        gene_nodes = set(main[main['target_type'] == 'protein']['target'])
        assert 'Gene1__P12345' in gene_nodes
        assert 'Gene1__Q99999' in gene_nodes
        assert not any(';' in n for n in gene_nodes)

    def test_frozenset_gene_separate_connectors(self):
        """frozenset target → one connector per AC (not a joined bare ID)."""
        df = _df(_row(
            target=frozenset({'Q99999', 'P12345'}),
            interaction_type='transport',
            resource='TCDB',
            locations=('e',),
        ))
        result = _call_format(df)
        conn = result[result['interaction_type'] == 'connector']
        sources = set(conn['source'])
        assert 'P12345' in sources
        assert 'Q99999' in sources
        assert 'P12345;Q99999' not in sources


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
        result = _call_format(df)
        main = result[result['interaction_type'] != 'connector']
        # 4 in → 4 out (no new rows)
        assert len(main) == 4

    def test_gem_all_rows_share_same_n(self):
        df = self._gem_rows()
        result = _call_format(df)
        main = result[result['interaction_type'] != 'connector']
        gene_cols = pd.concat([
            main.loc[main['source'].str.startswith('Gene', na=False), 'source'],
            main.loc[main['target'].str.startswith('Gene', na=False), 'target'],
        ])
        # All gene nodes should start with Gene1__
        assert all(n.startswith('Gene1__') for n in gene_cols)

    def test_gem_reverse_row_gets_rev_suffix(self):
        df = self._gem_rows()
        result = _call_format(df)
        main = result[result['interaction_type'] != 'connector']
        rev_rows = main[main['attrs'].apply(lambda a: a.get('reverse', False))]
        for _, r in rev_rows.iterrows():
            gene_node = r['source'] if r['source_type'] == 'protein' else r['target']
            assert gene_node.endswith('_rev'), f'{gene_node!r} should end with _rev'

    def test_gem_forward_row_no_rev_suffix(self):
        df = self._gem_rows()
        result = _call_format(df)
        main = result[result['interaction_type'] != 'connector']
        fwd_rows = main[~main['attrs'].apply(lambda a: a.get('reverse', False))]
        for _, r in fwd_rows.iterrows():
            gene_node = r['source'] if r['source_type'] == 'protein' else r['target']
            assert not gene_node.endswith('_rev'), f'{gene_node!r} should not end with _rev'

    def test_recon3d_treated_as_pre_expanded(self):
        df = self._gem_rows(resource='Recon3D')
        result = _call_format(df)
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 4


# ---------------------------------------------------------------------------
# Receptor and other (simple) rows
# ---------------------------------------------------------------------------

class TestFormatSimple:
    def test_receptor_bare_uniprot_target(self):
        df = _df(_row(
            source='CHEBI:9251',
            target='ENSG00000123456',
            interaction_type='ligand_receptor',
            resource='MRCLinksDB',
            locations=(),
        ))
        result = _call_format(df)
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 1
        assert main.iloc[0]['source'] == 'Metab__CHEBI:9251'
        assert main.iloc[0]['target'] == 'ENSG00000123456'

    def test_other_brenda_bare_uniprot(self):
        df = _df(_row(
            source='CHEBI:30031',
            target='ENSG00000099876',
            interaction_type='allosteric_regulation',
            resource='BRENDA',
            locations=('c',),
        ))
        result = _call_format(df)
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 1
        assert main.iloc[0]['source'] == 'Metab__CHEBI:30031_c'
        assert main.iloc[0]['target'] == 'ENSG00000099876'

    def test_receptor_no_rev_connector(self):
        """Receptors should not produce a *_rev gene connector."""
        df = _df(_row(
            source='CHEBI:9251',
            target='ENSG00000123456',
            interaction_type='ligand_receptor',
            resource='MRCLinksDB',
            locations=(),
        ))
        result = _call_format(df)
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
        result = _call_format(df)
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
    def _orphan_row(self, reverse=False, resource='Recon3D'):
        """One pre-expanded orphan transport row (met → reaction_id)."""
        return _row(
            source='CHEBI:15422',
            target='MAR99999',
            target_type='protein',
            id_type_b='reaction_id',
            interaction_type='transport',
            resource=resource,
            locations=('c',),
            attrs={'reaction_id': 'MAR99999', 'orphan': True, 'reverse': reverse},
        )

    def _orphan_df(self, **kw):
        return _df(self._orphan_row(**kw))

    def test_orphan_kept_by_default(self):
        result = _call_format(self._orphan_df())
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 1

    def test_orphan_dropped_when_excluded(self):
        result = _call_format(self._orphan_df(), include_orphans=False)
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 0

    def test_orphan_uses_gene_orphanreac_prefix(self):
        """Orphan node must use Gene{N}__orphanReac prefix."""
        result = _call_format(self._orphan_df())
        main = result[result['interaction_type'] != 'connector']
        target = main.iloc[0]['target']
        assert 'orphanReac' in target, f'Expected orphanReac in node ID, got: {target!r}'
        assert target.startswith('Gene'), f'Expected Gene prefix, got: {target!r}'

    def test_orphan_forward_node_id(self):
        result = _call_format(self._orphan_df())
        main = result[result['interaction_type'] != 'connector']
        assert main.iloc[0]['target'] == 'Gene1__orphanReacMAR99999'

    def test_orphan_reverse_node_has_rev_suffix(self):
        result = _call_format(self._orphan_df(reverse=True))
        main = result[result['interaction_type'] != 'connector']
        target = main.iloc[0]['target']
        assert target == 'Gene1__orphanReacMAR99999_rev'

    def test_orphan_connector_is_rxn_id_to_orphan_node(self):
        """Connector: bare reaction_id → Gene{N}__orphanReac<id> (plain string, no frozenset)."""
        result = _call_format(self._orphan_df())
        conn = result[result['interaction_type'] == 'connector']
        rxn_conn = conn[conn['source'] == 'MAR99999']
        assert len(rxn_conn) == 1
        assert rxn_conn.iloc[0]['target'] == 'Gene1__orphanReacMAR99999'

    def test_orphan_no_extra_plain_gene_connector(self):
        """No Gene__ connector without orphanReac should appear for an orphan row."""
        result = _call_format(self._orphan_df())
        conn = result[result['interaction_type'] == 'connector']
        plain_gene_conn = conn[
            conn['target'].str.startswith('Gene', na=False) &
            ~conn['target'].str.contains('orphanReac', na=False)
        ]
        assert len(plain_gene_conn) == 0

    def test_orphan_gem_transporter_uses_gene_orphanreac_prefix(self):
        """GEM_transporter orphan rows also use Gene{N}__orphanReac prefix."""
        result = _call_format(self._orphan_df(resource='GEM_transporter:Human-GEM'))
        main = result[result['interaction_type'] != 'connector']
        target = main.iloc[0]['target']
        assert target.startswith('Gene') and 'orphanReac' in target

    def test_non_orphan_gem_still_uses_gene_prefix(self):
        """Regression: non-orphan GEM rows must still produce Gene{N}__ nodes."""
        df = _df(_row(
            source='CHEBI:15422',
            target='ENSG00000141510',
            source_type='small_molecule',
            target_type='protein',
            interaction_type='transport',
            resource='GEM_transporter:Human-GEM',
            locations=('c',),
            attrs={'reaction_id': 'MAR001', 'reverse': False},
        ))
        result = _call_format(df)
        main = result[result['interaction_type'] != 'connector']
        target = main.iloc[0]['target']
        assert target.startswith('Gene'), f'Non-orphan must use Gene prefix: {target!r}'


# ---------------------------------------------------------------------------
# cosmos_formatted flag
# ---------------------------------------------------------------------------

class TestCosmosFormattedFlag:
    def test_flag_set_on_all_rows(self):
        df = _df(
            _row(interaction_type='transport', resource='TCDB', locations=('e', 'c')),
            _row(interaction_type='ligand_receptor', resource='MRCLinksDB', locations=()),
        )
        result = _call_format(df)
        for _, r in result.iterrows():
            assert r['attrs']['cosmos_formatted'] is True

    def test_empty_input(self):
        df = pd.DataFrame(columns=_INTERACTION_COLS)
        result = _call_format(df)
        assert result.empty


# ---------------------------------------------------------------------------
# Category-specific format wrappers
# ---------------------------------------------------------------------------

def _make_bundle(*rows):
    """Build a CosmosBundle from row dicts (post-translation Interaction records)."""
    from omnipath_metabo.datasets.cosmos._bundle import CosmosBundle
    from omnipath_metabo.datasets.cosmos._record import Interaction
    return CosmosBundle(
        network=[Interaction(**r) for r in rows],
    )


def _bundle_main_nodes(bundle) -> pd.DataFrame:
    """Return non-connector edges from a formatted bundle as a DataFrame."""
    from omnipath_metabo.datasets.cosmos._record import CosmosEdge
    df = pd.DataFrame(bundle.network)
    return df[df['interaction_type'] != 'connector']


class TestFormatTransporters:
    def _bundle(self):
        return _make_bundle(
            _row(interaction_type='transport', resource='TCDB', locations=('e', 'c')),
            _row(interaction_type='ligand_receptor', resource='MRCLinksDB', locations=()),
            _row(interaction_type='allosteric_regulation', resource='BRENDA', locations=('c',)),
        )

    def test_only_transporter_rows_in_output(self):
        main = _bundle_main_nodes(format_transporters(self._bundle()))
        resources = set(main['resource'])
        assert 'TCDB' in resources
        assert 'MRCLinksDB' not in resources
        assert 'BRENDA' not in resources

    def test_returns_cosmos_bundle(self):
        from omnipath_metabo.datasets.cosmos._bundle import CosmosBundle
        assert isinstance(format_transporters(self._bundle()), CosmosBundle)

    def test_plain_df_passed_through(self):
        """Plain DataFrame (not a bundle) bypasses filtering."""
        df = _df(_row(interaction_type='transport', resource='TCDB', locations=('e', 'c')))
        result = format_transporters(df)
        main = _bundle_main_nodes(result)
        # 2 locations × 4 rows = 8 main rows
        assert len(main) == 8

    def test_category_pure_bundle_is_noop(self):
        """Bundle already containing only transporters: same output as format_pkn."""
        bundle = _make_bundle(
            _row(interaction_type='transport', resource='TCDB', locations=('e', 'c')),
        )
        assert len(_bundle_main_nodes(format_transporters(bundle))) == 8


class TestFormatReceptorRow:
    """Receptor rows: bare UniProt target, one row per compartment."""

    def _run(self, locations=('e', 'c')):
        df = _df(_row(
            source='CHEBI:15422',
            target='P00533',
            source_type='small_molecule',
            target_type='protein',
            interaction_type='ligand_receptor',
            resource='MRCLinksDB',
            locations=locations,
        ))
        return _call_format(df)

    def test_two_rows_for_two_locations(self):
        """locations=('e','c') → 2 separate rows, one per compartment."""
        result = self._run(locations=('e', 'c'))
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 2

    def test_no_location_gives_one_row(self):
        result = self._run(locations=())
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 1

    def test_single_location_gives_one_row(self):
        result = self._run(locations=('e',))
        main = result[result['interaction_type'] != 'connector']
        assert len(main) == 1

    def test_protein_node_is_bare_uniprot(self):
        result = self._run(locations=('e', 'c'))
        main = result[result['interaction_type'] != 'connector']
        assert all(r == 'P00533' for r in main['target'])
        assert not any(r.startswith('Gene') for r in main['target'])

    def test_metabolite_nodes_separate_per_compartment(self):
        """Each compartment gets its own metabolite node ID."""
        result = self._run(locations=('e', 'c'))
        main = result[result['interaction_type'] != 'connector']
        sources = set(main['source'])
        assert 'Metab__CHEBI:15422_e' in sources
        assert 'Metab__CHEBI:15422_c' in sources
        assert 'Metab__CHEBI:15422_e;c' not in sources

    def test_single_location(self):
        result = self._run(locations=('e',))
        main = result[result['interaction_type'] != 'connector']
        assert main.iloc[0]['source'] == 'Metab__CHEBI:15422_e'

    def test_no_location_metabolite_has_no_suffix(self):
        result = self._run(locations=())
        main = result[result['interaction_type'] != 'connector']
        assert main.iloc[0]['source'] == 'Metab__CHEBI:15422'

    def test_cosmos_formatted_flag(self):
        result = self._run()
        main = result[result['interaction_type'] != 'connector']
        assert all(r['cosmos_formatted'] is True for r in main['attrs'])

    def test_no_gene_connector(self):
        """No Gene__ connector — protein node is already the bare UniProt AC."""
        result = self._run(locations=('e', 'c'))
        conn = result[result['interaction_type'] == 'connector']
        gene_conn = conn[conn['target'].str.startswith('Gene', na=False)]
        assert len(gene_conn) == 0


class TestFormatReceptors:
    def _bundle(self):
        return _make_bundle(
            _row(interaction_type='ligand_receptor', resource='MRCLinksDB', locations=()),
            _row(interaction_type='transport', resource='TCDB', locations=('e', 'c')),
            _row(interaction_type='allosteric_regulation', resource='BRENDA', locations=('c',)),
        )

    def test_only_receptor_rows_in_output(self):
        main = _bundle_main_nodes(format_receptors(self._bundle()))
        resources = set(main['resource'])
        assert 'MRCLinksDB' in resources
        assert 'TCDB' not in resources
        assert 'BRENDA' not in resources

    def test_returns_cosmos_bundle(self):
        from omnipath_metabo.datasets.cosmos._bundle import CosmosBundle
        assert isinstance(format_receptors(self._bundle()), CosmosBundle)


class TestFormatAllosteric:
    def _bundle(self):
        return _make_bundle(
            _row(interaction_type='allosteric_regulation', resource='BRENDA', locations=('c',)),
            _row(interaction_type='other', resource='STITCH', locations=()),
            _row(interaction_type='transport', resource='TCDB', locations=('e', 'c')),
            _row(interaction_type='ligand_receptor', resource='MRCLinksDB', locations=()),
        )

    def test_brenda_allosteric_included(self):
        main = _bundle_main_nodes(format_allosteric(self._bundle()))
        assert 'BRENDA' in set(main['resource'])

    def test_stitch_other_included(self):
        main = _bundle_main_nodes(format_allosteric(self._bundle()))
        assert 'STITCH' in set(main['resource'])

    def test_transporter_excluded(self):
        main = _bundle_main_nodes(format_allosteric(self._bundle()))
        assert 'TCDB' not in set(main['resource'])

    def test_receptor_excluded(self):
        main = _bundle_main_nodes(format_allosteric(self._bundle()))
        assert 'MRCLinksDB' not in set(main['resource'])

    def test_returns_cosmos_bundle(self):
        from omnipath_metabo.datasets.cosmos._bundle import CosmosBundle
        assert isinstance(format_allosteric(self._bundle()), CosmosBundle)


class TestFormatEnzymeMetabolite:
    def _bundle(self):
        return _make_bundle(
            _row(
                source='CHEBI:15422', target='ENSG00000141510',
                source_type='small_molecule', target_type='protein',
                id_type_a='chebi', id_type_b='ensg',
                interaction_type='catalysis', resource='GEM:Human-GEM',
                locations=('c',),
                attrs={'reaction_id': 'MAR001', 'reverse': False},
            ),
            _row(
                source='CHEBI:15422', target='ENSG00000141510',
                source_type='small_molecule', target_type='protein',
                id_type_a='chebi', id_type_b='ensg',
                interaction_type='catalysis', resource='GEM_transporter:Human-GEM',
                locations=('e',),
                attrs={'reaction_id': 'MAR002', 'reverse': False},
            ),
            _row(interaction_type='transport', resource='TCDB', locations=('e', 'c')),
        )

    def test_gem_metabolic_included(self):
        main = _bundle_main_nodes(format_enzyme_metabolite(self._bundle()))
        assert 'GEM:Human-GEM' in set(main['resource'])

    def test_gem_transporter_excluded(self):
        main = _bundle_main_nodes(format_enzyme_metabolite(self._bundle()))
        assert 'GEM_transporter:Human-GEM' not in set(main['resource'])

    def test_tcdb_excluded(self):
        main = _bundle_main_nodes(format_enzyme_metabolite(self._bundle()))
        assert 'TCDB' not in set(main['resource'])

    def test_returns_cosmos_bundle(self):
        from omnipath_metabo.datasets.cosmos._bundle import CosmosBundle
        assert isinstance(format_enzyme_metabolite(self._bundle()), CosmosBundle)
