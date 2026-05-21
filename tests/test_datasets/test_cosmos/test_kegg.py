#!/usr/bin/env python

"""Unit tests for the KEGG enzyme-metabolite COSMOS resource.

All tests are fast (no network I/O).  ``make_kegg_resource`` is patched to
return a mock whose ``reactions.raw()`` yields synthetic reaction dicts.

Module-level helpers ``_compound_id`` and ``_record_to_interactions`` are
tested directly; ``kegg_interactions`` is exercised through its public API.
"""

import sys
from unittest.mock import MagicMock, patch

import pytest

from omnipath_metabo.datasets.cosmos.resources.kegg import (
    _compound_id,
    _record_to_interactions,
    kegg_interactions,
)
from omnipath_metabo.datasets.cosmos._record import Interaction


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_ORGANISM_CODES = {9606: 'hsa', 10090: 'mmu', 10116: 'rno'}


def _rxn(
    rxn_id = 'R00001',
    uniprot = 'P12345',
    substrates = None,
    products = None,
):
    """Build a minimal synthetic KEGG reaction dict."""
    return {
        'Reaction':  rxn_id,
        'UniProt':   uniprot,
        'Substrate': substrates or [],
        'Product':   products or [],
    }


def _cpd(name = 'ATP', kegg_id = 'C00002', chebi = None):
    return {'name': name, 'kegg_id': kegg_id, 'chebi': chebi}


def _run(reaction_dicts, organism = 9606):
    """Call kegg_interactions with a fully mocked kegg_metabolic module.

    Patches ``sys.modules`` so the lazy import inside ``kegg_interactions``
    picks up our mock regardless of whether the real module is installed.
    """
    mock_resource = MagicMock()
    mock_resource.reactions.raw.return_value = iter(reaction_dicts)

    mock_module = MagicMock()
    mock_module.make_kegg_resource.return_value = mock_resource
    mock_module.kegg_organism_code.side_effect = lambda tax: _ORGANISM_CODES.get(tax)

    with patch.dict(sys.modules, {'pypath.inputs_v2.kegg_metabolic': mock_module}):
        return list(kegg_interactions(organism = organism))


# ---------------------------------------------------------------------------
# _compound_id — metabolite ID priority
# ---------------------------------------------------------------------------

class TestCompoundId:

    def test_chebi_takes_priority(self):
        cid, ctype = _compound_id(_cpd(chebi = 'chebi:30616'))
        assert ctype == 'chebi'
        assert cid == 'CHEBI:30616'

    def test_chebi_already_uppercase(self):
        cid, ctype = _compound_id(_cpd(chebi = 'CHEBI:30616'))
        assert cid == 'CHEBI:30616'
        assert ctype == 'chebi'

    def test_chebi_numeric_only(self):
        cid, ctype = _compound_id(_cpd(chebi = 'chebi:12345'))
        assert cid == 'CHEBI:12345'

    def test_kegg_id_without_prefix(self):
        cid, ctype = _compound_id({'name': 'ATP', 'kegg_id': 'C00002', 'chebi': None})
        assert ctype == 'kegg'
        assert cid == 'C00002'

    def test_kegg_id_with_cpd_prefix_stripped(self):
        cid, ctype = _compound_id({'name': 'ATP', 'kegg_id': 'cpd:C00002', 'chebi': None})
        assert ctype == 'kegg'
        assert cid == 'C00002'

    def test_name_fallback_when_no_chebi_or_kegg(self):
        cid, ctype = _compound_id({'name': 'ATP', 'kegg_id': '', 'chebi': None})
        assert ctype == 'synonym'
        assert cid == 'ATP'

    def test_empty_dict_returns_empty_synonym(self):
        cid, ctype = _compound_id({})
        assert ctype == 'synonym'
        assert cid == ''

    def test_chebi_empty_string_falls_through(self):
        cid, ctype = _compound_id({'name': 'Glc', 'kegg_id': 'C00031', 'chebi': ''})
        assert ctype == 'kegg'
        assert cid == 'C00031'


# ---------------------------------------------------------------------------
# _record_to_interactions — edge structure
# ---------------------------------------------------------------------------

class TestInteractionConversion:

    def test_substrate_to_enzyme_edge(self):
        record = _rxn(substrates = [_cpd(kegg_id = 'C00002')])
        edges = list(_record_to_interactions(record))
        sub_edges = [e for e in edges if e.source_type == 'small_molecule']
        assert len(sub_edges) == 1
        e = sub_edges[0]
        assert e.source == 'C00002'
        assert e.target == 'P12345'
        assert e.id_type_a == 'kegg'
        assert e.id_type_b == 'uniprot'

    def test_enzyme_to_product_edge(self):
        record = _rxn(products = [_cpd(kegg_id = 'C00003')])
        edges = list(_record_to_interactions(record))
        prod_edges = [e for e in edges if e.source_type == 'protein']
        assert len(prod_edges) == 1
        e = prod_edges[0]
        assert e.source == 'P12345'
        assert e.target == 'C00003'
        assert e.id_type_a == 'uniprot'
        assert e.id_type_b == 'kegg'

    def test_interaction_type_is_catalysis(self):
        record = _rxn(substrates = [_cpd()], products = [_cpd()])
        edges = list(_record_to_interactions(record))
        assert all(e.interaction_type == 'catalysis' for e in edges)

    def test_resource_is_kegg(self):
        record = _rxn(substrates = [_cpd()], products = [_cpd()])
        edges = list(_record_to_interactions(record))
        assert all(e.resource == 'KEGG' for e in edges)

    def test_mor_is_one(self):
        record = _rxn(substrates = [_cpd()], products = [_cpd()])
        edges = list(_record_to_interactions(record))
        assert all(e.mor == 1 for e in edges)

    def test_empty_uniprot_yields_nothing(self):
        record = _rxn(uniprot = '', substrates = [_cpd()])
        assert list(_record_to_interactions(record)) == []

    def test_whitespace_uniprot_yields_nothing(self):
        record = _rxn(uniprot = '   ', substrates = [_cpd()])
        assert list(_record_to_interactions(record)) == []

    def test_compound_without_id_skipped(self):
        # No chebi, no kegg_id, no name — _compound_id returns ('', ...)
        bad_cpd = {'name': '', 'kegg_id': '', 'chebi': None}
        record = _rxn(substrates = [bad_cpd], products = [bad_cpd])
        assert list(_record_to_interactions(record)) == []

    def test_chebi_compound_passes_through(self):
        cpd = _cpd(chebi = 'chebi:30616')
        record = _rxn(substrates = [cpd])
        edges = list(_record_to_interactions(record))
        sub = [e for e in edges if e.source_type == 'small_molecule']
        assert sub[0].id_type_a == 'chebi'
        assert sub[0].source == 'CHEBI:30616'

    def test_locations_is_empty_tuple(self):
        record = _rxn(substrates = [_cpd()], products = [_cpd()])
        edges = list(_record_to_interactions(record))
        assert all(e.locations == () for e in edges)


# ---------------------------------------------------------------------------
# attrs — reaction_id propagated
# ---------------------------------------------------------------------------

class TestAttrFields:

    def test_reaction_id_in_attrs(self):
        record = _rxn(rxn_id = 'R00099', substrates = [_cpd()], products = [_cpd()])
        edges = list(_record_to_interactions(record))
        assert all(e.attrs['reaction_id'] == 'R00099' for e in edges)

    def test_attrs_shared_dict_same_reaction(self):
        record = _rxn(substrates = [_cpd()], products = [_cpd()])
        edges = list(_record_to_interactions(record))
        assert len(edges) == 2
        # Both edges from the same reaction share the same attrs dict object.
        assert edges[0].attrs is edges[1].attrs


# ---------------------------------------------------------------------------
# Multiple UniProt ACs
# ---------------------------------------------------------------------------

class TestMultiUniprot:

    def test_two_uniprots_double_substrate_edges(self):
        record = _rxn(
            uniprot = 'P11111;Q22222',
            substrates = [_cpd(kegg_id = 'C00002')],
        )
        edges = list(_record_to_interactions(record))
        sub_edges = [e for e in edges if e.source_type == 'small_molecule']
        targets = {e.target for e in sub_edges}
        assert targets == {'P11111', 'Q22222'}

    def test_two_uniprots_double_product_edges(self):
        record = _rxn(
            uniprot = 'P11111;Q22222',
            products = [_cpd(kegg_id = 'C00003')],
        )
        edges = list(_record_to_interactions(record))
        prod_edges = [e for e in edges if e.source_type == 'protein']
        sources = {e.source for e in prod_edges}
        assert sources == {'P11111', 'Q22222'}

    def test_uniprot_with_spaces_trimmed(self):
        record = _rxn(
            uniprot = ' P11111 ; Q22222 ',
            substrates = [_cpd()],
        )
        edges = list(_record_to_interactions(record))
        targets = {e.target for e in edges if e.source_type == 'small_molecule'}
        assert 'P11111' in targets
        assert 'Q22222' in targets

    def test_multiple_compounds_multiple_enzymes(self):
        record = _rxn(
            uniprot = 'P11111;Q22222',
            substrates = [_cpd('ATP', 'C00002'), _cpd('NADH', 'C00004')],
        )
        edges = list(_record_to_interactions(record))
        sub_edges = [e for e in edges if e.source_type == 'small_molecule']
        # 2 compounds × 2 enzymes = 4 edges
        assert len(sub_edges) == 4


# ---------------------------------------------------------------------------
# kegg_interactions — organism routing
# ---------------------------------------------------------------------------

class TestOrganism:

    def test_human_proceeds(self):
        reactions = [_rxn(substrates = [_cpd()])]
        recs = _run(reactions, organism = 9606)
        assert len(recs) > 0

    def test_mouse_proceeds(self):
        reactions = [_rxn(substrates = [_cpd()])]
        recs = _run(reactions, organism = 10090)
        assert len(recs) > 0

    def test_rat_proceeds(self):
        reactions = [_rxn(substrates = [_cpd()])]
        recs = _run(reactions, organism = 10116)
        assert len(recs) > 0

    def test_unmapped_organism_yields_nothing(self, caplog):
        import logging
        with caplog.at_level(logging.WARNING):
            recs = _run([], organism = 9999)
        assert recs == []
        assert '9999' in caplog.text

    def test_empty_reactions_yields_nothing(self):
        assert _run([], organism = 9606) == []

    def test_result_are_interaction_instances(self):
        reactions = [_rxn(substrates = [_cpd()], products = [_cpd()])]
        recs = _run(reactions, organism = 9606)
        assert all(isinstance(r, Interaction) for r in recs)
