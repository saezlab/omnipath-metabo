"""Query-layer tests for /sets/metsigdb (cycle 010).

The query layer owns the SQL: which rows come back for a filter set, in which
order, and how many. It owns nothing about their shape.

Construction is tested without a database. Filter semantics need the built
substrate::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55435/omnipath \
        uv run --with pytest --with psycopg2-binary \
        pytest tests/test_sets_metsigdb_query.py -v
"""

from __future__ import annotations

import os
import re

import pytest

from omnipath_metabo.server.sets._metsigdb_projection import resolve_fields
from omnipath_metabo.server.sets._metsigdb_query import (
    DEFAULT_LIMIT,
    MAX_LIMIT,
    MetSigDBQuery,
    build_query,
    fetch,
)

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')


@pytest.fixture(scope='module')
def conn():
    psycopg2 = pytest.importorskip('psycopg2')
    if not DB_URL:
        pytest.skip('no OMNIPATH_DB_URL/DATABASE_URL')
    connection = psycopg2.connect(DB_URL)
    connection.autocommit = True
    try:
        yield connection
    finally:
        connection.close()


# ----------------------------------------------------------------- construction


def test_an_empty_spec_selects_the_whole_table_bounded():
    sql, params = build_query(MetSigDBQuery())
    assert 'WHERE' not in sql
    assert params['limit'] == DEFAULT_LIMIT
    assert params['offset'] == 0


def test_order_is_the_row_identity_so_paging_is_stable():
    sql, _ = build_query(MetSigDBQuery())
    assert 'ORDER BY resource, set_source_id, metabolite_entity_id' in sql


def test_every_query_is_bounded():
    """No filter combination may return the whole substrate in one response."""
    for spec in (MetSigDBQuery(), MetSigDBQuery(resource=('KEGG',))):
        sql, params = build_query(spec)
        assert 'LIMIT' in sql and 'OFFSET' in sql
        assert params['limit'] <= MAX_LIMIT


def test_resource_and_set_type_accept_several_values():
    """Exact equality, so `resource` keeps the primary key's index condition."""
    sql, params = build_query(
        MetSigDBQuery(resource=('KEGG', 'Reactome'), set_type=('pathway',))
    )
    assert 'resource = ANY(%(resource)s)' in sql
    assert 'set_type = ANY(%(set_type)s)' in sql
    assert params['resource'] == ['KEGG', 'Reactome']


def test_a_vocabulary_value_is_canonicalised_not_lowered():
    """Case-insensitivity is a normalisation at the edge, not a SQL expression.

    Lowering the column costs the primary key, which `resource` leads. Measured
    on this substrate, one five-row WikiPathways page went from 6 buffers to
    91,442 because the scan discarded 3,503,370 ClassyFire rows first.
    """
    sql, params = build_query(MetSigDBQuery(resource=('wikipathways',)))
    assert 'lower(' not in sql
    assert params['resource'] == ['WikiPathways']


def test_scalar_filters_use_equality():
    sql, params = build_query(
        MetSigDBQuery(set_source_id='R-HSA-1059683', organism=9606)
    )
    assert 'set_source_id = %(set_source_id)s' in sql
    assert 'organism = %(organism)s' in sql
    assert params['organism'] == 9606


def test_only_shared_columns_are_filterable():
    """The query layer offers the published filters and no others.

    `set_sub_type` joined them when MACdb's trait type and KEGG's overview maps
    became columns. The six identifier columns are still published and still
    not filterable.
    """
    # `fields` joined the dataclass in cycle 012 and is not a filter: it names
    # the response projection, which decides the columns the query reads. It is
    # excluded here for the same reason `limit` and `offset` are.
    filterable = {
        field
        for field in MetSigDBQuery.__dataclass_fields__
        if field not in {'limit', 'offset', 'fields'}
    }
    assert filterable == {
        'resource',
        'set_type',
        'set_sub_type',
        'organism',
        # Cycle 012 renamed two of them and kept the old names as aliases, so
        # both spellings are filters until a later cycle retires the originals.
        'set',
        'entity',
        'set_source_id',
        'metabolite_entity_id',
    }


def test_the_sub_type_filter_narrows_rows(conn):
    """MACdb calls every trait a disease; 116 of its 269 are not."""
    rows = fetch(conn, MetSigDBQuery(set_sub_type=('cancer',), limit=50))
    assert rows
    assert {row['resource'] for row in rows} == {'MACdb'}
    assert {row['set_sub_type'] for row in rows} == {'cancer'}


def test_a_page_knows_whether_more_follows(conn):
    from omnipath_metabo.server.sets._metsigdb_query import fetch_page

    rows, has_more = fetch_page(conn, MetSigDBQuery(resource=('KEGG',), limit=10))
    assert len(rows) == 10 and has_more

    rows, has_more = fetch_page(
        conn, MetSigDBQuery(set_source_id='rn00270', limit=1000)
    )
    assert len(rows) == 16 and not has_more


def test_the_count_ignores_the_page(conn):
    from omnipath_metabo.server.sets._metsigdb_query import count

    spec = MetSigDBQuery(resource=('KEGG',), limit=10)
    assert count(conn, spec) == 4969


# -------------------------------------------------------------------- semantics


def test_resource_filter_narrows_rows(conn):
    rows = fetch(conn, MetSigDBQuery(resource=('KEGG',), limit=50))
    assert rows
    assert {row['resource'] for row in rows} == {'KEGG'}


def test_set_type_filter_narrows_rows(conn):
    rows = fetch(conn, MetSigDBQuery(set_type=('disease',), limit=50))
    assert rows
    assert {row['set_type'] for row in rows} == {'disease'}
    assert {row['resource'] for row in rows} == {'MACdb'}


def test_organism_filter_matches_explicit_values_only(conn):
    """Null-organism rows stay in the dataset and never match the filter.

    `organism` is filterable but no longer in the default projection, so the
    assertion has to ask for the column it reads. Cycle 012 moved it to the
    nine `fields` reaches; filtering on a field and returning it are now two
    separate requests.
    """
    rows = fetch(
        conn,
        MetSigDBQuery(
            organism=9606, fields=resolve_fields(['organism']), limit=50
        ),
    )
    assert rows
    assert all(row['organism'] == 9606 for row in rows)
    assert {row['resource'] for row in rows} == {'Reactome'}


def test_set_source_id_filter_returns_one_set(conn):
    rows = fetch(conn, MetSigDBQuery(set_source_id='R-HSA-1059683', limit=500))
    assert rows
    assert {row['set_source_id'] for row in rows} == {'R-HSA-1059683'}
    assert rows[0]['set_size'] == len(rows)


def test_metabolite_filter_can_span_resources(conn):
    """One metabolite, every set it belongs to, wherever the set came from."""
    probe = fetch(conn, MetSigDBQuery(resource=('MACdb',), limit=1))[0]
    rows = fetch(
        conn,
        MetSigDBQuery(
            metabolite_entity_id=probe['metabolite_entity_id'], limit=MAX_LIMIT
        ),
    )
    assert len(rows) >= 1
    assert {row['metabolite_entity_id'] for row in rows} == {
        probe['metabolite_entity_id']
    }


def test_an_impossible_filter_pair_returns_no_rows(conn):
    """MACdb is the only disease producer, so KEGG has no disease sets."""
    assert fetch(conn, MetSigDBQuery(resource=('KEGG',), set_type=('disease',))) == []


def test_paging_neither_repeats_nor_skips(conn):
    first = fetch(conn, MetSigDBQuery(resource=('KEGG',), limit=100, offset=0))
    second = fetch(conn, MetSigDBQuery(resource=('KEGG',), limit=100, offset=100))
    whole = fetch(conn, MetSigDBQuery(resource=('KEGG',), limit=200, offset=0))

    def identity(rows):
        return [
            (row['resource'], row['set_source_id'], row['metabolite_entity_id'])
            for row in rows
        ]

    assert identity(first) + identity(second) == identity(whole)
    assert len(set(identity(whole))) == len(whole)


# ------------------------------------------------ identifier lookup (US3)


def test_a_value_is_recognised_by_its_shape():
    """A caller pastes what they hold; the route works out what it is.

    A type parameter would make the easy case harder and the hard case no
    easier: the namespace is often the thing the caller is asking about.
    """
    from omnipath_metabo.server.sets._metsigdb_query import identifier_columns

    assert identifier_columns('f0745119-5530-f7f4-a3df-77cf1d52f2fb') == (
        'metabolite_entity_id',
    )
    assert identifier_columns('HMDB0011757') == ('hmdb',)
    assert identifier_columns('CHEBI:21565') == ('chebi',)
    assert identifier_columns('C00001') == ('kegg',)
    assert identifier_columns('IHYJTAOFMMMOPX-LURJTMIESA-N') == ('inchikey',)


def test_a_bare_integer_is_chebi_or_pubchem():
    """Both namespaces are bare integers, so the value alone cannot decide.

    The union is the honest answer. Guessing one would silently drop the other's
    memberships from a result the caller believes is complete.
    """
    from omnipath_metabo.server.sets._metsigdb_query import identifier_columns

    assert set(identifier_columns('21565')) == {'chebi', 'pubchem'}


def test_an_unrecognisable_value_is_refused():
    from omnipath_metabo.server.sets._metsigdb_query import identifier_columns

    with pytest.raises(ValueError) as excinfo:
        identifier_columns('CCO')  # a SMILES, which has no recognisable shape

    assert 'CCO' in str(excinfo.value)


def test_smiles_is_not_an_accepted_entity_value():
    """Published and reachable through `fields`, never a way in.

    It has no distinguishing shape — any string could be a SMILES — so a rule
    for it would swallow every value the other rules declined. It is also the
    one identifier column left unindexed.
    """
    from omnipath_metabo.server.sets._metsigdb_query import IDENTIFIER_COLUMNS

    assert 'smiles' not in IDENTIFIER_COLUMNS


def test_the_entity_filter_finds_a_metabolite_by_external_id(conn):
    """The blocker cycle 010 left: an identifier the consumer already holds."""
    rows = fetch(conn, MetSigDBQuery(entity=('HMDB0011757',), limit=20))
    assert rows
    assert {row['hmdb'] for row in rows} == {'HMDB0011757'}


def test_the_entity_filter_still_takes_an_internal_id(conn):
    probe = fetch(conn, MetSigDBQuery(resource=('MACdb',), limit=1))[0]
    rows = fetch(
        conn,
        MetSigDBQuery(entity=(str(probe['metabolite_entity_id']),), limit=MAX_LIMIT),
    )
    assert {str(row['metabolite_entity_id']) for row in rows} == {
        str(probe['metabolite_entity_id'])
    }


def test_a_mixed_list_returns_the_union(conn):
    probe = fetch(conn, MetSigDBQuery(resource=('MACdb',), limit=1))[0]
    internal = str(probe['metabolite_entity_id'])

    both = fetch(
        conn,
        MetSigDBQuery(entity=(internal, 'HMDB0011757'), limit=MAX_LIMIT),
    )
    assert len(both) >= 1
    assert {str(row['metabolite_entity_id']) for row in both} >= {internal}


def test_a_set_value_is_never_read_as_an_identifier(conn):
    """MACdb set ids are bare integers and collide with ChEBI ids.

    Cycle 010's first set-name measurement reported 645 MACdb sets instead of
    269 for exactly this reason.
    """
    rows = fetch(conn, MetSigDBQuery(set=('1',), limit=50))
    assert rows
    assert {row['resource'] for row in rows} == {'MACdb'}
    assert {row['set_source_id'] for row in rows} == {'1'}


def test_identifier_lookup_is_resource_dependent(conn):
    """Reactome carries no HMDB identifier at all, so this is empty by design.

    The contract publishes the coverage table so an empty result here reads as
    documented sparsity rather than as a fault.
    """
    rows = fetch(
        conn,
        MetSigDBQuery(entity=('HMDB0011757',), resource=('Reactome',), limit=10),
    )
    assert rows == []


def test_the_query_reads_one_table_and_no_other():
    """The single-table rule cycle 010 rests on.

    Joining `entity_identifier_lookup` would make identifier lookup
    resource-independent and cost no index, and it was rejected: a result could
    then no longer be reproduced from the published rows. This is the test that
    notices if it comes back.
    """
    from omnipath_metabo.server.sets._metsigdb_query import (
        build_query,
        count_query,
        members_query,
        set_page_query,
    )

    specs = (
        MetSigDBQuery(),
        MetSigDBQuery(entity=('HMDB0011757', '21565'), limit=10),
        MetSigDBQuery(resource=('MACdb',), set=('1',)),
    )
    builders = (build_query, count_query, set_page_query)
    statements = [sql for spec in specs for sql, _ in (b(spec) for b in builders)]
    statements.append(members_query(specs[0], [('MACdb', '1')])[0])

    # Stated as the rule rather than by parsing SQL: the only schema-qualified
    # table is the substrate, and none of the core tables an identifier join
    # would reach for is mentioned at all.
    forbidden = (
        'entity_identifier_lookup',
        'identifier_evidence',
        'entity_evidence_resolution',
        'entity_ontology_term',
        'relation_evidence',
    )
    for sql in statements:
        qualified = set(re.findall(r'\bpublic\.\w+', sql))
        assert qualified <= {'public.metsigdb_membership'}, sql
        for table in forbidden:
            assert table not in sql, f'{table} reached from the serving layer'


# ------------------------------------------------------- no full scans (T050)


@pytest.mark.parametrize(
    ('label', 'spec'),
    [
        ('entity by hmdb', MetSigDBQuery(entity=('HMDB0011757',))),
        ('entity by inchikey', MetSigDBQuery(entity=('IHYJTAOFMMMOPX-LURJTMIESA-N',))),
        ('entity by kegg', MetSigDBQuery(entity=('C00001',))),
        ('entity by chebi prefixed', MetSigDBQuery(entity=('CHEBI:21565',))),
        ('entity by bare integer', MetSigDBQuery(entity=('21565',))),
        ('entity by uuid', MetSigDBQuery(
            entity=('f0745119-5530-f7f4-a3df-77cf1d52f2fb',))),
        ('entity mixed list', MetSigDBQuery(entity=('HMDB0011757', '21565'))),
        ('set', MetSigDBQuery(set=('R-HSA-1059683',))),
        ('set with resource', MetSigDBQuery(resource=('MACdb',), set=('1',))),
        ('resource, canonical', MetSigDBQuery(resource=('WikiPathways',))),
        ('resource, other case', MetSigDBQuery(resource=('wikipathways',))),
        ('set_type', MetSigDBQuery(set_type=('disease',))),
        ('set_sub_type', MetSigDBQuery(set_sub_type=('cancer',))),
        ('organism', MetSigDBQuery(organism=9606)),
        ('sub_type and resource', MetSigDBQuery(
            resource=('KEGG',), set_sub_type=('overview_map',))),
    ],
)
def test_no_accepted_filter_scans_the_whole_substrate(conn, label, spec):
    """FR-011 and SC-006, held by a test rather than by a one-off EXPLAIN.

    Every shape here is one a request can produce, and every one is paged: the
    route never issues an unbounded select.

    **A sequential scan is not the only way to read the whole table.** The
    `lower()` predicates this cycle briefly shipped planned as an *Index Only
    Scan* — and then discarded 3,503,370 rows in a filter to return five. Only
    the suite's runtime caught it. So the assertion is on rows discarded, not on
    the scan node: an index that the planner uses for ordering while filtering
    everything out is a full scan wearing a better name.

    `count(*)` is deliberately not covered. A caller asking `total=true` for a
    filter matching most of the substrate has asked for work no index can avoid,
    and the contract already says that costs.
    """
    sql, params = build_query(spec)
    with conn.cursor() as cur:
        cur.execute(f'EXPLAIN (ANALYZE, TIMING OFF, FORMAT TEXT) {sql}', params)
        plan = '\n'.join(row[0] for row in cur.fetchall())

    assert 'Seq Scan on metsigdb_membership' not in plan, (
        f'{label} scans the whole substrate:\n{plan}'
    )

    discarded = sum(
        int(match) for match in re.findall(r'Rows Removed by Filter: (\d+)', plan)
    )
    assert discarded < 100_000, (
        f'{label} discarded {discarded:,} rows to fill one page. The plan uses '
        f'an index for order and filters the rest, which reads the table:\n{plan}'
    )
