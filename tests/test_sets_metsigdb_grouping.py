"""Grouped responses: by resource, then by set, paged by set.

Cycle 012, user story 2. A flat page makes the reader rebuild a grouping the
data already has. The grouped shape hands it over instead, and pays for that
with a different page unit: `limit` counts **sets**, so no set straddles a page
boundary and every grouped response is coherent on its own.

Two rules the ceiling forces, and they pull in opposite directions:

- a set is never truncated, and
- a response is bounded at 50,000 rows.

Both hold because the page ends early rather than cutting a set, and because a
set larger than the ceiling returns complete and alone. ClassyFire's root class
has 145,937 members, so that second rule is not a hypothetical.

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55435/omnipath \
        uv run --with pytest --with psycopg2-binary \
        pytest tests/test_sets_metsigdb_grouping.py -v
"""

from __future__ import annotations

import os

import pytest

from omnipath_metabo.server.sets._metsigdb_projection import (
    MEMBER_FIELDS,
    SET_FIELDS,
)
from omnipath_metabo.server.sets._metsigdb_query import (
    ROW_CEILING,
    MetSigDBQuery,
    count_sets,
    fetch_groups,
)

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DB_URL,
    reason='No OMNIPATH_DB_URL/DATABASE_URL; grouping is measured against real sets',
)

# The two ClassyFire sets that exceed the ceiling on their own.
ROOT_CLASS = 'CHEMONTID:9999999'
ROOT_CLASS_SIZE = 145_937


@pytest.fixture(scope='module')
def conn():
    psycopg2 = pytest.importorskip('psycopg2')
    connection = psycopg2.connect(DB_URL)
    connection.autocommit = True
    try:
        yield connection
    finally:
        connection.close()


def _sets(groups):
    """Every set in a grouped response, flattened, in the order it arrived."""
    return [one_set for group in groups for one_set in group['sets']]


# ------------------------------------------------------------------ the shape


def test_a_group_carries_its_resource_and_its_sets(conn):
    groups, _ = fetch_groups(conn, MetSigDBQuery(resource=('MACdb',), limit=3))

    assert groups
    assert {group['resource'] for group in groups} == {'MACdb'}
    assert all(group['sets'] for group in groups)


def test_the_set_side_is_written_once_per_set(conn):
    """Not repeated on every member: that would undo the projection work."""
    groups, _ = fetch_groups(conn, MetSigDBQuery(resource=('MACdb',), limit=1))
    one_set = _sets(groups)[0]

    assert 'set_label' in one_set
    assert 'set_size' in one_set
    member = one_set['members'][0]
    for field in SET_FIELDS:
        assert field not in member, f'{field} belongs to the set, not the member'
    assert 'entity' in member
    assert set(member) <= set(MEMBER_FIELDS)


# ------------------------------------------------------------------- the page


def test_a_page_carries_at_most_the_requested_number_of_sets(conn):
    groups, _ = fetch_groups(conn, MetSigDBQuery(resource=('MACdb',), limit=5))
    assert len(_sets(groups)) <= 5


def test_no_set_in_a_response_is_partially_populated(conn):
    """The rule the whole page unit exists to keep."""
    groups, _ = fetch_groups(conn, MetSigDBQuery(resource=('MACdb',), limit=10))

    for one_set in _sets(groups):
        assert one_set['returned'] == one_set['set_size']


def test_the_grouped_total_counts_sets_not_rows(conn):
    """269 MACdb sets, not its 20,291 memberships."""
    assert count_sets(conn, MetSigDBQuery(resource=('MACdb',))) == 269


def test_paging_returns_each_set_exactly_once(conn):
    """Neither repeating nor skipping, across a full walk."""
    seen: list[str] = []
    offset = 0
    while True:
        groups, has_more = fetch_groups(
            conn, MetSigDBQuery(resource=('KEGG',), limit=40, offset=offset)
        )
        seen.extend(one_set['set'] for one_set in _sets(groups))
        if not has_more:
            break
        offset += 40

    assert len(seen) == len(set(seen)), 'a set came back on two pages'
    assert len(seen) == 176, 'KEGG publishes 176 sets'


def test_groups_arrive_in_resource_then_set_order(conn):
    """Stable and deterministic, so paging is repeatable."""
    groups, _ = fetch_groups(conn, MetSigDBQuery(limit=25))

    resources = [group['resource'] for group in groups]
    assert resources == sorted(resources)
    for group in groups:
        ids = [one_set['set'] for one_set in group['sets']]
        assert ids == sorted(ids)


# ---------------------------------------------------------------- the ceiling


def test_a_set_larger_than_the_ceiling_returns_complete_and_alone(conn):
    """FR-015a. The ceiling bounds how many sets share a page, never whether a
    published set can be returned at all."""
    groups, has_more = fetch_groups(
        conn, MetSigDBQuery(set_source_id=ROOT_CLASS, limit=10)
    )
    returned = _sets(groups)

    assert len(returned) == 1
    assert returned[0]['set'] == ROOT_CLASS
    assert returned[0]['returned'] == ROOT_CLASS_SIZE == returned[0]['set_size']
    assert len(returned[0]['members']) == ROOT_CLASS_SIZE


def test_a_page_ends_early_rather_than_exceeding_the_ceiling(conn):
    """Sets are added whole until the next would not fit.

    Offset by one, because ClassyFire's first set in sort order is
    `CHEMONTID:0000000` at 145,635 members, which triggers the oversized-alone
    rule instead and would never exercise this path.
    """
    groups, has_more = fetch_groups(
        conn, MetSigDBQuery(resource=('ClassyFire',), limit=2000, offset=1)
    )
    returned = _sets(groups)

    assert len(returned) < 2000, 'the ceiling never bit'
    assert has_more is True
    assert sum(one_set['returned'] for one_set in returned) <= ROW_CEILING


def test_an_oversized_first_set_is_the_only_exception(conn):
    """The two rules meet here: bounded page, and never a truncated set.

    When the first set of a page exceeds the ceiling on its own, it wins — it is
    returned whole and alone, and the response is over the ceiling. Any other
    outcome would make a published set unreachable.
    """
    groups, has_more = fetch_groups(
        conn, MetSigDBQuery(resource=('ClassyFire',), limit=2000)
    )
    returned = _sets(groups)

    assert len(returned) == 1
    assert returned[0]['set'] == 'CHEMONTID:0000000'
    assert returned[0]['returned'] == returned[0]['set_size'] > ROW_CEILING
    assert has_more is True


def test_the_ceiling_is_the_contract_value():
    assert ROW_CEILING == 50_000


# ------------------------------------------------------------ filtered groups


def test_a_filtered_set_reports_its_population_and_what_it_returned(conn):
    """FR-017. `set_size` is the set's published population; `returned` is what
    this response carries. Under a member-level filter the two differ, and
    reporting only one would make a filtered view look like a shrunken set."""
    probe, _ = fetch_groups(conn, MetSigDBQuery(resource=('MACdb',), limit=1))
    member = _sets(probe)[0]['members'][0]

    groups, _ = fetch_groups(
        conn,
        MetSigDBQuery(
            resource=('MACdb',), metabolite_entity_id=member['entity'], limit=50
        ),
    )
    filtered = _sets(groups)

    assert filtered
    for one_set in filtered:
        assert one_set['returned'] == len(one_set['members'])
        assert one_set['returned'] <= one_set['set_size']
    assert any(
        one_set['returned'] < one_set['set_size'] for one_set in filtered
    ), 'a one-metabolite filter should return fewer members than the set holds'


def test_a_filter_matching_nothing_groups_to_nothing(conn):
    groups, has_more = fetch_groups(
        conn, MetSigDBQuery(set_source_id='no-such-set', limit=10)
    )
    assert groups == []
    assert has_more is False
