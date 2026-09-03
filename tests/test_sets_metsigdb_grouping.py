"""Grouped responses: by resource, then by set, paged by set.

Empty until User Story 2. The module exists now so the phase that fills it has
somewhere to put its failing tests rather than creating the file and the tests
in one commit, which would make the red step invisible.

What it will assert, from the spec's FR-012 to FR-018:

- a page of `limit=N` carries at most N sets, and never a partial one;
- a set larger than the 50,000-row ceiling returns complete and alone;
- a page ends early rather than exceeding the ceiling, and says `has_more`;
- under a member-level filter a set reports `set_size` as its published
  population and `returned` as what the response carries;
- a full paged walk returns each matching set exactly once.

    DATABASE_URL=... uv run --with pytest --with psycopg2-binary \
        pytest tests/test_sets_metsigdb_grouping.py -v
"""

from __future__ import annotations

import os

import pytest

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DB_URL,
    reason='No OMNIPATH_DB_URL/DATABASE_URL; grouping is measured against real sets',
)
