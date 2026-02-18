#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos._record module."""

from omnipath_metabo.datasets.cosmos._record import Interaction


class TestInteraction:
    """Tests for the Interaction named tuple."""

    def test_fields(self):
        """Has all expected fields."""

        assert 'source' in Interaction._fields
        assert 'target' in Interaction._fields
        assert 'source_type' in Interaction._fields
        assert 'target_type' in Interaction._fields
        assert 'id_type_a' in Interaction._fields
        assert 'id_type_b' in Interaction._fields
        assert 'interaction_type' in Interaction._fields
        assert 'resource' in Interaction._fields
        assert 'mor' in Interaction._fields
        assert 'locations' in Interaction._fields

    def test_create_minimal(self):
        """Can create with required fields, locations defaults to ()."""

        rec = Interaction(
            source='CID12345',
            target='P12345',
            source_type='small_molecule',
            target_type='protein',
            id_type_a='pubchem',
            id_type_b='uniprot',
            interaction_type='signaling',
            resource='TestDB',
            mor=1,
        )

        assert rec.source == 'CID12345'
        assert rec.locations == ()

    def test_create_with_locations(self):
        """Can create with explicit locations."""

        rec = Interaction(
            source='CHEBI:12345',
            target='Q99999',
            source_type='small_molecule',
            target_type='protein',
            id_type_a='chebi',
            id_type_b='uniprot',
            interaction_type='transport',
            resource='SLC',
            mor=1,
            locations=('e', 'r'),
        )

        assert rec.locations == ('e', 'r')

    def test_is_tuple(self):
        """Interaction is a tuple subclass."""

        rec = Interaction(
            source='a', target='b',
            source_type='small_molecule', target_type='protein',
            id_type_a='x', id_type_b='y',
            interaction_type='t', resource='R', mor=0,
        )

        assert isinstance(rec, tuple)

    def test_indexable(self):
        """Fields accessible by index."""

        rec = Interaction(
            source='a', target='b',
            source_type='small_molecule', target_type='protein',
            id_type_a='x', id_type_b='y',
            interaction_type='t', resource='R', mor=0,
        )

        assert rec[0] == 'a'
        assert rec[1] == 'b'
