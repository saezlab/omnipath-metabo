#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos.location module."""

import pytest

from omnipath_metabo.datasets.cosmos.location import (
    locations_to_abbreviations,
    slc_locations,
    slc_routes,
    tcdb_locations,
    tcdb_routes,
)


class TestLocationsToAbbreviations:
    """Tests for locations_to_abbreviations function."""

    def test_single_location(self):
        """Map a single location to abbreviation."""

        from collections import namedtuple
        Loc = namedtuple('UniprotLocation', ['location', 'features'])

        mapping = tcdb_locations()
        locations = [Loc('Cytoplasm', None)]

        result = locations_to_abbreviations(locations, mapping)

        assert result == {'c'}

    def test_multiple_locations(self):
        """Map multiple locations to abbreviations."""

        from collections import namedtuple
        Loc = namedtuple('UniprotLocation', ['location', 'features'])

        mapping = tcdb_locations()
        locations = [
            Loc('Cytoplasm', None),
            Loc('Cell membrane', ('Multi-pass',)),
            Loc('Mitochondrion', None),
        ]

        result = locations_to_abbreviations(locations, mapping)

        assert 'c' in result
        assert 'e' in result
        assert 'm' in result

    def test_unknown_location(self):
        """Unknown locations are ignored."""

        from collections import namedtuple
        Loc = namedtuple('UniprotLocation', ['location', 'features'])

        mapping = tcdb_locations()
        locations = [Loc('Unknown compartment', None)]

        result = locations_to_abbreviations(locations, mapping)

        assert result == set()

    def test_mixed_known_unknown(self):
        """Mix of known and unknown locations."""

        from collections import namedtuple
        Loc = namedtuple('UniprotLocation', ['location', 'features'])

        mapping = tcdb_locations()
        locations = [
            Loc('Cytoplasm', None),
            Loc('Unknown compartment', None),
            Loc('Nucleus', None),
        ]

        result = locations_to_abbreviations(locations, mapping)

        assert 'c' in result
        assert 'n' in result
        assert len(result) == 2

    def test_empty_locations(self):
        """Empty locations list returns empty set."""

        mapping = tcdb_locations()
        locations = []

        result = locations_to_abbreviations(locations, mapping)

        assert result == set()

    def test_returns_set(self):
        """Result is a set type."""

        from collections import namedtuple
        Loc = namedtuple('UniprotLocation', ['location', 'features'])

        mapping = tcdb_locations()
        locations = [Loc('Cytoplasm', None)]

        result = locations_to_abbreviations(locations, mapping)

        assert isinstance(result, set)


class TestTcdbLocations:
    """Tests for tcdb_locations function."""

    def test_returns_dict(self):
        """Returns a dict."""

        result = tcdb_locations()

        assert isinstance(result, dict)

    def test_cytoplasm_mapping(self):
        """Cytoplasm maps to 'c'."""

        result = tcdb_locations()

        assert result['Cytoplasm'] == 'c'

    def test_mitochondrion_mapping(self):
        """Mitochondrion maps to 'm'."""

        result = tcdb_locations()

        assert result['Mitochondrion'] == 'm'

    def test_membrane_mapping(self):
        """Cell membrane maps to 'e'."""

        result = tcdb_locations()

        assert result['Cell membrane'] == 'e'

    def test_has_multiple_entries(self):
        """Contains multiple location entries."""

        result = tcdb_locations()

        assert len(result) > 50


class TestTcdbRoutes:
    """Tests for tcdb_routes function."""

    def test_returns_list(self):
        """Returns a list of tuples."""

        result = tcdb_routes()

        assert isinstance(result, list)
        assert all(isinstance(r, tuple) and len(r) == 2 for r in result)

    def test_routes_to_cytosol(self):
        """Has routes to cytosol."""

        result = tcdb_routes()

        assert any(to == 'c' for _, to in result)

    def test_extracellular_to_cytosol(self):
        """Has e -> c route."""

        result = tcdb_routes()

        assert ('e', 'c') in result


class TestSlcLocations:
    """Tests for slc_locations function."""

    def test_returns_dict(self):
        """Returns a dict."""

        result = slc_locations()

        assert isinstance(result, dict)

    def test_plasma_membrane_mapping(self):
        """Plasma membrane maps to 'e'."""

        result = slc_locations()

        assert result['Plasma membrane'] == 'e'

    def test_lysosome_mapping(self):
        """Lysosome maps to 'l'."""

        result = slc_locations()

        assert result['Lysosome'] == 'l'

    def test_has_compound_locations(self):
        """Contains compound location names with semicolons."""

        result = slc_locations()
        compound_locs = [loc for loc in result if ';' in loc]

        assert len(compound_locs) > 0


class TestSlcRoutes:
    """Tests for slc_routes function."""

    def test_returns_list(self):
        """Returns a list of tuples."""

        result = slc_routes()

        assert isinstance(result, list)
        assert all(isinstance(r, tuple) and len(r) == 2 for r in result)

    def test_cytosol_is_destination(self):
        """Cytosol is a transport destination."""

        result = slc_routes()

        assert any(to == 'c' for _, to in result)


class TestUniprotLocations:
    """Tests for uniprot_locations function (requires pypath)."""

    @pytest.mark.slow
    def test_returns_dict(self):
        """Returns a dict mapping UniProt IDs to sets."""

        from omnipath_metabo.datasets.cosmos.location import uniprot_locations

        result = uniprot_locations(organism=9606, reviewed=True)

        assert isinstance(result, dict)
        # Should have many proteins
        assert len(result) > 1000

    @pytest.mark.slow
    def test_location_has_namedtuple_attributes(self):
        """Location entries have location and features attributes."""

        from omnipath_metabo.datasets.cosmos.location import uniprot_locations

        result = uniprot_locations(organism=9606, reviewed=True)

        # Get first protein with locations
        for uniprot_id, locations in result.items():
            if locations:
                loc = next(iter(locations))
                assert hasattr(loc, 'location')
                assert hasattr(loc, 'features')
                break
