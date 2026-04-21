#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos._organisms module."""

from omnipath_metabo.datasets.cosmos._organisms import (
    available_gems,
    default_gem,
    needs_orthology,
    organism_resources,
)


class TestDefaultGem:
    """Tests for default_gem()."""

    def test_human(self):
        assert default_gem(9606) == 'Human-GEM'

    def test_mouse(self):
        assert default_gem(10090) == 'Mouse-GEM'

    def test_rat(self):
        assert default_gem(10116) == 'Rat-GEM'

    def test_zebrafish(self):
        assert default_gem(7955) == 'Zebrafish-GEM'

    def test_yeast(self):
        assert default_gem(4932) == 'yeast-GEM'

    def test_unknown_returns_none(self):
        assert default_gem(99999) is None

    def test_ecoli(self):
        assert default_gem(562) == 'Ecoli-GEM'


class TestAvailableGems:
    """Tests for available_gems()."""

    def test_returns_list(self):
        result = available_gems(9606)
        assert isinstance(result, list)

    def test_human_has_one_gem(self):
        result = available_gems(9606)
        assert result == ['Human-GEM']

    def test_mouse_has_one_gem(self):
        result = available_gems(10090)
        assert result == ['Mouse-GEM']

    def test_unknown_returns_empty_list(self):
        result = available_gems(99999)
        assert result == []

    def test_returns_copy(self):
        """Returned list should be a copy, not the internal reference."""
        result = available_gems(9606)
        result.append('Fake-GEM')
        assert available_gems(9606) == ['Human-GEM']


class TestOrganismResources:
    """Tests for organism_resources()."""

    def test_returns_dict(self):
        result = organism_resources(9606)
        assert isinstance(result, dict)

    def test_human_all_direct(self):
        result = organism_resources(9606)
        assert all(result.values()), (
            f'Expected all True for human, got: {result}'
        )

    def test_mouse_direct_resources(self):
        result = organism_resources(10090)
        # tcdb, brenda, stitch: 'all' -> True
        assert result['tcdb'] is True
        assert result['brenda'] is True
        assert result['stitch'] is True

    def test_mouse_mrclinksdb_direct(self):
        result = organism_resources(10090)
        assert result['mrclinksdb'] is True

    def test_mouse_gem_direct(self):
        """Mouse has a GEM -> direct support."""
        result = organism_resources(10090)
        assert result['gem'] is True

    def test_mouse_slc_not_direct(self):
        result = organism_resources(10090)
        assert result['slc'] is False

    def test_mouse_recon3d_not_direct(self):
        result = organism_resources(10090)
        assert result['recon3d'] is False

    def test_rat_gem_direct(self):
        result = organism_resources(10116)
        assert result['gem'] is True

    def test_unknown_organism_gem_not_direct(self):
        """Organism without a GEM -> gem is False."""
        result = organism_resources(99999)
        assert result['gem'] is False

    def test_unknown_organism_universal_resources_direct(self):
        """tcdb, brenda, stitch support all organisms."""
        result = organism_resources(99999)
        assert result['tcdb'] is True
        assert result['brenda'] is True
        assert result['stitch'] is True


class TestNeedsOrthology:
    """Tests for needs_orthology()."""

    def test_human_returns_empty(self):
        result = needs_orthology(9606)
        assert result == []

    def test_mouse_returns_slc_and_recon3d(self):
        result = needs_orthology(10090)
        assert sorted(result) == ['recon3d', 'slc']

    def test_rat_returns_slc_recon3d_mrclinksdb(self):
        """Rat lacks MRCLinksDB support too."""
        result = needs_orthology(10116)
        assert 'slc' in result
        assert 'recon3d' in result
        assert 'mrclinksdb' in result

    def test_unknown_organism(self):
        """Unknown organism needs orthology for slc, recon3d, mrclinksdb, gem."""
        result = needs_orthology(99999)
        assert 'slc' in result
        assert 'recon3d' in result
        assert 'gem' in result
