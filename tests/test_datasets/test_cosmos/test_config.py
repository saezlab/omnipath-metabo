#!/usr/bin/env python

"""Tests for omnipath_metabo.datasets.cosmos._config module."""

import pytest

from omnipath_metabo.datasets.cosmos._config import (
    _deep_merge,
    _expand_kwargs,
    config,
    default_config,
)


class TestDefaultConfig:
    """Tests for default_config."""

    def test_returns_dict(self):
        cfg = default_config()
        assert isinstance(cfg, dict)

    def test_has_organism(self):
        cfg = default_config()
        assert 'organism' in cfg
        assert cfg['organism'] == 9606

    def test_has_resources(self):
        cfg = default_config()
        assert 'resources' in cfg
        assert isinstance(cfg['resources'], dict)

    def test_default_resources(self):
        """All seven resources present by default."""

        cfg = default_config()
        expected = {'stitch', 'tcdb', 'slc', 'brenda', 'mrclinksdb', 'gem', 'recon3d'}
        assert set(cfg['resources']) == expected

    def test_stitch_defaults(self):
        cfg = default_config()
        stitch = cfg['resources']['stitch']
        assert stitch['score_threshold'] == 700
        assert 'activation' in stitch['mode']
        assert 'inhibition' in stitch['mode']
        assert 'binding' not in stitch['mode']
        assert stitch['a_is_acting'] is True

    def test_independent_copies(self):
        """Each call returns an independent copy."""

        cfg1 = default_config()
        cfg2 = default_config()
        cfg1['organism'] = 10090
        assert cfg2['organism'] == 9606


class TestConfig:
    """Tests for config with overrides."""

    def test_no_overrides(self):
        """Without overrides, equals default."""

        cfg = config()
        default = default_config()
        assert cfg == default

    def test_override_organism(self):
        cfg = config(organism=10090)
        assert cfg['organism'] == 10090

    def test_resource_shorthand(self):
        """Resource name as kwarg merges into resources."""

        cfg = config(stitch={'score_threshold': 400})
        assert cfg['resources']['stitch']['score_threshold'] == 400
        # Other stitch params preserved
        assert 'mode' in cfg['resources']['stitch']

    def test_disable_resource(self):
        """Setting resource to False disables it."""

        cfg = config(tcdb=False)
        assert cfg['resources']['tcdb'] is False

    def test_other_resources_preserved(self):
        """Disabling one resource doesn't affect others."""

        cfg = config(tcdb=False)
        assert isinstance(cfg['resources']['stitch'], dict)
        assert isinstance(cfg['resources']['slc'], dict)

    def test_dict_override(self):
        """Positional dict override works."""

        cfg = config({'organism': 10090})
        assert cfg['organism'] == 10090

    def test_yaml_file_override(self, tmp_path):
        """YAML file override works."""

        yaml_file = tmp_path / 'override.yaml'
        yaml_file.write_text('organism: 10090\n')

        cfg = config(yaml_file)
        assert cfg['organism'] == 10090

    def test_multiple_overrides_order(self):
        """Later overrides take precedence."""

        cfg = config(
            {'organism': 10090},
            organism=9606,
        )
        assert cfg['organism'] == 9606

    def test_resources_kwarg_replaces(self):
        """Explicit resources kwarg deep-merges."""

        cfg = config(resources={'stitch': {'score_threshold': 300}})
        assert cfg['resources']['stitch']['score_threshold'] == 300
        # mode preserved through deep merge
        assert 'mode' in cfg['resources']['stitch']


class TestDeepMerge:
    """Tests for _deep_merge."""

    def test_simple_override(self):
        base = {'a': 1, 'b': 2}
        _deep_merge(base, {'b': 3})
        assert base == {'a': 1, 'b': 3}

    def test_nested_merge(self):
        base = {'a': {'x': 1, 'y': 2}}
        _deep_merge(base, {'a': {'y': 3}})
        assert base == {'a': {'x': 1, 'y': 3}}

    def test_add_new_key(self):
        base = {'a': 1}
        _deep_merge(base, {'b': 2})
        assert base == {'a': 1, 'b': 2}

    def test_list_replaces(self):
        """Lists are replaced, not merged."""

        base = {'a': [1, 2]}
        _deep_merge(base, {'a': [3]})
        assert base == {'a': [3]}

    def test_scalar_replaces_dict(self):
        base = {'a': {'x': 1}}
        _deep_merge(base, {'a': False})
        assert base == {'a': False}

    def test_dict_replaces_scalar(self):
        base = {'a': False}
        _deep_merge(base, {'a': {'x': 1}})
        assert base == {'a': {'x': 1}}


class TestExpandKwargs:
    """Tests for _expand_kwargs."""

    def test_top_level_keys(self):
        result = _expand_kwargs({'organism': 10090})
        assert result == {'organism': 10090}

    def test_resource_shorthand(self):
        result = _expand_kwargs({'stitch': {'score_threshold': 400}})
        assert result == {'resources': {'stitch': {'score_threshold': 400}}}

    def test_mixed(self):
        result = _expand_kwargs({
            'organism': 10090,
            'stitch': {'score_threshold': 400},
        })
        assert result['organism'] == 10090
        assert result['resources']['stitch']['score_threshold'] == 400

    def test_false_resource(self):
        result = _expand_kwargs({'tcdb': False})
        assert result == {'resources': {'tcdb': False}}
