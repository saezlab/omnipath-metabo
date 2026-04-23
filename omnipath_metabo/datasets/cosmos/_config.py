#!/usr/bin/env python

#
# This file is part of the `omnipath_metabo` Python module
#
# Copyright 2026
# Heidelberg University Hospital
#
# File author(s): OmniPath Team (omnipathdb@gmail.com)
#
# Distributed under the BSD-3-Clause license
# See the file `LICENSE` or read a copy at
# https://opensource.org/license/bsd-3-clause
#

"""
COSMOS PKN configuration.

Temporary built-in configuration based on a default YAML file shipped
with the package. This will be replaced by a dedicated solution that
handles configs of multiple packages and loads/merges them from all
possible locations.
"""

from __future__ import annotations

__all__ = ['config', 'default_config']

import copy
from typing import TYPE_CHECKING

import yaml

from .data import data_path

if TYPE_CHECKING:
    from pathlib import Path


def default_config() -> dict:
    """
    Load the built-in default configuration.

    Returns:
        Nested dict with the full default config.
    """

    return _load_yaml(data_path('default_config.yaml'))


def config(
    *args: dict | Path | str,
    **kwargs,
) -> dict:
    """
    Build a configuration by merging defaults with overrides.

    Positional arguments are applied first (dicts or YAML file paths),
    then keyword arguments are merged as a final layer.  Keyword
    arguments that match a resource name are treated as parameter
    overrides for that resource (shorthand for nesting under
    ``resources``).

    The ``resources`` dict controls both which resources are active
    and their parameters.  A resource is active if its key is present
    and its value is not ``False``.  Set a resource to ``False`` to
    disable it.

    Args:
        *args:
            Dicts or paths to YAML files.  Later values take
            precedence over earlier ones.
        **kwargs:
            Config keys merged last.  Resource names are expanded
            into ``resources.<name>``.

    Returns:
        Merged configuration dict.

    Examples::

        # Use all defaults
        cfg = config()

        # Override STITCH score threshold (shorthand)
        cfg = config(stitch={'score_threshold': 400})

        # Run only STITCH and BRENDA
        cfg = config(tcdb=False, slc=False, mrclinksdb=False)

        # Load from a YAML file, then override
        cfg = config('my_config.yaml', stitch={'score_threshold': 400})
    """

    result = default_config()

    for arg in args:
        if isinstance(arg, dict):
            layer = arg
        else:
            layer = _load_yaml(arg)

        _deep_merge(result, layer)

    if kwargs:
        _deep_merge(result, _expand_kwargs(kwargs))

    _auto_select_gem(result)

    return result


def _expand_kwargs(kwargs: dict) -> dict:
    """
    Expand shorthand kwargs into the full config structure.

    Keys matching ``resources`` pass through as-is.  Any other key
    is treated as a resource name shorthand, e.g.
    ``stitch={...}`` becomes ``{'resources': {'stitch': {...}}}``.
    """

    top_keys = {
        'organism', 'resources', 'translate_ids', 'apply_blacklist',
        'orthology_translation',
    }
    expanded: dict = {}
    resource_overrides: dict = {}

    for key, value in kwargs.items():
        if key in top_keys:
            expanded[key] = value
        else:
            resource_overrides[key] = value

    if resource_overrides:
        expanded.setdefault('resources', {})
        _deep_merge(expanded['resources'], resource_overrides)

    return expanded


def _auto_select_gem(cfg: dict) -> None:
    """Auto-select GEM name based on organism if not explicitly set.

    If the user hasn't specified a GEM name in the config and the
    organism has a known default GEM, set it automatically.  If the
    organism has no GEM, disable the GEM resource.
    """

    organism = cfg.get('organism', 9606)
    gem_cfg = cfg.get('resources', {}).get('gem')

    if gem_cfg is False:
        return  # explicitly disabled

    if gem_cfg is None:
        gem_cfg = {}
        cfg.setdefault('resources', {})['gem'] = gem_cfg

    # Only auto-select if the user hasn't specified a GEM name
    if 'gem' not in gem_cfg or gem_cfg['gem'] == 'Human-GEM':
        from omnipath_metabo.datasets.cosmos._organisms import default_gem

        gem_name = default_gem(organism)

        if gem_name:
            gem_cfg['gem'] = gem_name
        elif organism != 9606:
            # No GEM for this organism — disable the resource
            cfg['resources']['gem'] = False

    # Disable human-only GEMs (Recon3D) for non-human organisms.
    # iMM1415 (mouse-only) relies on its own organism guard and is not
    # auto-disabled here, so that config() == default_config() for human builds.
    if organism != 9606:
        for key in ('recon3d', 'recon3d_metabolic'):
            if cfg.get('resources', {}).get(key) is not False:
                cfg.setdefault('resources', {})[key] = False


def _load_yaml(path: Path | str) -> dict:
    """Load a YAML file and return its contents as a dict."""

    with open(path) as f:
        return yaml.safe_load(f) or {}


def _deep_merge(base: dict, override: dict) -> dict:
    """
    Recursively merge *override* into *base* in place.

    For nested dicts, values are merged recursively. For all other
    types (including lists), the override value replaces the base.

    Returns:
        The mutated *base* dict.
    """

    for key, value in override.items():
        if (
            key in base
            and isinstance(base[key], dict)
            and isinstance(value, dict)
        ):
            _deep_merge(base[key], value)
        else:
            base[key] = copy.deepcopy(value)

    return base
