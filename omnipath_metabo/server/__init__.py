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

"""Litestar web service for the omnipath-metabo package."""

__all__ = ['create_app']


def create_app(*args, **kwargs):
    """Create the Litestar application (lazy import)."""
    from ._app import create_app as _create
    return _create(*args, **kwargs)
