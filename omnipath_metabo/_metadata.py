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

"""Package metadata (version, authors, etc)."""

__all__ = ['__version__', '__author__', '__license__']

import importlib.metadata

_FALLBACK_VERSION = '0.0.3'

try:
    __version__ = importlib.metadata.version('omnipath_metabo')
except importlib.metadata.PackageNotFoundError:
    # Package not installed (e.g. running from source checkout)
    __version__ = _FALLBACK_VERSION

__author__ = 'OmniPath Team'
__license__ = 'BSD-3-Clause'
