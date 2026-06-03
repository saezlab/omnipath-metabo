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
The omnipath-metabo post-build chemistry layer (Milestones E/F), computed
in-database via the RDKit cartridge against the main OmniPath Postgres.
"""

from __future__ import annotations

from omnipath_metabo.postbuild._postbuild import PostBuildStats, post_build_metabo

__all__ = ['PostBuildStats', 'post_build_metabo']
