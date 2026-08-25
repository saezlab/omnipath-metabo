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
Metabolite-set route families (cycle 010).

``/sets/metsigdb`` is the first, and in v1 the only one. It is a read-only
additive API over the ``metsigdb_membership`` table that ``omnipath-build``
publishes. The service performs no upstream retrieval and no identifier
remapping on the request path: everything the response needs is already in the
substrate.

Three modules, three responsibilities, so a failure localizes without reading
upstream source data:

- ``_metsigdb_query`` decides which rows come back, in which order, and how
  many.
- ``_metsigdb_projection`` decides their shape.
- ``_metsigdb_route`` validates the request and dispatches.
"""

from __future__ import annotations

__all__ = ['MetSigDBController']

from typing import Any


def __getattr__(name: str) -> Any:
    """Import the controller only when something asks for it.

    The route layer needs Litestar. The query and projection layers need
    neither Litestar nor a running service, and testing one boundary must not
    drag in another. An eager import here would make that impossible.
    """
    if name == 'MetSigDBController':
        from omnipath_metabo.server.sets._metsigdb_route import MetSigDBController

        return MetSigDBController
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
