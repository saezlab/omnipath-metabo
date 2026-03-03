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
COSMOS PKN output bundle.

:class:`CosmosBundle` is the top-level output object returned by the COSMOS
PKN build pipeline.  It groups the four typed record streams that together
constitute the PKN:

- **network** — the edge list as :class:`CosmosEdge` records
- **metabolites** — ChEBI provenance as :class:`CosmosMetabolite` records
- **proteins** — UniProt provenance as :class:`CosmosProtein` records
- **reactions** — GEM reaction metadata as :class:`CosmosReaction` records

Each component is stored as a plain list of namedtuples.  Convert any
component to a pandas DataFrame with ``pd.DataFrame(bundle.network)``, or
call :meth:`CosmosBundle.to_dataframes` to convert all at once.

Example::

    from omnipath_metabo.datasets.cosmos import build, format_pkn

    pkn_df   = build()
    bundle   = format_pkn(pkn_df)
    net_df   = pd.DataFrame(bundle.network)      # edge table
    met_df   = pd.DataFrame(bundle.metabolites)  # metabolite provenance
    prot_df  = pd.DataFrame(bundle.proteins)     # protein provenance
    rxn_df   = pd.DataFrame(bundle.reactions)    # GEM reaction metadata
"""

from __future__ import annotations

__all__ = ['CosmosBundle']

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd
    from ._record import CosmosEdge, CosmosMetabolite, CosmosProtein, CosmosReaction


@dataclass
class CosmosBundle:
    """
    Full output of the COSMOS PKN build and format pipeline.

    Attributes:
        network:
            Formatted PKN edges as :class:`~._record.CosmosEdge` records.
            Source and target are COSMOS node ID strings
            (e.g. ``'Metab__CHEBI:15422_c'``, ``'Gene1__P00533'``).
        metabolites:
            One :class:`~._record.CosmosMetabolite` per unique metabolite,
            mapping canonical ChEBI IDs back to original source identifiers.
        proteins:
            One :class:`~._record.CosmosProtein` per unique protein,
            mapping canonical UniProt ACs back to original source identifiers.
        reactions:
            One :class:`~._record.CosmosReaction` per unique GEM reaction
            that contributed edges to the network.  Empty for non-GEM
            resources.
    """

    network: list[CosmosEdge] = field(default_factory=list)
    metabolites: list[CosmosMetabolite] = field(default_factory=list)
    proteins: list[CosmosProtein] = field(default_factory=list)
    reactions: list[CosmosReaction] = field(default_factory=list)

    def to_dataframes(self) -> dict[str, pd.DataFrame]:
        """
        Convert all bundle components to pandas DataFrames.

        Returns:
            Dict with keys ``'network'``, ``'metabolites'``, ``'proteins'``,
            ``'reactions'``, each holding the corresponding DataFrame.
            Columns match the fields of the respective namedtuple class.
        """
        import pandas as pd

        return {
            'network': pd.DataFrame(self.network),
            'metabolites': pd.DataFrame(self.metabolites),
            'proteins': pd.DataFrame(self.proteins),
            'reactions': pd.DataFrame(self.reactions),
        }

    def __repr__(self) -> str:
        return (
            f'CosmosBundle('
            f'network={len(self.network)} edges, '
            f'metabolites={len(self.metabolites)}, '
            f'proteins={len(self.proteins)}, '
            f'reactions={len(self.reactions)})'
        )
