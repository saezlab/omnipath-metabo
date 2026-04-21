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
Orthology translation for COSMOS PKN bundles.

Translates protein IDs from human to a target organism using
orthologous gene pairs from omnipath-utils (via HTTP or DB).
Applied as a post-processing step after the standard build for
resources that only support human (SLC, Recon3D).
"""

from __future__ import annotations

__all__ = ['translate_bundle_by_orthology']

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

_log = logging.getLogger(__name__)

# Human-only resources whose protein IDs should be translated
# via orthology when building for a non-human organism.
_HUMAN_ONLY_RESOURCES = {'SLC', 'Recon3D'}


def _get_orthology_table(
    source: int,
    target: int,
    id_type: str = 'uniprot',
) -> dict[str, set[str]]:
    """Get the orthology translation table.

    Tries omnipath-utils DB mode first, falls back to omnipath-client HTTP.
    """

    # Try DB mode
    try:
        import os

        db_url = os.environ.get('OMNIPATH_UTILS_DB_URL')

        if db_url:
            from omnipath_utils.orthology import get_table

            return get_table(
                source=source,
                target=target,
                id_type=id_type,
            )
    except Exception:
        pass

    # HTTP fallback
    from omnipath_client.utils import orthology_dict

    return orthology_dict(
        source=source,
        target=target,
        id_type=id_type,
    )


def translate_bundle_by_orthology(
    df: 'pd.DataFrame',
    source_organism: int = 9606,
    target_organism: int = 10090,
    resources: set[str] | None = None,
) -> 'pd.DataFrame':
    """Translate protein IDs in human-only resources to target organism.

    Operates on the PKN DataFrame (after standard ID translation).
    Only modifies rows from resources listed in *resources* (defaults
    to :data:`_HUMAN_ONLY_RESOURCES`).  Protein IDs (UniProt ACs) are
    translated to orthologs; rows without an ortholog are dropped.

    Metabolite IDs (ChEBI) are not translated — they are organism-agnostic.

    Args:
        df: PKN DataFrame with translated IDs (source/target as str or frozenset).
        source_organism: Organism of the human-only resources (default: 9606).
        target_organism: Target organism for orthology translation.
        resources: Resource names to translate (default: SLC, Recon3D).

    Returns:
        DataFrame with translated protein IDs for the specified resources.
        Rows without orthologs are dropped.
    """

    if source_organism == target_organism:
        return df

    resources = resources or _HUMAN_ONLY_RESOURCES

    # Identify rows from human-only resources
    mask = df['resource'].isin(resources)
    n_human = mask.sum()

    if n_human == 0:
        _log.info(
            '[COSMOS] No human-only resource rows to translate by orthology.',
        )
        return df

    _log.info(
        '[COSMOS] Translating %d rows from %s via orthology (%d → %d)...',
        n_human,
        ', '.join(sorted(resources & set(df['resource'].unique()))),
        source_organism,
        target_organism,
    )

    # Get orthology table (human → target)
    orth_table = _get_orthology_table(
        source=source_organism,
        target=target_organism,
        id_type='uniprot',
    )

    if not orth_table:
        _log.warning(
            '[COSMOS] No orthology data for %d → %d. '
            'Human-only resources will be dropped.',
            source_organism,
            target_organism,
        )
        return df[~mask].reset_index(drop=True)

    _log.info(
        '[COSMOS] Orthology table: %d gene pairs (%d → %d)',
        len(orth_table),
        source_organism,
        target_organism,
    )

    # Split: rows to translate vs rows to keep as-is
    df_keep = df[~mask].copy()
    df_translate = df[mask].copy()

    # Translate protein columns
    translated_rows = []

    for _, row in df_translate.iterrows():
        src = row['source']
        tgt = row['target']

        # Translate protein columns (identified by entity type)
        new_src = _translate_id(src, row['source_type'], orth_table)
        new_tgt = _translate_id(tgt, row['target_type'], orth_table)

        if new_src is not None and new_tgt is not None:
            new_row = row.copy()
            new_row['source'] = new_src
            new_row['target'] = new_tgt
            translated_rows.append(new_row)

    if translated_rows:
        import pandas as pd

        df_translated = pd.DataFrame(translated_rows)
        result = pd.concat([df_keep, df_translated], ignore_index=True)
    else:
        result = df_keep.reset_index(drop=True)

    n_translated = len(translated_rows)
    n_dropped = n_human - n_translated
    _log.info(
        '[COSMOS] Orthology: %d/%d rows translated, %d dropped (no ortholog)',
        n_translated,
        n_human,
        n_dropped,
    )

    return result


def _translate_id(
    identifier: str | frozenset,
    entity_type: str,
    orth_table: dict[str, set[str]],
) -> str | frozenset | None:
    """Translate a single identifier via orthology.

    Only protein IDs are translated; metabolite IDs pass through.
    """

    if entity_type != 'protein':
        # Metabolites (ChEBI) are organism-agnostic — pass through
        return identifier

    if isinstance(identifier, frozenset):
        # Translate each component of the frozenset
        translated = set()
        for ac in identifier:
            orthologs = orth_table.get(ac, set())
            translated.update(orthologs)
        return frozenset(translated) if translated else None

    # Single string ID
    orthologs = orth_table.get(str(identifier), set())

    if not orthologs:
        return None

    if len(orthologs) == 1:
        return next(iter(orthologs))

    return frozenset(orthologs)
