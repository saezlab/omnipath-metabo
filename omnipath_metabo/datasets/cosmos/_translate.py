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
ID translation for COSMOS PKN interactions.

Translates all metabolite (source) IDs to ChEBI and all protein (target)
IDs to Ensembl gene IDs (ENSG).

Translation strategies by source id_type:
    - ``'chebi'``: identity (TCDB, SLC are already ChEBI)
    - ``'pubchem'``: UniChem PubChem→ChEBI bulk mapping (STITCH, MRCLinksDB)
    - ``'synonym'``: name→PubChem CID via PubChem REST API, then CID→ChEBI (BRENDA)

Translation strategies by target id_type:
    - ``'uniprot'``: pypath ``uniprot → ensg`` via BioMart (TCDB, SLC, MRCLinksDB)
    - ``'ensp'``: two-hop ``ensp → uniprot → ensg`` (STITCH)
    - ``'genesymbol'``: pypath ``genesymbol → ensg`` (BRENDA fallback proteins)
"""

from __future__ import annotations

__all__ = ['translate_pkn']

import json
import logging
import urllib.parse
import urllib.request
from functools import cache

import pandas as pd

_log = logging.getLogger(__name__)

_UNICHEM_PUBCHEM = 'PubChem'
_UNICHEM_CHEBI = 'ChEBI'
_PUBCHEM_NAME_URL = (
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON'
)


@cache
def _pubchem_to_chebi() -> dict[str, str]:
    """
    Build a PubChem CID → ChEBI ID mapping via UniChem.

    Downloaded once per session and cached in memory. When a PubChem CID maps
    to multiple ChEBI IDs the first one (arbitrary set ordering) is used.

    Returns:
        Dict mapping PubChem CID strings to ChEBI ID strings.
    """

    from pypath.inputs.unichem import unichem_mapping

    raw = unichem_mapping(_UNICHEM_PUBCHEM, _UNICHEM_CHEBI)
    return {cid: next(iter(chebis)) for cid, chebis in raw.items() if chebis}


@cache
def _name_to_chebi(name: str) -> str | None:
    """
    Translate a compound synonym/name to ChEBI.

    Resolves the name to a PubChem CID via the PubChem REST API, then maps
    that CID to ChEBI via the cached UniChem mapping. Results are cached per
    unique name to avoid redundant HTTP calls.

    Args:
        name: Compound name or synonym (e.g. ``'NAD+'``, ``'ATP'``).

    Returns:
        ChEBI ID string (e.g. ``'CHEBI:57540'``), or ``None`` if lookup fails.
    """

    url = _PUBCHEM_NAME_URL.format(urllib.parse.quote(name))

    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read())

        cids = data.get('IdentifierList', {}).get('CID', [])

        if not cids:
            return None

        return _pubchem_to_chebi().get(str(cids[0]))

    except Exception:
        return None


def _to_chebi(source_id: str, id_type: str) -> str | None:
    """
    Translate a metabolite identifier to ChEBI.

    Args:
        source_id: The metabolite identifier.
        id_type: The identifier type (``'chebi'``, ``'pubchem'``, or
            ``'synonym'``).

    Returns:
        ChEBI ID string, or ``None`` if translation is not possible.
    """

    if id_type == 'chebi':
        return source_id

    if id_type == 'pubchem':
        return _pubchem_to_chebi().get(str(source_id))

    if id_type == 'synonym':
        return _name_to_chebi(source_id)

    _log.debug('Unknown metabolite id_type %r, cannot translate to ChEBI.', id_type)
    return None


def _to_ensg(target_id: str, id_type: str, organism: int) -> str | None:
    """
    Translate a protein identifier to an Ensembl gene ID (ENSG).

    Args:
        target_id: The protein identifier.
        id_type: The identifier type (``'uniprot'``, ``'ensp'``, or
            ``'genesymbol'``).
        organism: NCBI taxonomy ID.

    Returns:
        Ensembl gene ID string (``ENSG...``), or ``None`` if translation
        is not possible.
    """

    import pypath.utils.mapping as mapping_mod

    if id_type == 'uniprot':
        result = mapping_mod.map_name(
            target_id,
            'uniprot',
            'ensg',
            ncbi_tax_id=organism,
        )

    elif id_type == 'ensp':
        # ENSP → UniProt → ENSG (direct ensp→ensg is not supported by pypath)
        uniprots = mapping_mod.map_name(
            target_id,
            'ensp',
            'uniprot',
            ncbi_tax_id=organism,
        )
        result = set()

        for u in uniprots:
            result |= mapping_mod.map_name(u, 'uniprot', 'ensg', ncbi_tax_id=organism)

    elif id_type == 'genesymbol':
        result = mapping_mod.map_name(
            target_id,
            'genesymbol',
            'ensg',
            ncbi_tax_id=organism,
        )

    else:
        _log.debug('Unknown protein id_type %r, cannot translate to ENSG.', id_type)
        return None

    if not result:
        return None

    if len(result) > 1:
        _log.debug(
            'Multiple ENSG IDs for %s (%s): %s — using first.',
            target_id,
            id_type,
            result,
        )

    return next(iter(result))


def translate_pkn(df: pd.DataFrame, organism: int = 9606) -> pd.DataFrame:
    """
    Translate COSMOS PKN IDs to unified types.

    Metabolite source IDs are translated to ChEBI.  Protein target IDs are
    translated to Ensembl gene IDs (ENSG).  Rows where either translation
    fails are dropped and a warning is logged.

    Args:
        df:
            PKN DataFrame as returned by :func:`~omnipath_metabo.datasets.cosmos.build`.
        organism:
            NCBI taxonomy ID (default: 9606 for human).

    Returns:
        DataFrame with ``source`` as ChEBI, ``target`` as ENSG, and updated
        ``id_type_a`` / ``id_type_b`` columns.  Index is reset.
    """

    df = df.copy()

    df['source'] = df.apply(
        lambda row: _to_chebi(row['source'], row['id_type_a']),
        axis=1,
    )
    df['target'] = df.apply(
        lambda row: _to_ensg(row['target'], row['id_type_b'], organism),
        axis=1,
    )

    n_failed_source = df['source'].isna().sum()
    n_failed_target = df['target'].isna().sum()

    if n_failed_source:
        _log.warning(
            'Dropped %d rows: metabolite ID could not be translated to ChEBI.',
            n_failed_source,
        )

    if n_failed_target:
        _log.warning(
            'Dropped %d rows: protein ID could not be translated to ENSG.',
            n_failed_target,
        )

    df = df.dropna(subset=['source', 'target'])
    df['id_type_a'] = 'chebi'
    df['id_type_b'] = 'ensg'

    return df.reset_index(drop=True)
