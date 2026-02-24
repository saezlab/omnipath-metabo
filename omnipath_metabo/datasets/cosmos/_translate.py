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

Translates all metabolite IDs to ChEBI and all protein IDs to Ensembl
gene IDs (ENSG).  Translation is direction-aware: GEM resources produce
enzyme→metabolite edges where the protein is the *source*, so
``source_type``/``target_type`` columns are checked rather than assuming
a fixed metabolite-source / protein-target convention.

Translation strategies by metabolite id_type:
    - ``'chebi'``: identity (TCDB, SLC are already ChEBI)
    - ``'pubchem'``: UniChem PubChem→ChEBI bulk mapping (STITCH, MRCLinksDB)
    - ``'synonym'``: name→PubChem CID via PubChem REST API, then CID→ChEBI (BRENDA)
    - ``'metatlas'``: GEM metabolites.tsv ``metsNoComp`` → ``metChEBIID``
      mapping; loaded once per GEM name from the ``resource`` column.
    - ``'bigg'``: Recon3D BiGG base metabolite ID → ChEBI, built from the
      annotation field of the BiGG JSON (``recon3d_metabolites()``).

Translation strategies by protein id_type:
    - ``'uniprot'``: pypath ``uniprot → ensg`` via BioMart (TCDB, SLC, MRCLinksDB)
    - ``'ensp'``: two-hop ``ensp → uniprot → ensg`` (STITCH)
    - ``'genesymbol'``: pypath ``genesymbol → ensg`` (BRENDA fallback proteins)
    - ``'ensembl'``: identity — GEM enzyme IDs are already Ensembl gene IDs
    - ``'entrez'``: Recon3D Entrez Gene ID → ENSG.  Tries the Ensembl gene
      cross-references embedded in the BiGG JSON first; falls back to pypath
      ``ncbigene → ensg`` BioMart mapping.  ``_ATN`` isoform suffixes must
      already be stripped by the caller (``recon3d.py`` does this at parse
      time).
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


@cache
def _bigg_to_chebi() -> dict[str, str]:
    """
    Build a BiGG base metabolite ID → ChEBI ID mapping from Recon3D.

    Reads the ChEBI cross-references embedded in the Recon3D BiGG JSON
    (``recon3d_metabolites()``) — no external API call required.  When a
    metabolite has multiple ChEBI annotations the first entry is used.
    Downloaded once and cached for the session.

    Returns:
        Dict mapping BiGG base IDs (e.g. ``'atp'``) to ChEBI IDs
        (e.g. ``'CHEBI:30616'``).
    """

    from pypath.inputs.recon3d._gem import recon3d_metabolites

    mapping: dict[str, str] = {}

    for met in recon3d_metabolites():
        base_id = met.get('base_id', '')
        chebis = met.get('chebi', [])

        if base_id and chebis and base_id not in mapping:
            mapping[base_id] = chebis[0]

    return mapping


@cache
def _entrez_to_ensg_bigg() -> dict[str, str]:
    """
    Build an Entrez Gene ID → ENSG mapping via BiGG gene symbols.

    The BiGG JSON gene objects carry a ``name`` field containing the HGNC
    gene symbol (e.g. ``'SLC25A21'``).  This function strips the ``_ATN``
    isoform suffix from the gene ``id`` to recover the base Entrez ID, then
    maps the gene symbol to ENSG using pypath's BioMart-backed
    ``genesymbol → ensg`` mapping.

    The genesymbol → ENSG mapping table is downloaded once by pypath and
    cached on disk; subsequent calls are fast dictionary lookups.  The final
    entrez → ENSG dict is cached in memory for the session.

    Returns:
        Dict mapping Entrez Gene ID strings (e.g. ``'89874'``) to ENSG IDs
        (e.g. ``'ENSG00000183032'``).
    """

    import re

    import pypath.utils.mapping as mapping_mod

    from pypath.inputs.recon3d._gem import recon3d_genes

    result: dict[str, str] = {}

    for gene in recon3d_genes():
        raw_id = gene.get('id', '')
        name = gene.get('name', '')

        if not raw_id or not name:
            continue

        entrez = re.sub(r'_AT\d+$', '', raw_id)

        if entrez in result:
            continue

        ensg_set = mapping_mod.map_name(name, 'genesymbol', 'ensg', ncbi_tax_id=9606)

        if ensg_set:
            result[entrez] = next(iter(ensg_set))

    return result


@cache
def _metatlas_to_chebi(gem: str) -> dict[str, str]:
    """
    Build a MetAtlas base metabolite ID → ChEBI mapping for a given GEM.

    Reads the GEM's ``metabolites.tsv`` (via
    :func:`pypath.inputs.metatlas.metatlas_gem_metabolites`) and extracts
    the ``metsNoComp`` (base ID without compartment) → ``metChEBIID``
    correspondence.  Entries with empty ChEBI values are skipped.

    Downloaded once per GEM name and cached for the session.

    Args:
        gem: GEM name (e.g. ``'Human-GEM'``).

    Returns:
        Dict mapping base MetAtlas IDs (e.g. ``'MAM00001'``) to ChEBI IDs
        (e.g. ``'CHEBI:15389'``).
    """

    from pypath.inputs.metatlas._gem import metatlas_gem_metabolites

    mapping: dict[str, str] = {}

    for row in metatlas_gem_metabolites(gem=gem):

        base_id = row.get('metsNoComp', '')
        chebi = row.get('metChEBIID', '')

        if base_id and chebi:
            mapping[base_id] = chebi

    return mapping


def _to_chebi(source_id: str, id_type: str, gem: str = '') -> str | None:
    """
    Translate a metabolite identifier to ChEBI.

    Args:
        source_id: The metabolite identifier.
        id_type: The identifier type (``'chebi'``, ``'pubchem'``,
            ``'synonym'``, or ``'metatlas'``).
        gem: GEM name, required when *id_type* is ``'metatlas'``
            (e.g. ``'Human-GEM'``).

    Returns:
        ChEBI ID string, or ``None`` if translation is not possible.
    """

    if id_type == 'chebi':
        return source_id

    if id_type == 'pubchem':
        return _pubchem_to_chebi().get(str(source_id))

    if id_type == 'synonym':
        return _name_to_chebi(source_id)

    if id_type == 'metatlas':
        return _metatlas_to_chebi(gem).get(source_id)

    if id_type == 'bigg':
        return _bigg_to_chebi().get(source_id)

    _log.debug('Unknown metabolite id_type %r, cannot translate to ChEBI.', id_type)
    return None


def _to_ensg(target_id: str, id_type: str, organism: int) -> str | None:
    """
    Translate a protein identifier to an Ensembl gene ID (ENSG).

    Args:
        target_id: The protein identifier.
        id_type: The identifier type — one of ``'uniprot'``, ``'ensp'``,
            ``'genesymbol'``, ``'ensembl'``, ``'entrez'``, or
            ``'reaction_id'`` (orphan pseudo-enzyme, passed through
            unchanged).
        organism: NCBI taxonomy ID.

    Returns:
        Ensembl gene ID string (``ENSG...``), or ``None`` if translation
        is not possible.
    """

    import pypath.utils.mapping as mapping_mod

    if id_type == 'ensembl':
        # Already an Ensembl gene ID — pass through as-is.
        # GEM enzyme IDs are ENSG... by construction.
        return target_id

    elif id_type == 'reaction_id':
        # Orphan reaction pseudo-enzyme: pass through the reaction ID unchanged.
        return target_id

    elif id_type == 'entrez':
        # Recon3D Entrez Gene IDs (_ATN suffixes already stripped by recon3d.py).
        # Try BiGG-embedded Ensembl cross-references first; fall back to
        # pypath BioMart mapping.
        ensg = _entrez_to_ensg_bigg().get(target_id)

        if ensg:
            return ensg

        result = mapping_mod.map_name(
            target_id,
            'ncbigene',
            'ensg',
            ncbi_tax_id=organism,
        )

        if not result:
            return None

        if len(result) > 1:
            _log.debug(
                'Multiple ENSG IDs for Entrez %s: %s — using first.',
                target_id,
                result,
            )

        return next(iter(result))

    elif id_type == 'uniprot':
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

    # Translation is direction-aware: source/target roles vary by resource.
    # GEM produces enzyme→metabolite edges (protein as source, metabolite as
    # target) in addition to the more common metabolite→protein direction.
    #
    # For MetAtlas metabolite IDs the GEM name is extracted from the
    # ``resource`` column (format: ``'GEM:<gem-name>'``).

    def _gem_name(resource: str) -> str:
        """Extract GEM name from 'GEM:Human-GEM' or 'GEM_transporter:Human-GEM'."""
        return resource.split(':', 1)[-1] if resource.startswith('GEM') and ':' in resource else ''

    def _translate_entity(entity_id, id_type, entity_type, organism, resource):

        if entity_type == 'small_molecule':
            return _to_chebi(entity_id, id_type, gem=_gem_name(resource))

        if entity_type == 'protein':
            return _to_ensg(entity_id, id_type, organism)

        return None

    df['source'] = df.apply(
        lambda row: _translate_entity(
            row['source'], row['id_type_a'], row['source_type'],
            organism, row['resource'],
        ),
        axis=1,
    )
    df['target'] = df.apply(
        lambda row: _translate_entity(
            row['target'], row['id_type_b'], row['target_type'],
            organism, row['resource'],
        ),
        axis=1,
    )

    n_failed_source = df['source'].isna().sum()
    n_failed_target = df['target'].isna().sum()

    if n_failed_source:
        _log.warning(
            'Dropped %d rows: source ID could not be translated.',
            n_failed_source,
        )

    if n_failed_target:
        _log.warning(
            'Dropped %d rows: target ID could not be translated.',
            n_failed_target,
        )

    df = df.dropna(subset=['source', 'target']).copy()

    _type_to_id = {'small_molecule': 'chebi', 'protein': 'ensg'}

    # Preserve 'reaction_id' entries (orphan pseudo-enzymes) — only update
    # rows where the id_type is a translatable type.
    mask_a = df['id_type_a'] != 'reaction_id'
    mask_b = df['id_type_b'] != 'reaction_id'

    df.loc[mask_a, 'id_type_a'] = df.loc[mask_a, 'source_type'].map(_type_to_id)
    df.loc[mask_b, 'id_type_b'] = df.loc[mask_b, 'target_type'].map(_type_to_id)

    return df.reset_index(drop=True)
