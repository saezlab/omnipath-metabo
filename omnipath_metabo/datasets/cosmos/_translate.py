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
    - ``'hmdb'``: UniChem HMDB→ChEBI mapping.  Old 5-digit HMDB IDs
      (e.g. ``HMDB00001``) are silently normalised to the current 7-digit
      format (e.g. ``HMDB0000001``) before lookup.

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

__all__ = ['translate_pkn', '_to_hmdb']

import json
import logging
import urllib.parse
import urllib.request
from functools import cache

import pandas as pd

_log = logging.getLogger(__name__)

_UNICHEM_PUBCHEM = 'PubChem'
_UNICHEM_CHEBI = 'ChEBI'
_UNICHEM_HMDB = 'HMDB'
_PUBCHEM_NAME_URL = (
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON'
)

_HMDB_DIGITS = 7  # current zero-padded width of the numeric part


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


def _normalise_hmdb(hmdb_id: str) -> str:
    """
    Normalise a HMDB identifier to the current 7-digit zero-padded format.

    HMDB IDs were historically distributed in a 5-digit format
    (e.g. ``HMDB00001``).  The current format uses 7 digits
    (e.g. ``HMDB0000001``).  This function converts old-format IDs
    transparently; IDs already in the new format are returned unchanged.

    Args:
        hmdb_id: HMDB identifier string (e.g. ``'HMDB00001'`` or
            ``'HMDB0000001'``).

    Returns:
        HMDB ID with 7-digit zero-padded numeric part
        (e.g. ``'HMDB0000001'``).
    """

    if not hmdb_id.upper().startswith('HMDB'):
        return hmdb_id

    return 'HMDB' + hmdb_id[4:].zfill(_HMDB_DIGITS)


@cache
def _hmdb_to_chebi() -> dict[str, str]:
    """
    Build a HMDB ID → ChEBI ID mapping via UniChem.

    Uses the UniChem HMDB→ChEBI mapping (source 18).  UniChem already
    stores HMDB IDs in the current 7-digit format; keys are normalised
    via :func:`_normalise_hmdb` for safety.  When a HMDB ID maps to
    multiple ChEBI IDs the first one is used.  Downloaded once per
    session and cached in memory.

    Returns:
        Dict mapping normalised HMDB IDs (7-digit) to ChEBI ID strings.
    """

    from pypath.inputs.unichem import unichem_mapping

    raw = unichem_mapping(_UNICHEM_HMDB, _UNICHEM_CHEBI)
    return {
        _normalise_hmdb(hmdb_id): next(iter(chebis))
        for hmdb_id, chebis in raw.items()
        if chebis
    }


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


@cache
def _pubchem_to_hmdb() -> dict[str, str]:
    """
    Build a PubChem CID → HMDB ID mapping via UniChem.

    Downloaded once per session and cached in memory.  When a PubChem CID
    maps to multiple HMDB IDs the first one is used.

    Returns:
        Dict mapping PubChem CID strings to HMDB ID strings.
    """
    from pypath.inputs.unichem import unichem_mapping

    raw = unichem_mapping(_UNICHEM_PUBCHEM, _UNICHEM_HMDB)
    return {cid: next(iter(hmdb_ids)) for cid, hmdb_ids in raw.items() if hmdb_ids}


@cache
def _chebi_to_hmdb() -> dict[str, str]:
    """
    Build a ChEBI ID → HMDB ID mapping via UniChem.

    Downloaded once per session and cached in memory.  When a ChEBI ID maps
    to multiple HMDB IDs the first one is used.

    Returns:
        Dict mapping ChEBI ID strings to HMDB ID strings.
    """
    from pypath.inputs.unichem import unichem_mapping

    raw = unichem_mapping(_UNICHEM_CHEBI, _UNICHEM_HMDB)
    return {chebi: next(iter(hmdb_ids)) for chebi, hmdb_ids in raw.items() if hmdb_ids}


@cache
def _metatlas_to_hmdb(gem: str) -> dict[str, str]:
    """
    Build a MetAtlas base metabolite ID → HMDB mapping for a given GEM.

    Reads the GEM's ``metabolites.tsv`` and extracts the ``metsNoComp`` →
    ``metHMDBID`` correspondence.  Entries with empty HMDB values are skipped.

    Downloaded once per GEM name and cached for the session.

    Args:
        gem: GEM name (e.g. ``'Human-GEM'``).

    Returns:
        Dict mapping base MetAtlas IDs to HMDB IDs, or an empty dict if the
        TSV does not contain an HMDB column.
    """
    from pypath.inputs.metatlas._gem import metatlas_gem_metabolites

    mapping: dict[str, str] = {}

    for row in metatlas_gem_metabolites(gem=gem):
        base_id = row.get('metsNoComp', '')
        hmdb = row.get('metHMDBID', '')
        if base_id and hmdb:
            mapping[base_id] = hmdb

    return mapping


@cache
def _bigg_to_hmdb() -> dict[str, str]:
    """
    Build a BiGG base metabolite ID → HMDB ID mapping from Recon3D.

    Reads the HMDB cross-references embedded in the Recon3D BiGG JSON.
    When a metabolite has multiple HMDB annotations the first entry is used.
    Downloaded once and cached for the session.

    Returns:
        Dict mapping BiGG base IDs to HMDB IDs.
    """
    from pypath.inputs.recon3d._gem import recon3d_metabolites

    mapping: dict[str, str] = {}

    for met in recon3d_metabolites():
        base_id = met.get('base_id', '')
        hmdb_ids = met.get('hmdb', [])

        if base_id and hmdb_ids and base_id not in mapping:
            mapping[base_id] = hmdb_ids[0]

    return mapping


def _to_hmdb(source_id: str, id_type: str, gem: str = '') -> str | None:
    """
    Translate a metabolite identifier to HMDB.

    Parallel to :func:`_to_chebi` but maps to HMDB instead of ChEBI.
    Intended for use in coverage comparison scripts — not used by the
    production ``translate_pkn`` pipeline.

    Args:
        source_id: The metabolite identifier.
        id_type: The identifier type (``'chebi'``, ``'pubchem'``,
            ``'bigg'``, ``'hmdb'``, ``'metatlas'``, or ``'synonym'``).
        gem: GEM name, required when *id_type* is ``'metatlas'``
            (e.g. ``'Human-GEM'``).

    Returns:
        HMDB ID string, or ``None`` if translation is not possible.
    """
    if id_type == 'hmdb':
        return _normalise_hmdb(source_id)

    if id_type == 'chebi':
        return _chebi_to_hmdb().get(source_id)

    if id_type == 'pubchem':
        return _pubchem_to_hmdb().get(str(source_id))

    if id_type == 'bigg':
        return _bigg_to_hmdb().get(source_id)

    if id_type == 'metatlas':
        return _metatlas_to_hmdb(gem).get(source_id) if gem else None

    if id_type == 'synonym':
        # name → pubchem CID → HMDB
        url = _PUBCHEM_NAME_URL.format(urllib.parse.quote(source_id))
        try:
            with urllib.request.urlopen(url, timeout=10) as resp:
                data = json.loads(resp.read())
            cids = data.get('IdentifierList', {}).get('CID', [])
            if not cids:
                return None
            return _pubchem_to_hmdb().get(str(cids[0]))
        except Exception:
            return None

    _log.debug('Unknown metabolite id_type %r, cannot translate to HMDB.', id_type)
    return None


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

    if id_type == 'hmdb':
        return _hmdb_to_chebi().get(_normalise_hmdb(source_id))

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


def _gem_name_series(resource_series: pd.Series) -> pd.Series:
    """
    Extract the GEM name from a ``resource`` column.

    Handles ``'GEM:Human-GEM'`` and ``'GEM_transporter:Human-GEM'`` formats,
    returning the part after the colon.  Non-GEM resources return ``''``.
    """
    gem_mask = resource_series.str.startswith('GEM') & resource_series.str.contains(':')
    result = pd.Series('', index=resource_series.index)
    result[gem_mask] = resource_series[gem_mask].str.split(':', n=1).str[1]
    return result


def _build_metab_mapping(
    id_type: str,
    ids: pd.Series,
    resource_series: pd.Series,
) -> dict[str, str | None]:
    """
    Build a {source_id → ChEBI} mapping for a group of metabolite IDs.

    Called once per ``(id_type, 'small_molecule')`` group.  Uses the cached
    bulk dictionaries wherever possible; falls back to per-name HTTP calls
    only for ``'synonym'`` (BRENDA) entries.

    Args:
        id_type: Metabolite identifier type (e.g. ``'pubchem'``, ``'metatlas'``).
        ids: Series of source IDs (same index as the group slice).
        resource_series: Series of resource names (same index).

    Returns:
        Dict mapping each unique source ID to its ChEBI ID (or ``None``).
    """
    unique_ids = ids.unique()

    if id_type == 'chebi':
        return {uid: uid for uid in unique_ids}

    if id_type == 'pubchem':
        mapping = _pubchem_to_chebi()
        return {uid: mapping.get(str(uid)) for uid in unique_ids}

    if id_type == 'bigg':
        mapping = _bigg_to_chebi()
        return {uid: mapping.get(uid) for uid in unique_ids}

    if id_type == 'hmdb':
        mapping = _hmdb_to_chebi()
        return {uid: mapping.get(_normalise_hmdb(uid)) for uid in unique_ids}

    if id_type == 'metatlas':
        # Each row may come from a different GEM — build per-GEM mappings.
        gem_names = _gem_name_series(resource_series)
        result: dict[str, str | None] = {}
        for uid in unique_ids:
            # Find the GEM name for this metabolite ID (use first occurrence).
            uid_mask = ids == uid
            gem = gem_names[uid_mask].iloc[0] if uid_mask.any() else ''
            result[uid] = _metatlas_to_chebi(gem).get(uid) if gem else None
        return result

    if id_type == 'synonym':
        # One HTTP call per unique name; @cache deduplicates repeated calls.
        return {uid: _name_to_chebi(uid) for uid in unique_ids}

    _log.debug('Unknown metabolite id_type %r, cannot translate to ChEBI.', id_type)
    return {uid: None for uid in unique_ids}


def _build_protein_mapping(
    id_type: str,
    ids: pd.Series,
    organism: int,
) -> dict[str, str | None]:
    """
    Build a {source_id → ENSG} mapping for a group of protein IDs.

    Called once per ``(id_type, 'protein')`` group.  Performs BioMart
    lookups in bulk by iterating unique IDs and mapping through pypath's
    cached mapping tables.

    Args:
        id_type: Protein identifier type (e.g. ``'uniprot'``, ``'ensp'``).
        ids: Series of source IDs (same index as the group slice).
        organism: NCBI taxonomy ID.

    Returns:
        Dict mapping each unique source ID to its ENSG ID (or ``None``).
    """
    import pypath.utils.mapping as mapping_mod

    unique_ids = ids.unique()

    if id_type == 'ensembl':
        return {uid: uid for uid in unique_ids}

    if id_type == 'reaction_id':
        return {uid: uid for uid in unique_ids}

    if id_type == 'entrez':
        bigg_map = _entrez_to_ensg_bigg()
        result: dict[str, str | None] = {}
        for uid in unique_ids:
            ensg = bigg_map.get(uid)
            if not ensg:
                res = mapping_mod.map_name(uid, 'ncbigene', 'ensg',
                                           ncbi_tax_id=organism)
                ensg = next(iter(res)) if res else None
            result[uid] = ensg
        return result

    if id_type == 'uniprot':
        result = {}
        for uid in unique_ids:
            res = mapping_mod.map_name(uid, 'uniprot', 'ensg',
                                       ncbi_tax_id=organism)
            result[uid] = next(iter(res)) if res else None
        return result

    if id_type == 'ensp':
        result = {}
        for uid in unique_ids:
            uniprots = mapping_mod.map_name(uid, 'ensp', 'uniprot',
                                            ncbi_tax_id=organism)
            ensg_set: set[str] = set()
            for u in uniprots:
                ensg_set |= mapping_mod.map_name(u, 'uniprot', 'ensg',
                                                 ncbi_tax_id=organism)
            result[uid] = next(iter(ensg_set)) if ensg_set else None
        return result

    if id_type == 'genesymbol':
        result = {}
        for uid in unique_ids:
            res = mapping_mod.map_name(uid, 'genesymbol', 'ensg',
                                       ncbi_tax_id=organism)
            result[uid] = next(iter(res)) if res else None
        return result

    _log.debug('Unknown protein id_type %r, cannot translate to ENSG.', id_type)
    return {uid: None for uid in unique_ids}


def _translate_column(
    df: pd.DataFrame,
    col: str,
    id_type_col: str,
    entity_type_col: str,
    organism: int,
) -> None:
    """
    Translate one column (``'source'`` or ``'target'``) in-place.

    Groups rows by ``(id_type, entity_type)`` and does bulk dict lookups
    via ``Series.map(dict)`` instead of per-row Python calls.

    Args:
        df: PKN DataFrame (modified in-place).
        col: Column to translate (``'source'`` or ``'target'``).
        id_type_col: Column holding the id_type for *col*.
        entity_type_col: Column holding the entity_type for *col*.
        organism: NCBI taxonomy ID.
    """
    for (id_type, entity_type), idx in df.groupby(
        [id_type_col, entity_type_col]
    ).groups.items():
        ids = df.loc[idx, col]

        if entity_type == 'small_molecule':
            mapping = _build_metab_mapping(id_type, ids, df.loc[idx, 'resource'])
        elif entity_type == 'protein':
            mapping = _build_protein_mapping(id_type, ids, organism)
        else:
            continue

        df.loc[idx, col] = ids.map(mapping)


def translate_pkn(df: pd.DataFrame, organism: int = 9606) -> pd.DataFrame:
    """
    Translate COSMOS PKN IDs to unified types.

    Metabolite source IDs are translated to ChEBI.  Protein target IDs are
    translated to Ensembl gene IDs (ENSG).  Rows where either translation
    fails are dropped and a warning is logged.

    Translation is vectorised: rows are grouped by ``(id_type, entity_type)``
    and bulk dict lookups are used via ``Series.map`` rather than per-row
    Python function calls.  This is significantly faster for large PKNs.

    Args:
        df:
            PKN DataFrame as returned by
            :func:`~omnipath_metabo.datasets.cosmos.build`.
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
    _translate_column(df, 'source', 'id_type_a', 'source_type', organism)
    _translate_column(df, 'target', 'id_type_b', 'target_type', organism)

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
