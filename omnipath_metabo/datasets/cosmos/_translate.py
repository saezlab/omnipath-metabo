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

Translates all metabolite IDs to ChEBI and all protein IDs to UniProt
accessions.  Translation is direction-aware: GEM resources produce
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
    - ``'uniprot'``: identity passthrough (TCDB, SLC, MRCLinksDB, BRENDA)
    - ``'ensp'``: pypath ``ensp → uniprot`` via BioMart (STITCH)
    - ``'genesymbol'``: pypath ``genesymbol → uniprot`` (BRENDA fallback)
    - ``'ensembl'``: pypath ``ensg → uniprot`` via BioMart (GEM enzyme IDs
      are ENSG; translated to UniProt for cross-network integration)
    - ``'entrez'``: Recon3D Entrez Gene ID → UniProt.  Tries BiGG-embedded
      gene symbol → UniProt first; falls back to pypath
      ``ncbigene → uniprot`` BioMart mapping.  ``_ATN`` isoform suffixes
      must already be stripped by the caller (``recon3d.py`` does this at
      parse time).
"""

from __future__ import annotations

__all__ = ['translate_pkn', '_to_hmdb', '_to_uniprot']

import logging
from functools import cache

import pandas as pd

_log = logging.getLogger(__name__)

_UNICHEM_PUBCHEM = 'PubChem'
_UNICHEM_CHEBI = 'ChEBI'
_UNICHEM_HMDB = 'HMDB'


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

    from pypath.utils.metabo import normalise_hmdb

    return normalise_hmdb(hmdb_id)


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


def _name_to_chebi(name: str) -> str | None:
    """
    Translate a compound synonym/name to ChEBI.

    Delegates to :func:`pypath.inputs.pubchem.pubchem_name_cids`, which
    caches each name→CIDs lookup on disk via pypath's curl infrastructure.
    Subsequent calls for the same name incur no HTTP request.

    Args:
        name: Compound name or synonym (e.g. ``'NAD+'``, ``'ATP'``).

    Returns:
        ChEBI ID string (e.g. ``'CHEBI:57540'``), or ``None`` if lookup fails.
    """

    from pypath.inputs.pubchem import pubchem_name_cids

    cids = pubchem_name_cids(name)

    if not cids:
        return None

    return _pubchem_to_chebi().get(next(iter(cids)))


@cache
def _bigg_to_chebi() -> dict[str, str]:
    """
    Build a BiGG base metabolite ID → ChEBI ID mapping from Recon3D.

    Two-step strategy:

    1. **Direct**: ChEBI cross-references embedded in the Recon3D BiGG
       JSON (``recon3d_metabolites()``).  Covers ~28% of all Recon3D
       metabolites and ~34% of transport metabolites.

    2. **MetaNetX bridge**: For metabolites that have no direct ChEBI
       annotation but do carry a ``metanetx.chemical`` cross-reference in
       the BiGG JSON, the MetaNetX MNXref ``chem_xref.tsv`` is used to
       map ``MNXM*`` → ChEBI.  Adds ~45 additional transport metabolites
       (+2.7 pp coverage).

    Downloaded once and cached for the session.

    Returns:
        Dict mapping BiGG base IDs (e.g. ``'atp'``) to ChEBI IDs
        (e.g. ``'CHEBI:30616'``).
    """

    from pypath.inputs.recon3d._gem import recon3d_metabolites

    mapping: dict[str, str] = {}
    mnx_pending: dict[str, str] = {}  # base_id → MNX ID, for unmapped

    for met in recon3d_metabolites():
        base_id = met.get('base_id', '')

        if not base_id:
            continue

        chebis = met.get('chebi', [])

        if chebis:
            if base_id not in mapping:
                mapping[base_id] = chebis[0]
        else:
            mnx_ids = met.get('metanetx', [])

            if mnx_ids and base_id not in mnx_pending:
                mnx_pending[base_id] = mnx_ids[0]

    n_direct = len(mapping)

    # MetaNetX bridge for metabolites without direct ChEBI
    if mnx_pending:
        try:
            from pypath.inputs.metanetx import metanetx_metabolite_chebi
            mnx_to_chebi = metanetx_metabolite_chebi()
            bridged = 0

            for base_id, mnx_id in mnx_pending.items():
                chebi = mnx_to_chebi.get(mnx_id)

                if chebi and base_id not in mapping:
                    mapping[base_id] = chebi
                    bridged += 1

        except Exception:
            bridged = 0

    return mapping


@cache
def _entrez_to_uniprot_bigg() -> dict[str, str]:
    """
    Build an Entrez Gene ID → UniProt AC mapping via BiGG gene symbols.

    The BiGG JSON gene objects carry a ``name`` field containing the HGNC
    gene symbol (e.g. ``'SLC25A21'``).  This function strips the ``_ATN``
    isoform suffix from the gene ``id`` to recover the base Entrez ID, then
    maps the gene symbol to UniProt using pypath's BioMart-backed
    ``genesymbol → uniprot`` mapping.

    The genesymbol → UniProt mapping table is downloaded once by pypath and
    cached on disk; subsequent calls are fast dictionary lookups.  The final
    entrez → UniProt dict is cached in memory for the session.

    Returns:
        Dict mapping Entrez Gene ID strings (e.g. ``'89874'``) to UniProt
        ACs (e.g. ``'Q9UBT5'``).
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

        uniprot_set = mapping_mod.map_name(name, 'genesymbol', 'uniprot', ncbi_tax_id=9606)

        if uniprot_set:
            result[entrez] = next(iter(uniprot_set))

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
def _hmdb_synonyms_chebi() -> dict[str, str]:
    """
    HMDB compound name/synonym → ChEBI ID mapping (lowercase keys).

    Delegates to :func:`pypath.inputs.hmdb.metabolites.synonyms_chebi`,
    which parses the HMDB XML once and is disk-cached by pypath curl.
    """

    from pypath.inputs.hmdb.metabolites import synonyms_chebi

    return synonyms_chebi()


@cache
def _ramp_synonyms_chebi() -> dict[str, str]:
    """
    RaMP synonym → ChEBI ID mapping (lowercase keys).

    Delegates to :func:`pypath.inputs.ramp._mapping.ramp_synonyms_chebi`,
    which inverts the RaMP ChEBI → synonym table.  RaMP aggregates
    HMDB, ChEBI, KEGG, WikiPathways, and Reactome synonyms.
    """

    from pypath.inputs.ramp._mapping import ramp_synonyms_chebi

    return ramp_synonyms_chebi()


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


def _to_uniprot(target_id: str, id_type: str, organism: int) -> str | None:
    """
    Translate a protein identifier to a UniProt accession.

    Args:
        target_id: The protein identifier.
        id_type: The identifier type — one of ``'uniprot'``, ``'ensp'``,
            ``'genesymbol'``, ``'ensembl'``, ``'entrez'``, or
            ``'reaction_id'`` (orphan pseudo-enzyme, passed through
            unchanged).
        organism: NCBI taxonomy ID.

    Returns:
        UniProt accession string, or ``None`` if translation is not possible.
    """

    import pypath.utils.mapping as mapping_mod

    if id_type == 'uniprot':
        # Already a UniProt AC — pass through as-is.
        return target_id

    elif id_type == 'reaction_id':
        # Orphan reaction pseudo-enzyme: pass through the reaction ID unchanged.
        return target_id

    elif id_type == 'entrez':
        # Recon3D Entrez Gene IDs (_ATN suffixes already stripped by recon3d.py).
        # Try BiGG-embedded gene symbol → UniProt first; fall back to
        # pypath BioMart ncbigene → uniprot mapping.
        uniprot = _entrez_to_uniprot_bigg().get(target_id)

        if uniprot:
            return uniprot

        result = mapping_mod.map_name(
            target_id,
            'ncbigene',
            'uniprot',
            ncbi_tax_id=organism,
        )

    elif id_type == 'ensp':
        # ENSP → UniProt (direct one-hop mapping via pypath BioMart)
        result = mapping_mod.map_name(
            target_id,
            'ensp',
            'uniprot',
            ncbi_tax_id=organism,
        )

    elif id_type == 'ensembl':
        # ENSG → UniProt (GEM enzyme IDs; translated for cross-network integration)
        result = mapping_mod.map_name(
            target_id,
            'ensg',
            'uniprot',
            ncbi_tax_id=organism,
        )

    elif id_type == 'genesymbol':
        result = mapping_mod.map_name(
            target_id,
            'genesymbol',
            'uniprot',
            ncbi_tax_id=organism,
        )

    else:
        _log.debug('Unknown protein id_type %r, cannot translate to UniProt.', id_type)
        return None

    if not result:
        return None

    if len(result) > 1:
        _log.debug(
            'Multiple UniProt ACs for %s (%s): %s — using first.',
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
        # Four-step fallback chain (bulk sources first, per-name HTTP last):
        #
        # 1. HMDB  — bulk XML, single download; covers common metabolomics
        #            compounds with direct name→ChEBI mapping.
        # 2. RaMP  — bulk SQLite; aggregates HMDB + ChEBI + KEGG + others;
        #            broader synonym coverage than HMDB alone.
        # 3. PubChem REST — per-name HTTP (disk-cached by pypath curl);
        #            covers drugs and synthetic substrates absent from HMDB/RaMP.
        # 4. PubChem CID → HMDB → ChEBI bridge — handles the case where
        #            PubChem knows the CID but UniChem lacks a ChEBI mapping.

        result: dict[str, str | None] = {uid: None for uid in unique_ids}

        # Step 1: HMDB bulk name lookup
        hmdb_map = _hmdb_synonyms_chebi()

        for uid in unique_ids:
            chebi = hmdb_map.get(uid.lower())
            if chebi:
                result[uid] = chebi

        # Step 2: RaMP synonym lookup for still-unresolved names
        unresolved = [uid for uid, v in result.items() if v is None]

        if unresolved:
            ramp_map = _ramp_synonyms_chebi()

            for uid in unresolved:
                chebi = ramp_map.get(uid.lower())
                if chebi:
                    result[uid] = chebi

        # Step 3: PubChem REST → ChEBI via UniChem
        unresolved = [uid for uid, v in result.items() if v is None]

        if unresolved:
            from pypath.inputs.pubchem import pubchem_names_cids

            name_to_cids = pubchem_names_cids(unresolved)
            pubchem_chebi = _pubchem_to_chebi()

            for uid, cids in name_to_cids.items():
                if cids and result[uid] is None:
                    result[uid] = pubchem_chebi.get(next(iter(cids)))

        # Step 4: PubChem CID → HMDB → ChEBI bridge
        # For names where PubChem returned a CID but UniChem has no ChEBI entry.
        unresolved = [uid for uid, v in result.items() if v is None]

        if unresolved:
            pub_hmdb = _pubchem_to_hmdb()
            hmdb_chebi = _hmdb_to_chebi()

            for uid in unresolved:
                cids = name_to_cids.get(uid, set())
                if cids:
                    hmdb_id = pub_hmdb.get(next(iter(cids)))
                    if hmdb_id:
                        result[uid] = hmdb_chebi.get(_normalise_hmdb(hmdb_id))

        return result

    _log.debug('Unknown metabolite id_type %r, cannot translate to ChEBI.', id_type)
    return {uid: None for uid in unique_ids}


def _build_protein_mapping(
    id_type: str,
    ids: pd.Series,
    organism: int,
) -> dict[str, str | None]:
    """
    Build a {source_id → UniProt} mapping for a group of protein IDs.

    Called once per ``(id_type, 'protein')`` group.  Performs BioMart
    lookups in bulk by iterating unique IDs and mapping through pypath's
    cached mapping tables.

    Args:
        id_type: Protein identifier type (e.g. ``'uniprot'``, ``'ensp'``).
        ids: Series of source IDs (same index as the group slice).
        organism: NCBI taxonomy ID.

    Returns:
        Dict mapping each unique source ID to its UniProt AC (or ``None``).
    """
    import pypath.utils.mapping as mapping_mod

    unique_ids = ids.unique()

    if id_type == 'uniprot':
        return {uid: uid for uid in unique_ids}

    if id_type == 'reaction_id':
        return {uid: uid for uid in unique_ids}

    try:
        from tqdm import tqdm as _tqdm
        _has_tqdm = True
    except ImportError:
        _has_tqdm = False

    def _progress(iterable, desc):
        if _has_tqdm:
            return _tqdm(iterable, desc = desc, unit = 'id')
        return iterable

    if id_type == 'entrez':
        bigg_map = _entrez_to_uniprot_bigg()
        result: dict[str, str | None] = {}
        for uid in _progress(unique_ids, 'entrez → uniprot'):
            uniprot = bigg_map.get(uid)
            if not uniprot:
                res = mapping_mod.map_name(uid, 'ncbigene', 'uniprot',
                                           ncbi_tax_id=organism)
                uniprot = next(iter(res)) if res else None
            result[uid] = uniprot
        return result

    if id_type == 'ensp':
        result = {}
        for uid in _progress(unique_ids, 'ensp → uniprot'):
            res = mapping_mod.map_name(uid, 'ensp', 'uniprot',
                                       ncbi_tax_id=organism)
            result[uid] = next(iter(res)) if res else None
        return result

    if id_type == 'ensembl':
        result = {}
        for uid in _progress(unique_ids, 'ensembl → uniprot'):
            res = mapping_mod.map_name(uid, 'ensg', 'uniprot',
                                       ncbi_tax_id=organism)
            result[uid] = next(iter(res)) if res else None
        return result

    if id_type == 'genesymbol':
        result = {}
        for uid in _progress(unique_ids, 'genesymbol → uniprot'):
            res = mapping_mod.map_name(uid, 'genesymbol', 'uniprot',
                                       ncbi_tax_id=organism)
            result[uid] = next(iter(res)) if res else None
        return result

    _log.debug('Unknown protein id_type %r, cannot translate to UniProt.', id_type)
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
    groups = df.groupby([id_type_col, entity_type_col]).groups

    try:
        from tqdm import tqdm
        groups = tqdm(
            groups.items(),
            total = len(groups),
            desc = f'translating {col}',
            unit = 'group',
        )
    except ImportError:
        groups = groups.items()

    for (id_type, entity_type), idx in groups:
        ids = df.loc[idx, col]

        if entity_type == 'small_molecule':
            mapping = _build_metab_mapping(id_type, ids, df.loc[idx, 'resource'])
        elif entity_type == 'protein':
            mapping = _build_protein_mapping(id_type, ids, organism)
        else:
            continue

        try:
            groups.set_description(
                f'translating {col}: {id_type} → '
                f'{"chebi" if entity_type == "small_molecule" else "uniprot"}'
                f' ({len(idx)} rows)'
            )
        except AttributeError:
            pass

        df.loc[idx, col] = ids.map(mapping)


def translate_pkn(df: pd.DataFrame, organism: int = 9606) -> pd.DataFrame:
    """
    Translate COSMOS PKN IDs to unified types.

    Metabolite source IDs are translated to ChEBI.  Protein target IDs are
    translated to UniProt accessions.  Rows where either translation
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
        DataFrame with ``source`` as ChEBI, ``target`` as UniProt, and updated
        ``id_type_a`` / ``id_type_b`` columns.  Index is reset.
    """
    df = df.copy()

    # Translation is direction-aware: source/target roles vary by resource.
    # GEM produces enzyme→metabolite edges (protein as source, metabolite as
    # target) in addition to the more common metabolite→protein direction.
    _translate_column(df, 'source', 'id_type_a', 'source_type', organism)
    _translate_column(df, 'target', 'id_type_b', 'target_type', organism)

    n_total = len(df)
    n_failed_source = df['source'].isna().sum()
    n_failed_target = df['target'].isna().sum()

    # Per-group totals — resource/id_type columns are unchanged by translation.
    total_by_source = df.groupby(['resource', 'id_type_a']).size()
    total_by_target = df.groupby(['resource', 'id_type_b']).size()

    if n_failed_source:
        pct_total = n_failed_source / n_total * 100 if n_total else 0
        _log.warning(
            'Dropped %d/%d (%.0f%%) rows: source ID could not be translated.',
            n_failed_source,
            n_total,
            pct_total,
        )
        breakdown = (
            df[df['source'].isna()]
            .groupby(['resource', 'id_type_a'])
            .size()
        )
        for (resource, id_type), count in breakdown.items():
            total = total_by_source.get((resource, id_type), count)
            pct = count / total * 100 if total else 0
            _log.warning(
                '  source: %d/%d (%.0f%%) from %s (id_type=%s)',
                count,
                total,
                pct,
                resource,
                id_type,
            )

    if n_failed_target:
        pct_total = n_failed_target / n_total * 100 if n_total else 0
        _log.warning(
            'Dropped %d/%d (%.0f%%) rows: target ID could not be translated.',
            n_failed_target,
            n_total,
            pct_total,
        )
        breakdown = (
            df[df['target'].isna()]
            .groupby(['resource', 'id_type_b'])
            .size()
        )
        for (resource, id_type), count in breakdown.items():
            total = total_by_target.get((resource, id_type), count)
            pct = count / total * 100 if total else 0
            _log.warning(
                '  target: %d/%d (%.0f%%) from %s (id_type=%s)',
                count,
                total,
                pct,
                resource,
                id_type,
            )

    df = df.dropna(subset=['source', 'target']).copy()

    _type_to_id = {'small_molecule': 'chebi', 'protein': 'uniprot'}

    # Preserve 'reaction_id' entries (orphan pseudo-enzymes) — only update
    # rows where the id_type is a translatable type.
    mask_a = df['id_type_a'] != 'reaction_id'
    mask_b = df['id_type_b'] != 'reaction_id'

    df.loc[mask_a, 'id_type_a'] = df.loc[mask_a, 'source_type'].map(_type_to_id)
    df.loc[mask_b, 'id_type_b'] = df.loc[mask_b, 'target_type'].map(_type_to_id)

    return df.reset_index(drop=True)
