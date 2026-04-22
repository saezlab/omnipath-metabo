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
enzyme->metabolite edges where the protein is the *source*, so
``source_type``/``target_type`` columns are checked rather than assuming
a fixed metabolite-source / protein-target convention.

Protein translations (ENSP, ENSG, Entrez, gene symbol -> UniProt) and
metabolite ID translations (PubChem, HMDB, LipidMaps, KEGG, MetaNetX,
BiGG -> ChEBI) are delegated to omnipath-utils (direct DB access) with
automatic fallback to the omnipath-client HTTP API.  The adapter in
:mod:`._mapping` handles backend selection transparently.

Name-based synonym resolution for BRENDA compounds still uses HMDB XML,
RaMP SQLite, and PubChem REST as before (these are name lookups, not ID
translations).
"""

from __future__ import annotations

__all__ = ['translate_pkn', '_to_hmdb']

import logging
import re
from functools import cache

import pandas as pd

from ._mapping import mapping_table, mapping_translate

_log = logging.getLogger(__name__)


# Patterns that indicate an experimental description rather than a chemical name.
_RE_NOT_CHEMICAL = re.compile(
    r'%'                        # percentages (e.g. "10.2% residual...")
    r'|\d+\s*(?:n|µ|u|m)M\b'  # concentrations (100 nM, 50 uM, 1 mM)
    r'|\b(?:activity|cells?|incubat|treatment|inhibitor|assay|compound'
    r'|concentration|after|protein|pathway|signaling|signalling)\b',
    re.IGNORECASE,
)
_MAX_CHEMICAL_NAME_LEN = 100


# UniChem source names used by _to_hmdb coverage analysis functions.
_UNICHEM_PUBCHEM = 'PubChem'
_UNICHEM_CHEBI = 'ChEBI'
_UNICHEM_HMDB = 'HMDB'
_PUBCHEM_NAME_URL = (
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON'
)


def _looks_like_chemical_name(name: str) -> bool:
    """
    Return True if *name* looks like a chemical name rather than an
    experimental description.

    Filters out strings containing percentages, concentration units, or
    assay-specific vocabulary that BRENDA includes as "synonyms" but are
    not resolvable compound names.
    """

    return len(name) <= _MAX_CHEMICAL_NAME_LEN and not _RE_NOT_CHEMICAL.search(name)


# -- Metabolite translation tables via omnipath-utils adapter ----------------


_BATCH_CHUNK_SIZE = 100  # max IDs per HTTP request (keep well under 120s timeout)


def _chunked_translate(
    ids: list[str],
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int,
    raw: bool = False,
) -> dict[str, set[str]]:
    """Translate IDs in chunks to avoid HTTP timeout for large batches."""

    result: dict[str, set[str]] = {}

    for i in range(0, len(ids), _BATCH_CHUNK_SIZE):
        chunk = ids[i:i + _BATCH_CHUNK_SIZE]
        result.update(mapping_translate(chunk, id_type, target_id_type, ncbi_tax_id, raw=raw))

    return result


def _batch_to_chebi(
    ids: list[str],
    id_type: str,
) -> dict[str, str | None]:
    """Batch-translate metabolite IDs to ChEBI via omnipath-utils.

    Uses ``mapping_translate`` (batch query) rather than downloading the
    full translation table — this scales with the actual IDs needed
    instead of the full table size, and works in HTTP-only mode.

    Large batches are automatically chunked to avoid HTTP timeouts.
    Returns ``{id: chebi_or_None}`` for every input ID.  When multiple
    ChEBI hits exist, the lexicographically smallest (deterministic) is
    picked.
    """

    if not ids:
        return {}

    raw = _chunked_translate(list(set(ids)), id_type, 'chebi', 0)

    result: dict[str, str | None] = {}

    for uid in ids:
        hits = raw.get(uid)
        result[uid] = min(hits) if hits else None

    return result


# -- Name-based synonym resolution (BRENDA) ---------------------------------


def _name_to_chebi(name: str) -> str | None:
    """
    Translate a compound synonym/name to ChEBI.

    Tries HMDB and RaMP bulk maps first. Falls back to PubChem REST only
    for names that pass :func:`_looks_like_chemical_name`.

    Args:
        name: Compound name or synonym (e.g. ``'NAD+'``, ``'ATP'``).

    Returns:
        ChEBI ID string (e.g. ``'CHEBI:57540'``), or ``None`` if not found.
    """

    key = name.lower()
    chebi = _hmdb_synonyms_chebi().get(key) or _ramp_synonyms_chebi().get(key)

    if chebi:
        return chebi

    if _looks_like_chemical_name(name):
        from pypath.inputs.pubchem import pubchem_name_cids

        cids = pubchem_name_cids(name)
        if cids:
            cid = next(iter(cids))
            return _batch_to_chebi([str(cid)], 'pubchem').get(str(cid))

    return None


@cache
def _hmdb_synonyms_chebi() -> dict[str, str]:
    """
    HMDB compound name/synonym -> ChEBI ID mapping (lowercase keys).

    Delegates to :func:`pypath.inputs.hmdb.metabolites.synonyms_chebi`,
    which parses the HMDB XML once and is disk-cached by pypath curl.
    Returns an empty dict if the download fails (e.g. Cloudflare block).
    """

    try:

        from pypath.inputs.hmdb.metabolites import synonyms_chebi

        result = synonyms_chebi()
        _log.info(f'[COSMOS] HMDB synonym->ChEBI map loaded: {len(result):,} entries')
        return result

    except Exception as e:

        _log.warning(f'[COSMOS] HMDB synonym->ChEBI unavailable ({type(e).__name__}), skipping')
        return {}


@cache
def _ramp_synonyms_chebi() -> dict[str, str]:
    """
    RaMP synonym -> ChEBI ID mapping (lowercase keys).

    Delegates to :func:`pypath.inputs.ramp._mapping.ramp_synonyms_chebi`,
    which inverts the RaMP ChEBI -> synonym table.  RaMP aggregates
    HMDB, ChEBI, KEGG, WikiPathways, and Reactome synonyms.
    Returns an empty dict if the download fails.
    """

    try:

        from pypath.inputs.ramp._mapping import ramp_synonyms_chebi

        result = ramp_synonyms_chebi()
        _log.info(f'[COSMOS] RaMP synonym->ChEBI map loaded: {len(result):,} entries')
        return result

    except Exception as e:

        _log.warning(f'[COSMOS] RaMP synonym->ChEBI unavailable ({type(e).__name__}), skipping')
        return {}


# -- Resource-specific protein mapping (BiGG) --------------------------------


@cache
def _entrez_to_uniprot_bigg() -> dict[str, str]:
    """
    Build an Entrez Gene ID -> UniProt AC mapping via BiGG gene symbols.

    The BiGG JSON gene objects carry a ``name`` field containing the HGNC
    gene symbol (e.g. ``'SLC25A21'``).  This function strips the ``_ATN``
    isoform suffix from the gene ``id`` to recover the base Entrez ID, then
    maps the gene symbol to UniProt using the omnipath-utils adapter.

    Returns:
        Dict mapping Entrez Gene ID strings (e.g. ``'89874'``) to UniProt
        ACs (e.g. ``'Q9UBT5'``).
    """

    import re

    from pypath.inputs.recon3d._gem import recon3d_genes

    _log.info('[COSMOS] Building Entrez->UniProt via BiGG gene symbols...')
    result: dict[str, str] = {}

    # Collect all gene symbols to translate in one batch
    id_symbol_pairs: list[tuple[str, str]] = []

    for gene in recon3d_genes():
        raw_id = gene.get('id', '')
        name = gene.get('name', '')

        if not raw_id or not name:
            continue

        entrez = re.sub(r'_AT\d+$', '', raw_id)

        if entrez not in result:
            id_symbol_pairs.append((entrez, name))

    if not id_symbol_pairs:
        return result

    # Batch translate all gene symbols at once
    symbols = list({sym for _, sym in id_symbol_pairs})
    sym_to_uniprot = _chunked_translate(symbols, 'genesymbol', 'uniprot', 9606)

    for entrez, name in id_symbol_pairs:
        if entrez in result:
            continue
        targets = sym_to_uniprot.get(name)
        if targets:
            result[entrez] = next(iter(targets))

    _log.info('[COSMOS] Entrez->UniProt (BiGG): %d entries', len(result))
    return result


# -- MetAtlas (GEM) metabolite -> ChEBI --------------------------------------


@cache
def _metatlas_to_chebi(gem: str) -> dict[str, str]:
    """
    Build a MetAtlas base metabolite ID -> ChEBI mapping for a given GEM.

    Reads the GEM's ``metabolites.tsv`` and builds a ``metsNoComp``
    (base ID without compartment) -> ChEBI mapping using a six-step
    fallback chain for rows where ``metChEBIID`` is absent:

    1. ``metChEBIID``    -- direct annotation in the GEM (authoritative).
    2. ``metMetaNetXID`` -- MetaNetX->ChEBI via omnipath-utils DB.
    3. ``metLipidMapsID`` -- LipidMaps->ChEBI via omnipath-utils DB.
    4. ``metPubChemID``  -- PubChem->ChEBI via omnipath-utils DB.
    5. ``metKEGGID``     -- KEGG->ChEBI via omnipath-utils DB.
    6. ``metHMDBID``     -- HMDB->ChEBI via omnipath-utils DB.

    Downloaded once per GEM name and cached for the session.

    Args:
        gem: GEM name (e.g. ``'Human-GEM'``).

    Returns:
        Dict mapping base MetAtlas IDs (e.g. ``'MAM00001'``) to ChEBI IDs
        (e.g. ``'CHEBI:15389'``).
    """

    from pypath.inputs.metatlas._gem import metatlas_gem_metabolites

    _log.info('[COSMOS] Loading MetAtlas->ChEBI mapping for %s...', gem)
    rows = list(metatlas_gem_metabolites(gem=gem))

    # Step 1: direct ChEBI annotation
    mapping: dict[str, str] = {}
    missing: list[dict] = []

    for row in rows:
        base_id = row.get('metsNoComp', '')
        chebi = row.get('metChEBIID', '')
        if not base_id:
            continue
        if chebi:
            mapping[base_id] = chebi
        else:
            missing.append(row)

    n_direct = len(mapping)
    _log.info(
        '[COSMOS] MetAtlas->ChEBI (%s): %d direct; %d missing, trying fallbacks...',
        gem, n_direct, len(missing),
    )

    # Collect cross-reference IDs across all missing rows, then batch-translate
    # each external ID type in one round-trip via the adapter.  This scales
    # with the actual IDs needed rather than the size of the reference tables.
    mnx_ids: set[str] = set()
    lm_ids: set[str] = set()
    pc_ids: set[str] = set()
    kegg_ids: set[str] = set()
    hmdb_ids: set[str] = set()

    for row in missing:
        for part in (row.get('metMetaNetXID') or '').split(';'):
            part = part.strip()
            if part:
                mnx_ids.add(part)
        if row.get('metLipidMapsID'):
            lm_ids.add(row['metLipidMapsID'])
        if row.get('metPubChemID'):
            pc_ids.add(row['metPubChemID'])
        if row.get('metKEGGID'):
            kegg_ids.add(row['metKEGGID'])
        if row.get('metHMDBID'):
            hmdb_ids.add(row['metHMDBID'])

    mnx_map = _batch_to_chebi(list(mnx_ids), 'metanetx') if mnx_ids else {}
    lm_map = _batch_to_chebi(list(lm_ids), 'lipidmaps') if lm_ids else {}
    pc_map = _batch_to_chebi(list(pc_ids), 'pubchem') if pc_ids else {}
    kegg_map = _batch_to_chebi(list(kegg_ids), 'kegg') if kegg_ids else {}
    hmdb_map = _batch_to_chebi(list(hmdb_ids), 'hmdb') if hmdb_ids else {}

    n_mnx = n_lm = n_pc = n_kegg = n_hmdb = 0

    for row in missing:
        base_id = row['metsNoComp']
        if base_id in mapping:
            continue

        # Step 2: MetaNetX (may be semicolon-separated)
        mnx_raw = row.get('metMetaNetXID', '')
        if mnx_raw:
            for mnx_id in mnx_raw.split(';'):
                chebi = mnx_map.get(mnx_id.strip())
                if chebi:
                    mapping[base_id] = chebi
                    n_mnx += 1
                    break
            if base_id in mapping:
                continue

        # Step 3: LipidMaps
        lm_id = row.get('metLipidMapsID', '')
        if lm_id and lm_map.get(lm_id):
            mapping[base_id] = lm_map[lm_id]
            n_lm += 1
            continue

        # Step 4: PubChem
        pc_id = row.get('metPubChemID', '')
        if pc_id and pc_map.get(pc_id):
            mapping[base_id] = pc_map[pc_id]
            n_pc += 1
            continue

        # Step 5: KEGG compound
        kegg_id = row.get('metKEGGID', '')
        if kegg_id and kegg_map.get(kegg_id):
            mapping[base_id] = kegg_map[kegg_id]
            n_kegg += 1
            continue

        # Step 6: HMDB
        hmdb_id = row.get('metHMDBID', '')
        if hmdb_id and hmdb_map.get(hmdb_id):
            mapping[base_id] = hmdb_map[hmdb_id]
            n_hmdb += 1

    _log.info(
        '[COSMOS] MetAtlas->ChEBI (%s): %d total '
        '(+%d MetaNetX, +%d LipidMaps, +%d PubChem, +%d KEGG, +%d HMDB fallback)',
        gem, len(mapping), n_mnx, n_lm, n_pc, n_kegg, n_hmdb,
    )
    return mapping


# -- Coverage analysis helpers (_to_hmdb) ------------------------------------


@cache
def _pubchem_to_hmdb() -> dict[str, str]:
    """
    Build a PubChem CID -> HMDB ID mapping via UniChem.

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
    Build a ChEBI ID -> HMDB ID mapping via UniChem.

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
    Build a MetAtlas base metabolite ID -> HMDB mapping for a given GEM.

    Reads the GEM's ``metabolites.tsv`` and extracts the ``metsNoComp`` ->
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
    Build a BiGG base metabolite ID -> HMDB ID mapping from Recon3D.

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

    Intended for use in coverage comparison scripts -- not used by the
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
        return source_id

    if id_type == 'chebi':
        return _chebi_to_hmdb().get(source_id)

    if id_type == 'pubchem':
        return _pubchem_to_hmdb().get(str(source_id))

    if id_type == 'bigg':
        return _bigg_to_hmdb().get(source_id)

    if id_type == 'metatlas':
        return _metatlas_to_hmdb(gem).get(source_id) if gem else None

    if id_type == 'synonym':
        # name -> pubchem CID -> HMDB
        import json
        import urllib.parse
        import urllib.request

        url = _PUBCHEM_NAME_URL.format(urllib.parse.quote(source_id))
        try:
            with urllib.request.urlopen(url, timeout = 10) as resp:
                data = json.loads(resp.read())
            cids = data.get('IdentifierList', {}).get('CID', [])
            if not cids:
                return None
            return _pubchem_to_hmdb().get(str(cids[0]))
        except Exception:
            return None

    _log.debug('Unknown metabolite id_type %r, cannot translate to HMDB.', id_type)
    return None


# -- Bulk mapping builders ---------------------------------------------------


def _normalize_chebi(chebi_id: str | None) -> str | None:
    """Ensure a ChEBI ID has the ``CHEBI:`` prefix (e.g. ``'3647'`` → ``'CHEBI:3647'``)."""
    if chebi_id is None:
        return None
    s = str(chebi_id)
    return s if s.startswith('CHEBI:') else f'CHEBI:{s}'


def _build_metab_mapping(
    id_type: str,
    ids: pd.Series,
    resource_series: pd.Series,
) -> dict[str, str | None]:
    """
    Build a {source_id -> ChEBI} mapping for a group of metabolite IDs.

    Called once per ``(id_type, 'small_molecule')`` group.  Uses the cached
    bulk dictionaries from omnipath-utils wherever possible; falls back to
    per-name HTTP calls only for ``'synonym'`` (BRENDA) entries.

    All returned ChEBI IDs are normalised to have the ``CHEBI:`` prefix
    (some upstream sources omit it).

    Args:
        id_type: Metabolite identifier type (e.g. ``'pubchem'``, ``'metatlas'``).
        ids: Series of source IDs (same index as the group slice).
        resource_series: Series of resource names (same index).

    Returns:
        Dict mapping each unique source ID to its ChEBI ID (or ``None``).
    """
    unique_ids = ids.unique()

    if id_type == 'chebi':
        return {uid: _normalize_chebi(uid) for uid in unique_ids}

    if id_type in ('pubchem', 'bigg', 'hmdb'):
        result = _batch_to_chebi(
            [str(uid) for uid in unique_ids],
            id_type,
        )
        n_resolved = sum(1 for v in result.values() if v is not None)
        _log.info(
            '[COSMOS] %s->ChEBI: %d/%d unique IDs resolved',
            id_type,
            n_resolved,
            len(unique_ids),
        )
        return {k: _normalize_chebi(v) for k, v in result.items()}

    if id_type == 'metatlas':
        # Each row may come from a different GEM -- build per-GEM mappings.
        gem_names = _gem_name_series(resource_series)
        result: dict[str, str | None] = {}
        for uid in unique_ids:
            # Find the GEM name for this metabolite ID (use first occurrence).
            uid_mask = ids == uid
            gem = gem_names[uid_mask].iloc[0] if uid_mask.any() else ''
            result[uid] = _metatlas_to_chebi(gem).get(uid) if gem else None
        return {k: _normalize_chebi(v) for k, v in result.items()}

    if id_type == 'synonym':
        # Three-step fallback chain (bulk sources first, per-name HTTP last):
        #
        # 1. HMDB    -- bulk XML, single download; covers common metabolomics
        #              compounds with direct name->ChEBI mapping.
        # 2. RaMP    -- bulk SQLite; aggregates HMDB + ChEBI + KEGG + others;
        #              broader synonym coverage than HMDB alone.
        # 3. PubChem -- per-name HTTP (disk-cached); only attempted for names
        #              that pass _looks_like_chemical_name() to skip BRENDA's
        #              experimental-description pseudo-synonyms.

        result: dict[str, str | None] = {uid: None for uid in unique_ids}

        # Step 1: HMDB bulk name lookup
        hmdb_map = _hmdb_synonyms_chebi()

        for uid in unique_ids:
            chebi = hmdb_map.get(uid.lower())
            if chebi:
                result[uid] = chebi

        # Step 2: RaMP synonym lookup for still-unresolved names
        unresolved = [uid for uid, v in result.items() if v is None]
        n_after_hmdb = len(unique_ids) - len(unresolved)
        _log.info('[COSMOS] synonym->ChEBI: %d/%d resolved by HMDB; trying RaMP for %d...', n_after_hmdb, len(unique_ids), len(unresolved))

        if unresolved:
            ramp_map = _ramp_synonyms_chebi()

            for uid in unresolved:
                chebi = ramp_map.get(uid.lower())
                if chebi:
                    result[uid] = chebi

        # Step 3: PubChem REST -> ChEBI for plausible chemical names only
        candidates = [
            uid for uid, v in result.items()
            if v is None and _looks_like_chemical_name(uid)
        ]
        n_after_ramp = len(unique_ids) - len([v for v in result.values() if v is None])
        filtered_out = len([uid for uid, v in result.items() if v is None and not _looks_like_chemical_name(uid)])
        _log.info(
            '[COSMOS] synonym->ChEBI: %d/%d resolved after RaMP; %d candidates for PubChem (%d filtered as non-chemical)',
            n_after_ramp, len(unique_ids), len(candidates), filtered_out,
        )

        if candidates:
            from pypath.inputs.pubchem import pubchem_names_cids

            name_to_cids = pubchem_names_cids(candidates)

            # Collect all CIDs and batch-translate in one call
            all_cids: set[str] = set()
            for cids in name_to_cids.values():
                if cids:
                    all_cids.add(str(next(iter(cids))))

            cid_to_chebi = _batch_to_chebi(list(all_cids), 'pubchem')

            for uid, cids in name_to_cids.items():
                if cids and result[uid] is None:
                    cid = str(next(iter(cids)))
                    result[uid] = cid_to_chebi.get(cid)

        n_final = len([v for v in result.values() if v is not None])
        _log.info('[COSMOS] synonym->ChEBI: %d/%d resolved total', n_final, len(unique_ids))
        return {k: _normalize_chebi(v) for k, v in result.items()}

    _log.debug('Unknown metabolite id_type %r, cannot translate to ChEBI.', id_type)
    return {uid: None for uid in unique_ids}


def _build_protein_mapping(
    id_type: str,
    ids: pd.Series,
    organism: int,
) -> dict[str, frozenset | None]:
    """
    Build a {source_id -> frozenset[UniProt]} mapping for a group of protein IDs.

    Called once per ``(id_type, 'protein')`` group.  Uses omnipath-utils
    batch translation with automatic fallback for deprecated IDs.

    Args:
        id_type: Protein identifier type (e.g. ``'uniprot'``, ``'ensp'``).
        ids: Series of source IDs (same index as the group slice).
        organism: NCBI taxonomy ID.

    Returns:
        Dict mapping each unique source ID to a frozenset of UniProt ACs
        (or ``None``).
    """

    unique_ids = ids.unique()

    if id_type == 'uniprot':
        return {uid: frozenset({uid}) for uid in unique_ids}

    if id_type == 'reaction_id':
        return {uid: frozenset({uid}) for uid in unique_ids}

    if id_type == 'entrez':
        _log.info('[COSMOS] Translating %d Entrez IDs -> UniProt...', len(unique_ids))

        # Try BiGG-embedded gene symbol -> UniProt first (resource-specific)
        bigg_map = _entrez_to_uniprot_bigg()
        result: dict[str, frozenset | None] = {}
        remaining = []

        for uid in unique_ids:
            uniprot = bigg_map.get(uid)
            if uniprot:
                result[uid] = frozenset({uniprot})
            else:
                remaining.append(uid)
                result[uid] = None

        # Batch translate remaining via omnipath-utils
        if remaining:
            batch = _chunked_translate(remaining, 'entrez', 'uniprot', organism)
            for uid in remaining:
                targets = batch.get(uid)
                if targets:
                    result[uid] = frozenset(targets)

        return result

    if id_type == 'ensp':
        _log.info('[COSMOS] Translating %d ENSP IDs -> UniProt...', len(unique_ids))

        batch = _chunked_translate(list(unique_ids), 'ensp', 'uniprot', organism)
        return {
            uid: frozenset(targets) if (targets := batch.get(uid)) else None
            for uid in unique_ids
        }

    if id_type == 'ensembl':
        _log.info('[COSMOS] Translating %d Ensembl IDs -> UniProt...', len(unique_ids))

        # Separate simple and compound IDs
        simple_ids = []
        compound_ids = []
        for uid in unique_ids:
            if '_' in uid:
                compound_ids.append(uid)
            else:
                simple_ids.append(uid)

        # Batch translate all simple IDs + all parts of compound IDs
        all_parts = set(simple_ids)
        for uid in compound_ids:
            all_parts.update(uid.split('_'))

        batch = _chunked_translate(list(all_parts), 'ensg', 'uniprot', organism)

        result: dict[str, frozenset | None] = {}

        # Simple IDs: direct lookup
        for uid in simple_ids:
            targets = batch.get(uid)
            result[uid] = frozenset(targets) if targets else None

        # Compound IDs: union of all component results
        for uid in compound_ids:
            all_res: set[str] = set()
            for part in uid.split('_'):
                targets = batch.get(part)
                if targets:
                    all_res.update(targets)
            result[uid] = frozenset(all_res) if all_res else None

        return result

    if id_type == 'genesymbol':
        _log.info('[COSMOS] Translating %d gene symbols -> UniProt...', len(unique_ids))

        batch = _chunked_translate(list(unique_ids), 'genesymbol', 'uniprot', organism)
        return {
            uid: frozenset(targets) if (targets := batch.get(uid)) else None
            for uid in unique_ids
        }

    _log.debug('Unknown protein id_type %r, cannot translate to UniProt.', id_type)
    return {uid: None for uid in unique_ids}


# -- Orchestration -----------------------------------------------------------


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
        target_type = 'chebi' if entity_type == 'small_molecule' else 'uniprot'

        try:
            groups.set_description(
                f'translating {col}: {id_type} -> {target_type} ({len(idx)} rows)'
            )
        except AttributeError:
            pass

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

    # Arrow string dtype (inferred by pandas from Interaction namedtuples) rejects
    # frozenset values -- cast source/target to object before translation.
    for _col in ('source', 'target'):
        if df[_col].dtype != object:
            df[_col] = df[_col].astype(object)

    # Translation is direction-aware: source/target roles vary by resource.
    # GEM produces enzyme->metabolite edges (protein as source, metabolite as
    # target) in addition to the more common metabolite->protein direction.
    _translate_column(df, 'source', 'id_type_a', 'source_type', organism)
    _translate_column(df, 'target', 'id_type_b', 'target_type', organism)

    n_total = len(df)
    n_failed_source = df['source'].isna().sum()
    n_failed_target = df['target'].isna().sum()

    # Per-group totals -- resource/id_type columns are unchanged by translation.
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

    # Preserve 'reaction_id' entries (orphan pseudo-enzymes) -- only update
    # rows where the id_type is a translatable type.
    mask_a = df['id_type_a'] != 'reaction_id'
    mask_b = df['id_type_b'] != 'reaction_id'

    df.loc[mask_a, 'id_type_a'] = df.loc[mask_a, 'source_type'].map(_type_to_id)
    df.loc[mask_b, 'id_type_b'] = df.loc[mask_b, 'target_type'].map(_type_to_id)

    return df.reset_index(drop=True)
