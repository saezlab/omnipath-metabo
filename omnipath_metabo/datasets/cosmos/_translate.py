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
    - ``'synonym'``: name→ChEBI via HMDB bulk XML, RaMP bulk SQLite, then
      PubChem REST (filtered to plausible chemical names) (BRENDA)
    - ``'metatlas'``: GEM metabolites.tsv ``metsNoComp`` → ``metChEBIID``
      mapping; loaded once per GEM name from the ``resource`` column.
    - ``'bigg'``: Recon3D BiGG base metabolite ID → ChEBI, built from the
      annotation field of the BiGG JSON (``recon3d_metabolites()``).
    - ``'hmdb'``: UniChem HMDB→ChEBI mapping.  Old 5-digit HMDB IDs
      (e.g. ``HMDB00001``) are silently normalised to the current 7-digit
      format (e.g. ``HMDB0000001``) before lookup.

Translation strategies by protein id_type:
    - ``'uniprot'``: identity passthrough (TCDB, SLC, MRCLinksDB, BRENDA)
    - ``'ensp'``: pypath ``ensp → uniprot`` via BioMart (STITCH); deprecated
      ENSPs not in pypath fall back to UniProt ID Mapping REST API
    - ``'genesymbol'``: pypath ``genesymbol → uniprot`` (BRENDA fallback)
    - ``'ensembl'``: pypath ``ensg → uniprot`` via BioMart (GEM enzyme IDs
      are ENSG; translated to UniProt for cross-network integration).
      Compound ``_``-joined IDs (enzyme complexes, e.g. ``ENSG1_ENSG2``)
      are split and each component tried independently.  Deprecated single
      ENSGs fall back to UniProt ID Mapping REST API (``Ensembl → UniProtKB``)
    - ``'entrez'``: Recon3D Entrez Gene ID → UniProt.  Tries BiGG-embedded
      gene symbol → UniProt first; falls back to pypath
      ``ncbigene → uniprot`` BioMart mapping.  ``_ATN`` isoform suffixes
      must already be stripped by the caller (``recon3d.py`` does this at
      parse time).
"""

from __future__ import annotations

__all__ = ['translate_pkn', '_to_hmdb', '_to_uniprot', '_lipidmaps_to_chebi', '_metatlas_to_chebi', '_ensp_to_uniprot_rest', '_ensg_to_uniprot_rest']

import json
import logging
import re
import time
from functools import cache

import pandas as pd
import requests

_log = logging.getLogger(__name__)

_UNICHEM_PUBCHEM = 'PubChem'
_UNICHEM_CHEBI = 'ChEBI'
_UNICHEM_HMDB = 'HMDB'

_PUBCHEM_NAME_URL = (
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON'
)
_PUBCHEM_CID_XREFS_URL = (
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/xrefs/RegistryID/JSON'
)
# BiGG requires HTTP not HTTPS
_BIGG_METABOLITE_URL = 'http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'
# Maximum number of unique CIDs for which the per-CID PubChem REST fallback
# is attempted.  Above this threshold the residual set is likely irreducible
# (drugs/xenobiotics without ChEBI) and the REST calls would waste build time.
# MRCLinksDB leaves ≤ 4 CIDs unresolved after UniChem+RaMP — well within cap.
# STITCH leaves ~4,797 — skipped.
_PUBCHEM_REST_MAX_CIDS = 50

# UniProt ID Mapping REST API endpoints (batch ENSP → UniProtKB)
_UNIPROT_IDMAP_RUN = 'https://rest.uniprot.org/idmapping/run'
_UNIPROT_IDMAP_STATUS = 'https://rest.uniprot.org/idmapping/status/{}'
_UNIPROT_IDMAP_RESULTS = 'https://rest.uniprot.org/idmapping/results/{}?size=500'
# How long to wait between polls when the job is still RUNNING
_UNIPROT_IDMAP_POLL_INTERVAL = 2.0  # seconds
_UNIPROT_IDMAP_MAX_WAIT = 120.0      # seconds

# Patterns that indicate an experimental description rather than a chemical name.
_RE_NOT_CHEMICAL = re.compile(
    r'%'                        # percentages (e.g. "10.2% residual...")
    r'|\d+\s*(?:n|µ|u|m)M\b'  # concentrations (100 nM, 50 µM, 1 mM)
    r'|\b(?:activity|cells?|incubat|treatment|inhibitor|assay|compound'
    r'|concentration|after|protein|pathway|signaling|signalling)\b',
    re.IGNORECASE,
)
_MAX_CHEMICAL_NAME_LEN = 100


def _looks_like_chemical_name(name: str) -> bool:
    """
    Return True if *name* looks like a chemical name rather than an
    experimental description.

    Filters out strings containing percentages, concentration units, or
    assay-specific vocabulary that BRENDA includes as "synonyms" but are
    not resolvable compound names.
    """

    return len(name) <= _MAX_CHEMICAL_NAME_LEN and not _RE_NOT_CHEMICAL.search(name)


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

    _log.info('[COSMOS] Loading PubChem→ChEBI via UniChem (bulk)...')
    raw = unichem_mapping(_UNICHEM_PUBCHEM, _UNICHEM_CHEBI)
    result = {cid: next(iter(chebis)) for cid, chebis in raw.items() if chebis}
    _log.info('[COSMOS] PubChem→ChEBI loaded: %d entries', len(result))
    return result


@cache
def _pubchem_to_chebi_ramp() -> dict[str, str]:
    """
    Build a PubChem CID → ChEBI ID mapping via RaMP.

    Supplements :func:`_pubchem_to_chebi` (UniChem) by querying the RaMP
    SQLite database.  RaMP aggregates cross-references from HMDB, ChEBI,
    KEGG, WikiPathways, and Reactome, covering a different compound set
    than UniChem.  Downloaded once and cached for the session.

    Returns an empty dict if RaMP is unavailable.

    Returns:
        Dict mapping PubChem CID strings to ChEBI ID strings
        (e.g. ``{'5793': 'CHEBI:17234'}``).
    """

    try:

        from pypath.inputs.ramp._mapping import ramp_mapping

        _log.info('[COSMOS] Loading PubChem→ChEBI via RaMP...')
        # curies=False (default) strips prefixes: 'pubchem:5793' → '5793',
        # 'chebi:17234' → '17234'.  Re-add CHEBI: prefix to match UniChem format.
        raw = ramp_mapping('pubchem', 'chebi')
        result = {
            cid: 'CHEBI:' + next(iter(chebis))
            for cid, chebis in raw.items()
            if chebis
        }
        _log.info('[COSMOS] PubChem→ChEBI (RaMP): %d entries', len(result))
        return result

    except Exception as e:

        _log.warning(
            '[COSMOS] RaMP PubChem→ChEBI unavailable (%s), skipping',
            type(e).__name__,
        )
        return {}


@cache
def _pubchempy_cids_to_chebi(cids: tuple[str, ...]) -> dict[str, str | None]:
    """
    Look up ChEBI IDs for PubChem CIDs via the PubChem PUG REST API.

    Last-resort fallback for CIDs absent from both UniChem and RaMP.
    Queries the ``xrefs/RegistryID`` endpoint per CID and filters for
    ``CHEBI:``-prefixed entries.  Rate-limited to ≤ 5 requests/second per
    PubChem policy.

    Results are cached for the session; the argument must be a ``tuple``
    to satisfy ``@cache`` hashability.

    Args:
        cids: Tuple of PubChem CID strings.

    Returns:
        Dict mapping CID strings to ChEBI IDs (or ``None`` if not found).
    """

    import json
    import time
    import urllib.request

    result: dict[str, str | None] = {}

    for cid in cids:

        url = _PUBCHEM_CID_XREFS_URL.format(cid)

        try:

            with urllib.request.urlopen(url, timeout = 10) as resp:
                data = json.loads(resp.read())

            chebi = None

            for entry in data.get('InformationList', {}).get('Information', []):
                for rid in entry.get('RegistryID', []):
                    if rid.startswith('CHEBI:'):
                        chebi = rid
                        break
                if chebi:
                    break

            result[cid] = chebi

        except Exception:
            result[cid] = None

        time.sleep(0.2)  # ≤ 5 req/s per PubChem policy

    n_found = sum(1 for v in result.values() if v is not None)
    _log.info(
        '[COSMOS] PubChem CID→ChEBI (REST): %d/%d CIDs resolved',
        n_found,
        len(cids),
    )
    return result


def _normalise_hmdb(hmdb_id: str) -> str:

    from pypath.utils.metabo import normalise_hmdb

    return normalise_hmdb(hmdb_id)


@cache
def _lipidmaps_to_chebi() -> dict[str, str]:
    """
    Build a LipidMaps ID → ChEBI ID mapping via UniChem.

    Uses the UniChem LipidMaps→ChEBI bulk mapping.  When a LipidMaps ID
    maps to multiple ChEBI IDs the first one is used.  Downloaded once per
    session and cached in memory.

    Returns:
        Dict mapping LipidMaps IDs (e.g. ``'LMFA01010001'``) to ChEBI ID
        strings (e.g. ``'CHEBI:15756'``).
    """

    from pypath.inputs.unichem import unichem_mapping

    _UNICHEM_LIPIDMAPS = 'LIPID MAPS\u00ae'
    _log.info('[COSMOS] Loading LipidMaps→ChEBI via UniChem (bulk)...')
    raw = unichem_mapping(_UNICHEM_LIPIDMAPS, _UNICHEM_CHEBI)
    result = {lm_id: next(iter(chebis)) for lm_id, chebis in raw.items() if chebis}
    _log.info('[COSMOS] LipidMaps→ChEBI loaded: %d entries', len(result))
    return result


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

    _log.info('[COSMOS] Loading HMDB→ChEBI via UniChem (bulk)...')
    raw = unichem_mapping(_UNICHEM_HMDB, _UNICHEM_CHEBI)
    result = {
        _normalise_hmdb(hmdb_id): next(iter(chebis))
        for hmdb_id, chebis in raw.items()
        if chebis
    }
    _log.info('[COSMOS] HMDB→ChEBI loaded: %d entries', len(result))
    return result


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
            return _pubchem_to_chebi().get(next(iter(cids)))

    return None


@cache
def _bigg_to_chebi() -> dict[str, str]:
    """
    Build a BiGG base metabolite ID → ChEBI ID mapping from Recon3D.

    Uses ChEBI cross-references embedded in the Recon3D BiGG JSON
    (``recon3d_metabolites()``).  Metabolites without a direct ChEBI
    cross-reference in the JSON have no ChEBI annotation in any public
    database (BiGG API and MetaNetX bridge both return nothing for them) —
    the ~31% drop rate is irreducible.

    Downloaded once and cached for the session.

    Returns:
        Dict mapping BiGG base IDs (e.g. ``'atp'``) to ChEBI IDs
        (e.g. ``'CHEBI:30616'``).
    """

    from pypath.inputs.recon3d._gem import recon3d_metabolites

    _log.info('[COSMOS] Building BiGG→ChEBI from Recon3D JSON...')
    mapping: dict[str, str] = {}

    for met in recon3d_metabolites():
        base_id = met.get('base_id', '')
        chebis = met.get('chebi', [])

        if base_id and chebis and base_id not in mapping:
            mapping[base_id] = chebis[0]

    _log.info('[COSMOS] BiGG→ChEBI: %d entries', len(mapping))
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

    _log.info('[COSMOS] Building Entrez→UniProt via BiGG gene symbols...')
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

    _log.info('[COSMOS] Entrez→UniProt (BiGG): %d entries', len(result))
    return result


@cache
def _metatlas_to_chebi(gem: str) -> dict[str, str]:
    """
    Build a MetAtlas base metabolite ID → ChEBI mapping for a given GEM.

    Reads the GEM's ``metabolites.tsv`` (via
    :func:`pypath.inputs.metatlas.metatlas_gem_metabolites`) and builds
    a ``metsNoComp`` (base ID without compartment) → ChEBI mapping using
    a four-step fallback chain for rows where ``metChEBIID`` is absent:

    1. ``metChEBIID``   — direct annotation in the GEM (authoritative).
    2. ``metMetaNetXID`` → :func:`_metanetx_to_chebi` — covers ~25% of gap;
       may be semicolon-separated (all values tried).
    3. ``metLipidMapsID`` → :func:`_lipidmaps_to_chebi` — good for lipids
       (~83% hit rate on LM IDs present).
    4. ``metPubChemID``  → :func:`_pubchem_to_chebi` /
       :func:`_pubchem_to_chebi_ramp` — supplementary coverage.
    5. ``metHMDBID``     → :func:`_hmdb_to_chebi` — minor additional hits.

    Downloaded once per GEM name and cached for the session.

    Args:
        gem: GEM name (e.g. ``'Human-GEM'``).

    Returns:
        Dict mapping base MetAtlas IDs (e.g. ``'MAM00001'``) to ChEBI IDs
        (e.g. ``'CHEBI:15389'``).
    """

    from pypath.inputs.metatlas._gem import metatlas_gem_metabolites
    from pypath.inputs.metanetx import metanetx_metabolite_chebi

    _log.info('[COSMOS] Loading MetAtlas→ChEBI mapping for %s...', gem)
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
        '[COSMOS] MetAtlas→ChEBI (%s): %d direct; %d missing, trying fallbacks...',
        gem, n_direct, len(missing),
    )

    # Pre-load all bulk maps once (cached after first call)
    mnx_map = metanetx_metabolite_chebi()
    lm_map = _lipidmaps_to_chebi()
    pc_map = _pubchem_to_chebi()
    pc_ramp = _pubchem_to_chebi_ramp()
    hmdb_map = _hmdb_to_chebi()

    n_mnx = n_lm = n_pc = n_hmdb = 0

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
        if lm_id:
            chebi = lm_map.get(lm_id)
            if chebi:
                mapping[base_id] = chebi
                n_lm += 1
                continue

        # Step 4: PubChem (UniChem then RaMP)
        pc_id = row.get('metPubChemID', '')
        if pc_id:
            chebi = pc_map.get(pc_id) or pc_ramp.get(pc_id)
            if chebi:
                mapping[base_id] = chebi
                n_pc += 1
                continue

        # Step 5: HMDB
        hmdb_id = row.get('metHMDBID', '')
        if hmdb_id:
            chebi = hmdb_map.get(_normalise_hmdb(hmdb_id))
            if chebi:
                mapping[base_id] = chebi
                n_hmdb += 1

    _log.info(
        '[COSMOS] MetAtlas→ChEBI (%s): %d total '
        '(+%d MetaNetX, +%d LipidMaps, +%d PubChem, +%d HMDB fallback)',
        gem, len(mapping), n_mnx, n_lm, n_pc, n_hmdb,
    )
    return mapping


@cache
def _hmdb_synonyms_chebi() -> dict[str, str]:
    """
    HMDB compound name/synonym → ChEBI ID mapping (lowercase keys).

    Delegates to :func:`pypath.inputs.hmdb.metabolites.synonyms_chebi`,
    which parses the HMDB XML once and is disk-cached by pypath curl.
    Returns an empty dict if the download fails (e.g. Cloudflare block).
    """

    try:

        from pypath.inputs.hmdb.metabolites import synonyms_chebi

        result = synonyms_chebi()
        _log.info(f'[COSMOS] HMDB synonym→ChEBI map loaded: {len(result):,} entries')
        return result

    except Exception as e:

        _log.warning(f'[COSMOS] HMDB synonym→ChEBI unavailable ({type(e).__name__}), skipping')
        return {}


@cache
def _ramp_synonyms_chebi() -> dict[str, str]:
    """
    RaMP synonym → ChEBI ID mapping (lowercase keys).

    Delegates to :func:`pypath.inputs.ramp._mapping.ramp_synonyms_chebi`,
    which inverts the RaMP ChEBI → synonym table.  RaMP aggregates
    HMDB, ChEBI, KEGG, WikiPathways, and Reactome synonyms.
    Returns an empty dict if the download fails.
    """

    try:

        from pypath.inputs.ramp._mapping import ramp_synonyms_chebi

        result = ramp_synonyms_chebi()
        _log.info(f'[COSMOS] RaMP synonym→ChEBI map loaded: {len(result):,} entries')
        return result

    except Exception as e:

        _log.warning(f'[COSMOS] RaMP synonym→ChEBI unavailable ({type(e).__name__}), skipping')
        return {}


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
        cid = str(source_id)
        return (
            _pubchem_to_chebi().get(cid)
            or _pubchem_to_chebi_ramp().get(cid)
        )

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
        if not result:
            # Fallback for deprecated ENSPs via UniProt ID Mapping REST API
            rest_map = _ensp_to_uniprot_rest([target_id])
            uniprot = rest_map.get(target_id)
            if uniprot:
                return uniprot

    elif id_type == 'ensembl':
        # ENSG → UniProt (GEM enzyme IDs; translated for cross-network integration)
        result = mapping_mod.map_name(
            target_id,
            'ensg',
            'uniprot',
            ncbi_tax_id=organism,
        )
        if not result and '_' in target_id:
            # Compound enzyme-complex ID: try each component
            for part in target_id.split('_'):
                result = mapping_mod.map_name(part, 'ensg', 'uniprot',
                                              ncbi_tax_id=organism)
                if result:
                    break
        if not result and '_' not in target_id:
            # Single deprecated ENSG: REST fallback
            rest_map = _ensg_to_uniprot_rest([target_id])
            uniprot = rest_map.get(target_id)
            if uniprot:
                return uniprot

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
        unichem = _pubchem_to_chebi()
        ramp = _pubchem_to_chebi_ramp()
        result: dict[str, str | None] = {}
        residual: list[str] = []

        for uid in unique_ids:
            cid = str(uid)
            chebi = unichem.get(cid) or ramp.get(cid)
            result[uid] = chebi
            if chebi is None:
                residual.append(cid)

        if residual and len(residual) <= _PUBCHEM_REST_MAX_CIDS:
            _log.info(
                '[COSMOS] pubchem→ChEBI: %d/%d unique CIDs unresolved after '
                'UniChem+RaMP; trying PubChem REST for %d...',
                len(residual),
                len(unique_ids),
                len(residual),
            )
            rest_map = _pubchempy_cids_to_chebi(tuple(sorted(residual)))

            for uid in unique_ids:
                if result[uid] is None:
                    result[uid] = rest_map.get(str(uid))

        elif residual:
            _log.info(
                '[COSMOS] pubchem→ChEBI: %d/%d unique CIDs unresolved after '
                'UniChem+RaMP; skipping PubChem REST (residual %d > cap %d).',
                len(residual),
                len(unique_ids),
                len(residual),
                _PUBCHEM_REST_MAX_CIDS,
            )

        n_resolved = sum(1 for v in result.values() if v is not None)
        _log.info(
            '[COSMOS] pubchem→ChEBI: %d/%d unique CIDs resolved',
            n_resolved,
            len(unique_ids),
        )
        return result

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
        # Three-step fallback chain (bulk sources first, per-name HTTP last):
        #
        # 1. HMDB    — bulk XML, single download; covers common metabolomics
        #              compounds with direct name→ChEBI mapping.
        # 2. RaMP    — bulk SQLite; aggregates HMDB + ChEBI + KEGG + others;
        #              broader synonym coverage than HMDB alone.
        # 3. PubChem — per-name HTTP (disk-cached); only attempted for names
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
        _log.info('[COSMOS] synonym→ChEBI: %d/%d resolved by HMDB; trying RaMP for %d...', n_after_hmdb, len(unique_ids), len(unresolved))

        if unresolved:
            ramp_map = _ramp_synonyms_chebi()

            for uid in unresolved:
                chebi = ramp_map.get(uid.lower())
                if chebi:
                    result[uid] = chebi

        # Step 3: PubChem REST → ChEBI for plausible chemical names only
        candidates = [
            uid for uid, v in result.items()
            if v is None and _looks_like_chemical_name(uid)
        ]
        n_after_ramp = len(unique_ids) - len([v for v in result.values() if v is None])
        filtered_out = len([uid for uid, v in result.items() if v is None and not _looks_like_chemical_name(uid)])
        _log.info(
            '[COSMOS] synonym→ChEBI: %d/%d resolved after RaMP; %d candidates for PubChem (%d filtered as non-chemical)',
            n_after_ramp, len(unique_ids), len(candidates), filtered_out,
        )

        if candidates:
            from pypath.inputs.pubchem import pubchem_names_cids

            name_to_cids = pubchem_names_cids(candidates)
            pubchem_chebi = _pubchem_to_chebi()

            for uid, cids in name_to_cids.items():
                if cids and result[uid] is None:
                    result[uid] = pubchem_chebi.get(next(iter(cids)))

        n_final = len([v for v in result.values() if v is not None])
        _log.info('[COSMOS] synonym→ChEBI: %d/%d resolved total', n_final, len(unique_ids))
        return result

    _log.debug('Unknown metabolite id_type %r, cannot translate to ChEBI.', id_type)
    return {uid: None for uid in unique_ids}


def _uniprot_idmap_batch(ids: list[str], from_db: str) -> dict[str, str]:
    """
    Batch ID mapping via the UniProt ID Mapping REST API.

    Submits all *ids* in one job (``from_db → UniProtKB``), polls until
    complete, and returns a ``{input_id → UniProt_AC}`` dict for resolved
    entries only.  Any HTTP or job error returns ``{}`` silently.

    Args:
        ids: Input identifiers to map.
        from_db: UniProt source database name (e.g. ``'Ensembl_Protein'``,
            ``'Ensembl'``).

    Returns:
        Dict mapping each resolved input ID to its first UniProt accession.
    """

    if not ids:
        return {}

    result: dict[str, str] = {}

    try:
        resp = requests.post(
            _UNIPROT_IDMAP_RUN,
            data={'from': from_db, 'to': 'UniProtKB', 'ids': ','.join(ids)},
            timeout=30,
        )
        resp.raise_for_status()
        job_id = resp.json()['jobId']
    except Exception as exc:
        _log.warning('[COSMOS] UniProt idmap (%s): submit failed (%s)', from_db, exc)
        return {}

    deadline = time.monotonic() + _UNIPROT_IDMAP_MAX_WAIT

    while time.monotonic() < deadline:
        try:
            status_resp = requests.get(
                _UNIPROT_IDMAP_STATUS.format(job_id),
                timeout=30,
                allow_redirects=False,
            )
            if status_resp.status_code in (301, 302, 303):
                break
            status_data = status_resp.json()
            job_status = status_data.get('jobStatus', '')
            if job_status == 'FINISHED':
                break
            if job_status == 'ERROR':
                _log.warning('[COSMOS] UniProt idmap (%s): job error', from_db)
                return {}
        except Exception as exc:
            _log.warning('[COSMOS] UniProt idmap (%s): poll failed (%s)', from_db, exc)
            return {}
        time.sleep(_UNIPROT_IDMAP_POLL_INTERVAL)

    try:
        results_resp = requests.get(
            _UNIPROT_IDMAP_RESULTS.format(job_id),
            timeout=30,
        )
        results_resp.raise_for_status()
        data = results_resp.json()
    except Exception as exc:
        _log.warning('[COSMOS] UniProt idmap (%s): fetch failed (%s)', from_db, exc)
        return {}

    for entry in data.get('results', []):
        src = entry.get('from', '')
        to = entry.get('to', '')
        if src and to and src not in result:
            # 'to' is a string AC or a dict with 'primaryAccession'
            uniprot = to.get('primaryAccession', '') if isinstance(to, dict) else str(to)
            if uniprot:
                result[src] = uniprot

    _log.info(
        '[COSMOS] UniProt idmap (%s): %d/%d resolved',
        from_db, len(result), len(ids),
    )
    return result


def _ensp_to_uniprot_rest(ensp_ids: list[str]) -> dict[str, str]:
    """
    Resolve Ensembl protein IDs (ENSP) to UniProt accessions via the
    UniProt ID Mapping REST API (``Ensembl_Protein → UniProtKB``).

    Used as a fallback for ENSP IDs not found in pypath's BioMart table
    (typically deprecated/retired identifiers from older Ensembl releases).
    IDs that fail to resolve (truly retired with no UniProt record) are
    absent from the returned dict.

    Args:
        ensp_ids: List of ENSP identifiers to resolve.

    Returns:
        Dict mapping each resolved ENSP to its UniProt accession.
    """

    return _uniprot_idmap_batch(ensp_ids, 'Ensembl_Protein')


def _ensg_to_uniprot_rest(ensg_ids: list[str]) -> dict[str, str]:
    """
    Resolve Ensembl gene IDs (ENSG) to UniProt accessions via the
    UniProt ID Mapping REST API (``Ensembl → UniProtKB``).

    Used as a fallback for single ENSG IDs not found in pypath's BioMart
    table (typically deprecated identifiers from older Ensembl releases).
    Does **not** handle compound ``_``-joined IDs — split those before
    calling this function.

    Args:
        ensg_ids: List of single ENSG identifiers to resolve.

    Returns:
        Dict mapping each resolved ENSG to its UniProt accession.
    """

    return _uniprot_idmap_batch(ensg_ids, 'Ensembl')


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
        _log.info('[COSMOS] Translating %d Entrez IDs → UniProt...', len(unique_ids))
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

        # Fallback: deprecated ENSPs not in pypath's BioMart table
        missing = [uid for uid, v in result.items() if v is None]
        if missing:
            rest_map = _ensp_to_uniprot_rest(missing)
            for uid, uniprot in rest_map.items():
                result[uid] = uniprot

        return result

    if id_type == 'ensembl':
        _log.info('[COSMOS] Translating %d Ensembl IDs → UniProt...', len(unique_ids))
        result = {}
        for uid in _progress(unique_ids, 'ensembl → uniprot'):
            res = mapping_mod.map_name(uid, 'ensg', 'uniprot',
                                       ncbi_tax_id=organism)
            if res:
                result[uid] = next(iter(res))
            elif '_' in uid:
                # Compound enzyme-complex ID (e.g. ENSG1_ENSG2_ENSG3).
                # Try each component; take the first that resolves.
                for part in uid.split('_'):
                    res = mapping_mod.map_name(part, 'ensg', 'uniprot',
                                               ncbi_tax_id=organism)
                    if res:
                        result[uid] = next(iter(res))
                        break
                else:
                    result[uid] = None
            else:
                result[uid] = None

        # REST fallback for any still-missing single ENSGs
        missing_single = [
            uid for uid, v in result.items()
            if v is None and '_' not in uid
        ]
        if missing_single:
            rest_map = _ensg_to_uniprot_rest(missing_single)
            for uid, uniprot in rest_map.items():
                result[uid] = uniprot

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
        target_type = 'chebi' if entity_type == 'small_molecule' else 'uniprot'

        try:
            groups.set_description(
                f'translating {col}: {id_type} → {target_type} ({len(idx)} rows)'
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
