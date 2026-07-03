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
CellPhoneDB v5.0 non-peptidic ligand and PPI interactions for COSMOS PKN (US6).

Data acquisition is delegated to the raw CSV layer of
:mod:`pypath.inputs_v2.cellphonedb` (``resource.interactions.raw()`` /
``resource.complexes.raw()``, plain dicts) -- the same pattern every other
inputs_v2-native COSMOS resource uses (see ``kegg.py``). A newer typed
``Entity`` layer exists upstream but is not yet pulled into this package's
pinned ``pypath`` version; see the spec's research.md R10 for the rationale.

``interaction_input.csv`` (2,911 rows total) splits into three slices via
the ``is_ppi`` column and the shape of ``partner_a``:

- **Non-peptidic ligand** (1,029 rows, ``is_ppi=False`` and ``partner_a``
  matches :data:`SYNTH_RE`): a synthetic label such as
  ``'Adenosine_byNT5E_and_SLC29A1'`` encodes the metabolite name
  (``'Adenosine'``) and the producing-cell machinery gene symbols
  (``['NT5E', 'SLC29A1']``). The machinery is stored in
  ``attrs['producing_machinery']`` for reference only -- no
  enzyme->metabolite edges are yielded (GEM and KEGG already cover that
  layer). ``partner_b`` (the receiving-cell protein or complex) is routed
  to ``'transport'`` or ``'ligand_receptor'`` via CellPhoneDB's own
  ``protein_input.csv`` curation (:func:`_build_protein_type_map`) --
  not the generic cross-database classifier other resources
  (``mrclinksdb.py``) reuse from ``stitch.py``, since CellPhoneDB curates
  its own receptor/transporter roles for exactly these proteins.
- **PPI** (1,878 rows, ``is_ppi=True``): protein->protein, with MOR derived
  from ``modulatory_effect``.
- **Out-of-scope** (4 rows, ``is_ppi=False`` but ``partner_a`` does NOT
  match :data:`SYNTH_RE` -- the ``hCGB1/2/3/7_complex -> P22888`` glycoprotein
  hormone rows): logged individually and skipped, never force-included in
  either slice (FR-028).

Both slices resolve ``partner_a``/``partner_b`` complex names via
``complex_input.csv``, expanding to one edge per constituent UniProt AC;
values that are neither a known complex nor a bare UniProt accession are
logged and skipped.
"""

from __future__ import annotations

__all__ = ['cellphonedb_nonpeptidic_interactions', 'cellphonedb_ppi_interactions']

import logging
import re
from collections.abc import Generator

from .._record import Interaction

_log = logging.getLogger(__name__)

# Matches CellPhoneDB's synthetic-metabolite partner_a label, e.g.
# 'Adenosine_byNT5E_and_SLC29A1' -> name 'Adenosine', machinery
# 'NT5E_and_SLC29A1'. Rows that don't match (e.g. 'hCGB1_complex') are not
# metabolite systems -- see FR-028.
SYNTH_RE = re.compile(r'^(.+?)_by([A-Za-z0-9_]+)$')

_UNIPROT_RE = re.compile(r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$')


# ---------------------------------------------------------------------------
# Data loading and shared parsing helpers
# ---------------------------------------------------------------------------

def _load_cellphonedb_data() -> tuple[list[dict], list[dict]]:
    """
    Load raw CellPhoneDB interaction and complex rows.

    Returns:
        ``(interaction_rows, complex_rows)`` -- plain dicts from
        ``interaction_input.csv`` and ``complex_input.csv`` respectively.
    """
    from pypath.inputs_v2 import cellphonedb

    resource = cellphonedb.resource
    interaction_rows = list(resource.interactions.raw())
    complex_rows = list(resource.complexes.raw())
    return interaction_rows, complex_rows


def _load_cellphonedb_protein_types() -> list[dict]:
    """
    Load raw CellPhoneDB ``protein_input.csv`` rows.

    Not yet wired into ``resource.proteins`` (only ``.interactions``/
    ``.complexes`` are registered `Dataset`s in the pinned pypath version),
    so this fetches and parses it directly via the same ``Download`` +
    CSV-parsing primitives the registered datasets use internally.
    """
    from pypath.inputs_v2 import cellphonedb

    opener = cellphonedb.download_proteins.open()
    return list(cellphonedb.iter_csv(opener))


def _build_protein_type_map(protein_rows: list[dict]) -> dict[str, str]:
    """
    Build a ``{uniprot: 'receptor' | 'transporter'}`` map from CellPhoneDB's
    own ``protein_input.csv`` curation (R2, revised 2026-07-03).

    Uses CellPhoneDB's own ``receptor`` column and ``'Transporter'`` tag --
    curated specifically for these proteins in this dataset -- rather than
    a generic cross-database classifier (OmniPath Intercell + TCDB + Guide
    to Pharmacology, as reused from ``stitch.py`` in the original R2
    decision). Proteins flagged as neither are absent from the map (caller
    defaults to ``'other'``).
    """
    result: dict[str, str] = {}
    for row in protein_rows:
        uniprot = row.get('uniprot', '')
        if not uniprot:
            continue
        if row.get('receptor') == 'TRUE':
            result[uniprot] = 'receptor'
        elif 'Transporter' in (row.get('tags') or ''):
            result[uniprot] = 'transporter'
    return result


def _build_complex_map(complex_rows: list[dict]) -> dict[str, list[str]]:
    """
    Build a ``{complex_name: [uniprot, ...]}`` map from complex rows.

    Synthetic metabolite-system rows (``complex_name`` contains ``'_by'``)
    are excluded -- these describe producing-cell machinery, not protein
    complexes.
    """
    result: dict[str, list[str]] = {}
    for row in complex_rows:
        name = row.get('complex_name', '')
        if not name or '_by' in name:
            continue
        members = [
            row[f'uniprot_{i}']
            for i in range(1, 6)
            if row.get(f'uniprot_{i}')
        ]
        if members:
            result[name] = members
    return result


def _parse_synthetic_metabolite(partner_a: str) -> tuple[str, list[str]] | None:
    """
    Extract the metabolite name and producing-cell machinery from a
    CellPhoneDB synthetic ``partner_a`` label.

    Returns:
        ``(metabolite_name, machinery_gene_symbols)``, or ``None`` if
        *partner_a* does not match the synthetic-metabolite pattern
        (FR-028 -- e.g. the ``hCGB*_complex`` rows).
    """
    match = SYNTH_RE.match(partner_a)
    if match is None:
        return None
    name = match.group(1)
    machinery = match.group(2).split('_and_')
    return name, machinery


def _resolve_entity(value: str, complex_map: dict[str, list[str]]) -> list[str]:
    """
    Resolve a ``partner_a``/``partner_b`` value to one or more UniProt ACs.

    - Known complex name -> its constituent UniProt members.
    - UniProt AC shape -> itself (single-element list).
    - Neither -> empty list (callers log and skip).
    """
    if value in complex_map:
        return complex_map[value]
    if _UNIPROT_RE.match(value):
        return [value]
    return []


def _mor_from_modulatory_effect(modulatory_effect: str) -> int:
    """``'Inhibitory'`` -> ``-1``; absent/empty/anything else -> ``+1``."""
    return -1 if modulatory_effect == 'Inhibitory' else 1


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def cellphonedb_nonpeptidic_interactions(
    organism: int = 9606,
    **_kwargs,
) -> Generator[Interaction, None, None]:
    """
    Yield CellPhoneDB v5.0 non-peptidic ligand interactions.

    Metabolite/small-molecule source (extracted from the synthetic
    ``partner_a`` label) to protein target (``partner_b``, UniProt AC or
    complex, expanded via ``complex_input.csv``). Producing-cell machinery
    genes are stored in ``attrs['producing_machinery']`` only -- never
    yielded as edges.

    Args:
        organism: NCBI taxonomy ID. CellPhoneDB is human-curated only;
            any organism other than 9606 yields nothing.

    Yields:
        :class:`~.._record.Interaction` records with
        ``source_type='small_molecule'``, ``id_type_a='name'``,
        ``interaction_type`` routed to ``'transport'`` or
        ``'ligand_receptor'`` via CellPhoneDB's own ``protein_input.csv``
        curation (:func:`_build_protein_type_map`, R2 revised 2026-07-03),
        ``resource='CellPhoneDB'``, ``mor=1``.
    """
    if organism != 9606:
        _log.info(
            '[COSMOS] CellPhoneDB non-peptidic: human-curated only, '
            'skipping organism %d.', organism,
        )
        return

    interaction_rows, complex_rows = _load_cellphonedb_data()
    complex_map = _build_complex_map(complex_rows)
    protein_types = _build_protein_type_map(_load_cellphonedb_protein_types())

    n_yielded = 0
    n_skipped_nonmatch = 0
    n_skipped_unresolved = 0

    for row in interaction_rows:
        if row.get('is_ppi') != 'False':
            continue

        partner_a = row.get('partner_a', '')
        parsed = _parse_synthetic_metabolite(partner_a)

        if parsed is None:
            n_skipped_nonmatch += 1
            _log.info(
                "[COSMOS] CellPhoneDB non-peptidic: partner_a '%s' does not match "
                'the synthetic-metabolite pattern, skipping (FR-028).', partner_a,
            )
            continue

        metabolite_name, machinery = parsed
        partner_b = row.get('partner_b', '')
        targets = _resolve_entity(partner_b, complex_map)

        if not targets:
            n_skipped_unresolved += 1
            _log.info(
                "[COSMOS] CellPhoneDB non-peptidic: partner_b '%s' is not a known "
                'complex or UniProt AC, skipping.', partner_b,
            )
            continue

        attrs = {'producing_machinery': machinery}

        for target in targets:
            ptype = protein_types.get(target, 'other')
            interaction_type = 'transport' if ptype == 'transporter' else 'ligand_receptor'

            n_yielded += 1
            yield Interaction(
                source=metabolite_name,
                target=target,
                source_type='small_molecule',
                target_type='protein',
                id_type_a='name',
                id_type_b='uniprot',
                interaction_type=interaction_type,
                resource='CellPhoneDB',
                mor=1,
                locations=(),
                attrs=attrs,
            )

    _log.info(
        '[COSMOS] CellPhoneDB non-peptidic: %d interactions yielded '
        '(%d partner_a non-matches, %d partner_b unresolved).',
        n_yielded, n_skipped_nonmatch, n_skipped_unresolved,
    )


def cellphonedb_ppi_interactions(
    organism: int = 9606,
    **_kwargs,
) -> Generator[Interaction, None, None]:
    """
    Yield CellPhoneDB v5.0 protein-protein interactions.

    Both ``partner_a`` and ``partner_b`` are resolved to UniProt ACs
    (expanding complexes via ``complex_input.csv``, one edge per member on
    each expanded side); MOR is derived from ``modulatory_effect``.

    Args:
        organism: NCBI taxonomy ID. CellPhoneDB is human-curated only;
            any organism other than 9606 yields nothing.

    Yields:
        :class:`~.._record.Interaction` records with
        ``source_type='protein'``, ``target_type='protein'``,
        ``id_type_a='uniprot'``, ``id_type_b='uniprot'``,
        ``interaction_type='ligand_receptor'``, ``resource='CellPhoneDB'``.
    """
    if organism != 9606:
        _log.info(
            '[COSMOS] CellPhoneDB PPI: human-curated only, skipping organism %d.',
            organism,
        )
        return

    interaction_rows, complex_rows = _load_cellphonedb_data()
    complex_map = _build_complex_map(complex_rows)

    n_yielded = 0
    n_skipped = 0

    for row in interaction_rows:
        if row.get('is_ppi') != 'True':
            continue

        partner_a = row.get('partner_a', '')
        partner_b = row.get('partner_b', '')
        sources = _resolve_entity(partner_a, complex_map)
        targets = _resolve_entity(partner_b, complex_map)

        if not sources or not targets:
            n_skipped += 1
            _log.info(
                "[COSMOS] CellPhoneDB PPI: could not resolve '%s' -> '%s', skipping.",
                partner_a, partner_b,
            )
            continue

        mor = _mor_from_modulatory_effect(row.get('modulatory_effect', ''))

        for source in sources:
            for target in targets:
                n_yielded += 1
                yield Interaction(
                    source=source,
                    target=target,
                    source_type='protein',
                    target_type='protein',
                    id_type_a='uniprot',
                    id_type_b='uniprot',
                    interaction_type='ligand_receptor',
                    resource='CellPhoneDB',
                    mor=mor,
                    locations=(),
                    attrs={},
                )

    _log.info(
        '[COSMOS] CellPhoneDB PPI: %d interactions yielded (%d unresolved rows skipped).',
        n_yielded, n_skipped,
    )
