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
Recon3D transporter interactions for COSMOS PKN.

Identifies transport reactions from Recon3D by detecting metabolites that
appear on both the reactant and product side of a reaction with different
compartment codes — the hallmark of membrane transport.  Only these
transport events are included; general metabolic reactions in Recon3D are
not covered here (Human-GEM via gem.py is the source for metabolic edges).

For each transport event two directed edges are generated per enzyme:

- *metabolite[in_comp] → enzyme*: metabolite enters the transporter.
- *enzyme → metabolite[out_comp]*: metabolite exits on the other side.

This produces a condensed representation: only metabolites that physically
cross a membrane generate edges.  Cofactors that are consumed and produced
within the same compartment (e.g. ATP→ADP in the cytoplasm) are excluded
even if they appear on both sides of the reaction.

Metabolite IDs are BiGG base IDs (compartment suffix stripped; id_type
``'bigg'``).  Enzyme IDs are Entrez Gene IDs with ``_ATN`` isoform suffixes
stripped (id_type ``'entrez'``).  Both are translated downstream by
:func:`~omnipath_metabo.datasets.cosmos._translate.translate_pkn`.

.. note::

    **Human-only resource.**  Recon3D is a human (*Homo sapiens*) metabolic
    reconstruction — no equivalent model exists for other species.  For
    multi-species transporter coverage, use ``gem.py`` with a Chalmers
    Sysbio / MetAtlas GEM (e.g. ``Mouse-GEM``, ``Rat-GEM``) and filter by
    ``resource='GEM_transporter:<gem-name>'``; those models carry an
    explicit ``'Transport reactions'`` subsystem annotation.

Transport direction is stored in ``attrs['transport_from']`` and
``attrs['transport_to']``.  The ``locations`` tuple carries the metabolite
compartment for each edge (source compartment for met→enzyme edges; target
compartment for enzyme→met edges).

Reference implementations:
    - ``recon3D_BIGG/script/recon3DBIGG_to_df.R`` — per-metabolite enzyme
      labelling for pure transport reactions (all reactants == all products
      by base ID).
    - ``Generation-PKN-COSMOS/src/final_functions_PKN_COSMOS.R``
      (``.format_GSMM_COSMOS``) — transporter detection by self-join on
      metabolite base IDs within each enzyme's edge set.
    - OmnipathR ``chalmers_gem.R`` (``chalmers_gem_network``) — self-join
      on (ri, source, target, reverse) to mark ``transporter = TRUE`` on
      edges where the same metabolite appears on both sides of a reaction.
"""

from __future__ import annotations

__all__ = ['recon3d_transporter_interactions']

import re
from collections import Counter, defaultdict
from collections.abc import Generator

from .._record import Interaction


def _strip_compartment(met_id: str) -> tuple[str, str]:
    """
    Split a BiGG metabolite ID into (base_id, compartment).

    BiGG encodes compartment as a single-letter suffix after the last
    underscore: ``'atp_c'`` → ``('atp', 'c')``.  If the last segment is
    not a single letter the full ID is returned with an empty string.

    Args:
        met_id: BiGG metabolite ID including compartment suffix.

    Returns:
        Tuple of (base_id, compartment_code).
    """
    parts = met_id.rsplit('_', 1)

    if len(parts) == 2 and len(parts[1]) == 1 and parts[1].isalpha():
        return parts[0], parts[1]

    return met_id, ''


def _parse_gene_rule(rule: str) -> list[str]:
    """
    Parse a Recon3D gene_reaction_rule into clean Entrez enzyme identifiers.

    Strips ``_ATN`` isoform suffixes from gene IDs before building complex
    subunit strings.  OR-separated parts become separate isoenzymes; AND-
    separated parts become underscore-joined complex subunit IDs.

    Args:
        rule: Gene-reaction rule string
            (e.g. ``'1234_AT1 or (5678_AT1 and 9012_AT2)'``).

    Returns:
        List of enzyme identifier strings with ``_ATN`` suffixes removed.
        Empty list for orphan reactions (empty or whitespace-only rule).
    """
    if not rule or not rule.strip():
        return []

    rule = re.sub(r'[()]', '', rule).strip()
    enzymes = []

    for or_part in re.split(r'\bor\b', rule, flags=re.IGNORECASE):
        genes = [
            re.sub(r'_AT\d+$', '', g.strip())
            for g in re.split(r'\band\b', or_part, flags=re.IGNORECASE)
            if g.strip()
        ]
        genes = [g for g in genes if g]

        if not genes:
            continue

        enzymes.append('_'.join(sorted(genes)) if len(genes) > 1 else genes[0])

    return enzymes


def recon3d_transporter_interactions(
    organism: int = 9606,
    metab_max_degree: int = 400,
    include_reverse: bool = True,
    include_orphans: bool = True,
) -> Generator[Interaction, None, None]:
    """
    Yield Recon3D transporter interactions as uniform Interaction records.

    Transport reactions are identified by detecting metabolites whose
    BiGG base ID (compartment suffix stripped) appears on both the reactant
    and product side of a reaction with *different* compartment codes.  This
    is the molecular signature of membrane transport.

    Only transported metabolites generate edges.  Metabolites that appear on
    both sides of a reaction but within the *same* compartment (e.g. a
    cofactor regenerated in the same location) are excluded.

    Two directed edges are generated per transported metabolite per enzyme:

    - ``met[in_comp] → enzyme``: the metabolite enters the transporter.
    - ``enzyme → met[out_comp]``: the metabolite exits on the other side.

    Reversible reactions produce an additional reversed pair with
    ``attrs['reverse'] = True`` where ``transport_from`` and ``transport_to``
    are swapped.

    High-degree metabolites are filtered in a two-pass approach: all
    candidate edges are collected first, metabolite degrees are counted, then
    only edges involving metabolites with degree ≤ *metab_max_degree* are
    yielded.

    Args:
        organism:
            NCBI taxonomy ID.  Accepted for API consistency with other
            resource processors; Recon3D is a human-only reconstruction
            so this parameter does not filter data.  Passing a non-human
            taxon ID does not produce species-specific output.
        metab_max_degree:
            Maximum number of edges a metabolite may participate in.
            Metabolites exceeding this threshold are treated as cofactors
            and removed.  Default: 400.
        include_reverse:
            If ``True``, include reversed edges for reversible transport
            reactions (``attrs['reverse'] = True``).  Default: ``True``.
        include_orphans:
            If ``True``, include transport reactions with no gene rule,
            using the reaction ID as a pseudo-enzyme node with
            ``id_type = 'reaction_id'`` and ``attrs['orphan'] = True``.
            Default: ``True``.

    Yields:
        :class:`~omnipath_metabo.datasets.cosmos._record.Interaction`
        records.  For metabolite → enzyme edges, ``source_type`` is
        ``'small_molecule'`` and ``id_type_a`` is ``'bigg'``; ``target_type``
        is ``'protein'`` and ``id_type_b`` is ``'entrez'``.  Roles are
        swapped for enzyme → metabolite edges.  ``interaction_type`` is
        ``'transport'`` and ``resource`` is ``'Recon3D'``.
    """
    from pypath.inputs.recon3d._gem import recon3d_reactions

    reactions = recon3d_reactions()

    # Two-pass: collect all candidate edge tuples first, then filter by
    # metabolite degree.  Tuple layout (11 fields):
    #   src, src_type, src_comp,
    #   tgt, tgt_type, tgt_comp,
    #   rxn_id, is_reverse, transport_from, transport_to, is_complex
    raw: list[tuple] = []
    n_transport_rxn = 0

    for rxn in reactions:
        enzymes = _parse_gene_rule(rxn['gene_reaction_rule'])

        if not enzymes and not include_orphans:
            continue

        mets: dict = rxn['metabolites']
        lb = rxn['lower_bound']
        ub = rxn['upper_bound']
        rxn_id = rxn['id']
        reversible = rxn['reversible']
        direction = 1 if lb + ub >= 0 else -1

        # Partition metabolites into reactants (consumed) and products
        # (produced), respecting the reaction direction sign.
        reactant_comps: defaultdict[str, list[str]] = defaultdict(list)
        product_comps: defaultdict[str, list[str]] = defaultdict(list)

        for met_full, coef in mets.items():
            base, comp = _strip_compartment(met_full)

            if coef * direction < 0:
                reactant_comps[base].append(comp)
            elif coef * direction > 0:
                product_comps[base].append(comp)

        # Transported metabolites: same base_id on both sides AND at least
        # one (in_comp, out_comp) pair with different compartments.
        transported: list[tuple[str, str, str]] = []

        for base in set(reactant_comps) & set(product_comps):
            for in_comp in reactant_comps[base]:
                for out_comp in product_comps[base]:
                    if in_comp != out_comp:
                        transported.append((base, in_comp, out_comp))

        if not transported:
            continue

        if not enzymes:
            # Orphan transport reaction: use reaction ID as pseudo-enzyme.
            n_transport_rxn += 1

            for base_id, in_comp, out_comp in transported:
                raw.append((
                    base_id, 'metabolite', in_comp,
                    rxn_id, 'reaction', '',
                    rxn_id, False, in_comp, out_comp, False,
                ))
                raw.append((
                    rxn_id, 'reaction', '',
                    base_id, 'metabolite', out_comp,
                    rxn_id, False, in_comp, out_comp, False,
                ))

                if reversible and include_reverse:
                    raw.append((
                        base_id, 'metabolite', out_comp,
                        rxn_id, 'reaction', '',
                        rxn_id, True, out_comp, in_comp, False,
                    ))
                    raw.append((
                        rxn_id, 'reaction', '',
                        base_id, 'metabolite', in_comp,
                        rxn_id, True, out_comp, in_comp, False,
                    ))

            continue

        n_transport_rxn += 1

        for base_id, in_comp, out_comp in transported:
            for enzyme in enzymes:
                is_complex = '_' in enzyme

                # Forward: met[in_comp] → enzyme → met[out_comp]
                raw.append((
                    base_id, 'metabolite', in_comp,
                    enzyme, 'protein', '',
                    rxn_id, False, in_comp, out_comp, is_complex,
                ))
                raw.append((
                    enzyme, 'protein', '',
                    base_id, 'metabolite', out_comp,
                    rxn_id, False, in_comp, out_comp, is_complex,
                ))

                if reversible and include_reverse:
                    # Reverse: met[out_comp] → enzyme → met[in_comp]
                    raw.append((
                        base_id, 'metabolite', out_comp,
                        enzyme, 'protein', '',
                        rxn_id, True, out_comp, in_comp, is_complex,
                    ))
                    raw.append((
                        enzyme, 'protein', '',
                        base_id, 'metabolite', in_comp,
                        rxn_id, True, out_comp, in_comp, is_complex,
                    ))

    # Count each metabolite's total degree across all candidate edges.
    metab_degree: Counter = Counter()

    for item in raw:
        if item[1] == 'metabolite':
            metab_degree[item[0]] += 1

        if item[4] == 'metabolite':
            metab_degree[item[3]] += 1

    n_filtered = 0

    for (
        src, src_type, src_comp,
        tgt, tgt_type, tgt_comp,
        rxn_id, is_reverse, transport_from, transport_to, is_complex,
    ) in raw:
        met = src if src_type == 'metabolite' else tgt

        if metab_degree[met] > metab_max_degree:
            n_filtered += 1
            continue

        is_orphan = src_type == 'reaction' or tgt_type == 'reaction'

        if src_type == 'metabolite':
            source = src
            source_type = 'small_molecule'
            id_type_a = 'bigg'
            target = tgt
            target_type = 'protein'
            id_type_b = 'reaction_id' if is_orphan else 'entrez'
            locations = (src_comp,) if src_comp else ()

        else:
            source = src
            source_type = 'protein'
            id_type_a = 'reaction_id' if is_orphan else 'entrez'
            target = tgt
            target_type = 'small_molecule'
            id_type_b = 'bigg'
            locations = (tgt_comp,) if tgt_comp else ()

        attrs = {
            'reverse': is_reverse,
            'reaction_id': rxn_id,
            'enzyme_complex': is_complex,
            'transport_from': transport_from,
            'transport_to': transport_to,
        }

        if is_orphan:
            attrs['orphan'] = True

        yield Interaction(
            source=source,
            target=target,
            source_type=source_type,
            target_type=target_type,
            id_type_a=id_type_a,
            id_type_b=id_type_b,
            interaction_type='transport',
            resource='Recon3D',
            mor=0,
            locations=locations,
            attrs=attrs,
        )
