# COSMOS PKN Formatter — Implementation Notes

This document captures design constraints for the future COSMOS formatter step.
The PKN builder (``build()``, ``build_transporters()``, etc.) is already complete;
the formatter applies node-ID prefixes/suffixes required by the COSMOS R package.

---

## Reverse-edge convention

Both forward **and** reverse edges are already present in the PKN output when
``include_reverse=True`` (the default for ``gem_interactions`` and
``recon3d_transporter_interactions``).

Reverse edges are identified by ``attrs['reverse'] == True``.

**The formatter MUST NOT re-generate reverse edges.**  Doing so would double
the edge count and produce an incorrect network.  The correct approach is:

1. Transform node IDs for all rows (both forward and reverse).
2. Use ``attrs['reverse']`` only to apply directional labels or annotations,
   **not** to create new edges.
3. Set ``attrs['cosmos_formatted'] = True`` after formatting to prevent
   accidental double-application.

---

## Required node-ID transformations

| Entity type | Transformation |
|---|---|
| Metabolite (ChEBI) | Add ``Metab__`` prefix and ``__<compartment>`` suffix |
| Protein (ENSG) | Add ``Gene<N>__`` prefix where N is a **sequential reaction index** (1, 2, 3, …) assigned per unique ``reaction_id``; add ``_rev`` suffix when ``attrs['reverse'] == True`` |
| Orphan pseudo-enzyme (reaction_id) | Retain for now; removable via ``include_orphans`` flag |

The reaction index makes the same gene a distinct node in each reaction it
catalyses, enabling the PKN to represent the same protein acting in multiple
independent reactions without creating spurious loops.

Example (reaction index 42):
- ``CHEBI:15422`` in cytoplasm → ``Metab__CHEBI:15422__c``
- ``ENSG00000141510`` (forward) → ``Gene42__ENSG00000141510``
- ``ENSG00000141510`` (reverse) → ``Gene42__ENSG00000141510_rev``

---

## Compartment information

Compartment codes are stored in:
- ``locations`` tuple on each ``Interaction`` record (one or two entries)
- ``attrs['transport_from']`` / ``attrs['transport_to']`` for Recon3D edges
- GEM metabolite IDs carry a compartment suffix pre-stripping (accessible via
  the ``locations`` column)

---

## Implementation checklist (when formatter is written)

- [ ] Apply node-ID prefixes/suffixes to ``source`` and ``target`` columns
- [ ] Set ``attrs['cosmos_formatted'] = True`` on each row
- [ ] Do **not** call ``include_reverse`` or expand edges — the PKN already
      has both directions
- [ ] Validate that no ``(source, target, reaction_id, reverse)`` key appears
      more than once after formatting (the safety test in ``test_gem.py``
      confirms this pre-formatting)
- [ ] Handle orphan pseudo-enzymes (``id_type == 'reaction_id'``) — decide
      whether to include, drop, or relabel them in the formatted network
