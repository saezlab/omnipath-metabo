# COSMOS PKN — Future Features

---

## 1. Intercellular / cell-cell communication parameter

Add a parameter to `build_transporters()` (and possibly `build()`) to filter
for proteins involved in intercellular cell-cell communication — i.e. proteins
localised at the cell surface / extracellular space.

**Filter logic:** keep only edges where `'e'` appears in the `locations` tuple.
Both extracellular-facing and plasma membrane proteins are annotated as `'e'`
in the BiGG compartment convention used by the `locations` field.

**Proposed API:**

```python
transporters = cosmos.build_transporters(cell_surface_only=True)
```

**Implementation:** add `cell_surface_only: bool = False` parameter to
`build_transporters()`; when `True`, post-filter to keep only rows where
`'e' in row.locations`.

---

## 2. Remove duplicated entries

Deduplicate edges in the PKN output. The same `(source, target, mor, resource)`
combination can appear multiple times when a metabolite-protein pair is covered
by multiple reactions or annotation entries within the same resource.

**Implementation:** after collecting all edges in `build()` and after ID
translation, deduplicate the DataFrame by `(source, target, mor, resource)`,
keeping the first occurrence. Consider whether `locations` and `attrs` should
be merged (union of locations, merged attrs) or simply dropped on duplicates.

**Important:** Recon3D and GEM transport reactions are bidirectional — both
forward and reverse edges are present by default (`include_reverse=True`).
Reverse edges are identified by `attrs['reverse'] == True`. Deduplication
must treat forward and reverse edges as distinct entries and must NOT collapse
them. The deduplication key should include the `reverse` flag:
`(source, target, mor, resource, attrs.get('reverse', False))`.
