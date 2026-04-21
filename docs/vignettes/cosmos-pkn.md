# COSMOS PKN

This tutorial shows how to build and explore the COSMOS prior-knowledge network
(PKN) using omnipath-metabo. The PKN covers metabolite-protein interactions,
protein-protein signaling (PPI), and gene regulatory networks (GRN) from nine
curated databases and is ready to use directly with the
[cosmosR](https://github.com/saezlab/cosmosR) R package or via the
[metabo.omnipathdb.org](https://metabo.omnipathdb.org) web service.

## Setup

```python
import pandas as pd
from omnipath_metabo.datasets import cosmos
```

---

## Build and inspect each category

Build each PKN category separately and inspect before formatting.
Each `build_*()` call returns a `CosmosBundle` with translated ChEBI / UniProt IDs.

### Transporters

```python
# Human -- no orphans, cell-surface only
transporters = cosmos.build_transporters(
    recon3d={'include_orphans': False},
    gem={'include_orphans': False},
    cell_surface_only=True,
)
df_t = pd.DataFrame(transporters.network)
print(df_t.groupby('resource').size())
```

```python
# Mouse -- include orphans, cell-surface only
transporters_mouse = cosmos.build_transporters(
    recon3d={'include_orphans': True},
    gem={'include_orphans': True, 'gem': 'Mouse-GEM'},
    organism=10090,
    cell_surface_only=True,
)
df_tm = pd.DataFrame(transporters_mouse.network)
```

### Receptors

```python
receptors = cosmos.build_receptors(cell_surface_only=True)
df_r = pd.DataFrame(receptors.network)
print(df_r.groupby('resource').size())
```

### Allosteric regulation

```python
# Default: BRENDA + STITCH-other
allosteric = cosmos.build_allosteric()
df_a = pd.DataFrame(allosteric.network)
print(df_a.groupby(['resource', 'interaction_type']).size())

# Higher STITCH confidence threshold
allosteric_strict = cosmos.build_allosteric(stitch={'score_threshold': 900})

# BRENDA only
allosteric_brenda = cosmos.build_allosteric(stitch=False)
```

### Enzyme-metabolite (GEM stoichiometric reactions)

```python
enzyme_met = cosmos.build_enzyme_metabolite(gem={'include_orphans': False})
df_e = pd.DataFrame(enzyme_met.network)
print(df_e.groupby('resource').size())
```

### Protein-protein interactions (PPI)

```python
from omnipath_metabo.datasets.cosmos._build import build_ppi

ppi = build_ppi()
df_ppi = pd.DataFrame(ppi.network)
print(f'PPI edges: {len(df_ppi):,}')
```

### Gene regulatory network (GRN)

```python
from omnipath_metabo.datasets.cosmos._build import build_grn

grn = build_grn()
df_grn = pd.DataFrame(grn.network)
print(f'GRN edges: {len(df_grn):,}')
```

---

## Format for cosmosR

Format each bundle with COSMOS node ID conventions:

- Metabolites: `Metab__CHEBI:XXXX_<comp1;comp2>`
- Proteins (forward): `Gene<N>__<UniProtAC>`
- Proteins (reverse): `Gene<N>__<UniProtAC>_rev`
- Orphan reactions: `Gene<N>__orphanReac<reaction_id>`
- Receptor / allosteric proteins: bare `<UniProtAC>`

```python
from omnipath_metabo.datasets.cosmos import (
    format_transporters,
    format_receptors,
    format_allosteric,
    format_enzyme_metabolite,
)

fmt_t = format_transporters(transporters)
fmt_r = format_receptors(receptors)
fmt_a = format_allosteric(allosteric)
fmt_e = format_enzyme_metabolite(enzyme_met)

for name, bundle in [('transporters', fmt_t), ('receptors', fmt_r),
                     ('allosteric', fmt_a), ('enzyme_met', fmt_e)]:
    df = pd.DataFrame(bundle.network)
    main = df[df['interaction_type'] != 'connector']
    print(f'{name}: {len(main):,} edges, {len(df):,} with connectors')
```

Export to CSV:

```python
for name, bundle in [('transporters', fmt_t), ('receptors', fmt_r),
                     ('allosteric', fmt_a), ('enzyme_met', fmt_e)]:
    df = pd.DataFrame(bundle.network)
    df[['source', 'target', 'mor']].rename(columns={'mor': 'sign'}).to_csv(
        f'cosmos_pkn_{name}.csv', index=False,
    )
```

---

## CLI usage

```bash
# Full human PKN (3-column CSV: source, target, sign)
cosmos-pkn

# Transporters only, TSV output
cosmos-pkn --subset transporters --output cosmos_transporters.tsv

# Mouse, no STITCH, no connector edges
cosmos-pkn --organism 10090 --no-stitch --no-connector-edges

# All columns
cosmos-pkn --all-columns --output cosmos_pkn_full.csv
```

Available subsets: `all` (default), `transporters`, `receptors`, `allosteric`,
`enzyme_metabolite`, `ppi`, `grn`.

---

## Customise the build

**Change organism (mouse):**

```python
receptors_mouse = cosmos.build_receptors(organism=10090, cell_surface_only=True)
```

**Lower the STITCH confidence threshold:**

```python
allosteric_loose = cosmos.build_allosteric(stitch={'score_threshold': 500})
```

**Disable GEMs for a faster transporter build:**

```python
transporters_fast = cosmos.build_transporters(gem=False, recon3d=False)
```

---

## Multi-organism builds

The `build()` function and all `build_*()` category builders accept an `organism`
parameter (NCBI taxonomy ID). Resources that natively support the organism
query it directly; human-only resources (SLC, Recon3D) are translated via
orthology.

```python
from omnipath_metabo.datasets.cosmos import build

mouse_pkn = build(organism=10090)
df_mouse = pd.DataFrame(mouse_pkn.network)
print(f'Mouse PKN: {len(df_mouse):,} edges')
```

Supported organisms with GEMs: human (9606), mouse (10090), rat (10116),
zebrafish (7955), fruit fly (7227), worm (6239), yeast (4932), *E. coli* (562).

---

## Web service

Pre-built PKNs are available from the
[metabo.omnipathdb.org](https://metabo.omnipathdb.org) web service:

```
GET https://metabo.omnipathdb.org/cosmos/pkn?organism=9606&categories=all
GET https://metabo.omnipathdb.org/cosmos/pkn?organism=10090&categories=transporters,ppi
GET https://metabo.omnipathdb.org/cosmos/categories
GET https://metabo.omnipathdb.org/cosmos/organisms
GET https://metabo.omnipathdb.org/cosmos/resources
```

---

## Python client (omnipath-client)

The [omnipath-client](https://github.com/saezlab/omnipath-client) package
provides a high-level interface to the web service:

```python
import omnipath_client as oc

# Full human PKN
df = oc.cosmos.get_pkn('human')

# Mouse transporters and PPI only
df = oc.cosmos.get_pkn(10090, categories=['transporters', 'ppi'])

# As an AnnNet graph
g = oc.cosmos.get_pkn('human', format='annnet')

# Explore available data
print(oc.cosmos.categories())
print(oc.cosmos.organisms())
print(oc.cosmos.resources())
```

---

## Default resources

| Resource | Category | Species |
|----------|----------|---------|
| TCDB | Transporters | Multi-species |
| SLC | Transporters | Human only |
| MRCLinksDB | Transporters, receptors | Human, mouse |
| STITCH | Receptors, allosteric | Multi-species |
| BRENDA | Allosteric regulation | Multi-species |
| Human-GEM | Transporters, enzyme-metabolite | Human |
| Recon3D | Transporters | Human |
| OmniPath PPI | Protein-protein signaling | Multi-species |
| OmniPath GRN | Gene regulation (CollecTRI/DoRothEA) | Multi-species |
