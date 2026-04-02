# COSMOS PKN Quickstart

This tutorial shows how to build and explore the COSMOS prior-knowledge network
(PKN) using `omnipath-metabo`. The PKN covers metabolite-protein interactions
from seven curated databases and is ready to use directly with the
[cosmosR](https://github.com/saezlab/cosmosR) R package.

## Prerequisites

Run all commands from the `omnipath-metabo` project directory:

```bash
cd ~/biotools/omnipath-metabo
```

Execute scripts with:

```bash
uv run python your_script.py
```

Or start an interactive session:

```bash
uv run python
```

---

## 1. Setup

```python
import logging
import warnings

import pandas as pd

from omnipath_metabo.datasets import cosmos
from omnipath_metabo.datasets.cosmos._record import Interaction

warnings.filterwarnings('ignore', module='paramiko')
warnings.filterwarnings('ignore', module='rdata')
logging.basicConfig(level=logging.INFO, format='%(message)s')

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
```

---

## 2. Build and inspect each category

Build each PKN category separately and inspect before formatting.
Each `build_*()` call returns a `CosmosBundle` with translated ChEBI / UniProt IDs.

### Transporters

```python
# Human — no orphans, cell-surface only
transporters = cosmos.build_transporters(
    recon3d={'include_orphans': False},
    gem={'include_orphans': False},
    cell_surface_only=True,
)
df_t = pd.DataFrame(transporters.network, columns=Interaction._fields)

print(df_t.groupby('resource').size())
print(df_t[df_t['resource'] == 'TCDB'])
print(df_t[df_t['resource'] == 'SLC'])
print(df_t[df_t['resource'] == 'Recon3D'])
print(df_t[df_t['resource'] == 'GEM_transporter:Human-GEM'])
print(df_t[df_t['resource'] == 'MRCLinksDB'])

# Inspect GEM transporter column types
df_gem = df_t[df_t['resource'] == 'GEM_transporter:Human-GEM']
for col in ['source_type', 'target_type', 'id_type_a', 'id_type_b', 'interaction_type']:
    print(df_gem[col].value_counts())
```

```python
# Mouse — include orphans, cell-surface only
transporters_mouse = cosmos.build_transporters(
    recon3d={'include_orphans': True},
    gem={'include_orphans': True, 'gem': 'Mouse-GEM'},
    organism=10090,
    cell_surface_only=True,
)
df_tm = pd.DataFrame(transporters_mouse.network, columns=Interaction._fields)

# Inspect orphan transport reactions (no gene rule)
df_gem_mouse = df_tm[df_tm['resource'] == 'GEM_transporter:Mouse-GEM']
orphans = df_gem_mouse[
    df_gem_mouse['attrs'].apply(lambda a: isinstance(a, dict) and a.get('orphan', False)) &
    (df_gem_mouse['source_type'] == 'small_molecule')
]
print(f'Orphan transport reactions (Mouse-GEM): {orphans["attrs"].apply(lambda a: a["reaction_id"]).nunique()}')
print(orphans[['source', 'target', 'locations', 'attrs']].head())
```

### Receptors

```python
# Human
receptors = cosmos.build_receptors(cell_surface_only=True)
df_r = pd.DataFrame(receptors.network, columns=Interaction._fields)
print(df_r.groupby('resource').size())
print(df_r[df_r['resource'] == 'STITCH'])
```

```python
# Mouse
receptors_mouse = cosmos.build_receptors(organism=10090, cell_surface_only=True)
df_rm = pd.DataFrame(receptors_mouse.network, columns=Interaction._fields)
print(df_rm.groupby('resource').size())
```

### Allosteric regulation

```python
# Default: BRENDA + STITCH-other
allosteric = cosmos.build_allosteric()
df_a = pd.DataFrame(allosteric.network, columns=Interaction._fields)
print(df_a.groupby(['resource', 'interaction_type']).size())

# Higher STITCH confidence threshold
allosteric_strict = cosmos.build_allosteric(stitch={'score_threshold': 900})

# BRENDA only
allosteric_brenda = cosmos.build_allosteric(stitch=False)
```

### Enzyme–metabolite (GEM stoichiometric reactions)

```python
enzyme_met = cosmos.build_enzyme_metabolite(gem={'include_orphans': False})
df_e = pd.DataFrame(enzyme_met.network, columns=Interaction._fields)
print(df_e.groupby('resource').size())
print(f'Unique reactions: {df_e["attrs"].apply(lambda a: a.get("reaction_id", "")).nunique()}')
```

---

## 3. Format for cosmosR

After inspecting, format each bundle separately using the category-specific
format functions.  Node IDs are converted to COSMOS R package conventions:

- Metabolites: `Metab__CHEBI:XXXX_<compartment>`
- Proteins (forward): `Gene<N>__<UniProtAC>`
- Proteins (reverse): `Gene<N>__<UniProtAC>_rev`
- Orphan reactions (no gene): `Rxn<N>__<reaction_id>`

```python
from omnipath_metabo.datasets.cosmos import (
    format_transporters,
    format_receptors,
    format_allosteric,
    format_enzyme_metabolite,
)

fmt_t  = format_transporters(transporters)
fmt_r  = format_receptors(receptors)
fmt_a  = format_allosteric(allosteric)
fmt_e  = format_enzyme_metabolite(enzyme_met)

for name, bundle in [('transporters', fmt_t), ('receptors', fmt_r),
                     ('allosteric', fmt_a), ('enzyme_met', fmt_e)]:
    df = pd.DataFrame(bundle.network)
    main = df[df['interaction_type'] != 'connector']
    print(f'{name}: {len(main):,} edges, {len(df):,} with connectors')
```

Export each category — ready to load in `cosmosR::preprocess_COSMOS_*`:

```python
for name, bundle in [('transporters', fmt_t), ('receptors', fmt_r),
                     ('allosteric', fmt_a), ('enzyme_met', fmt_e)]:
    df = pd.DataFrame(bundle.network)
    df[['source', 'target', 'mor']].to_csv(f'cosmos_pkn_{name}.csv', index=False)
```

---

## 4. Customise the build

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

**TCDB and SLC only, without ID translation:**

```python
transporters_minimal = cosmos.build_transporters(
    gem=False,
    recon3d=False,
    mrclinksdb=False,
    translate_ids=False,
)
```

---

## 5. Provenance lookup

The bundle tracks the mapping between translated canonical IDs and the original
identifiers from each source database.

```python
# Metabolites: ChEBI ← original source ID
df_met = pd.DataFrame(transporters.metabolites)
# columns: chebi | original_id | id_type | resource | name
print(df_met.head())

# Proteins: UniProt ← original ID (ENSP, Ensembl, Entrez, ...)
df_prot = pd.DataFrame(transporters.proteins)
# columns: uniprot | original_id | id_type | resource | gene_symbol
print(df_prot.head())

# GEM reactions
df_rxn = pd.DataFrame(transporters.reactions)
# columns: reaction_id | gem | subsystem | genes | metabolites
print(df_rxn.head())
```

---

## Default resources

| Resource | Category | Species |
|---|---|---|
| TCDB | Transporters | Multi-species |
| SLC | Transporters | Human only |
| MRCLinksDB | Transporters, receptors | Human, mouse |
| STITCH | Receptors, allosteric | Multi-species |
| BRENDA | Allosteric regulation | Multi-species |
| Human-GEM | Transporters, enzyme-metabolite | Human |
| Recon3D | Transporters | Human |
