# COSMOS PKN Quickstart

This tutorial shows how to build and explore the COSMOS prior-knowledge network
(PKN) using `omnipath-metabo`. The PKN covers metabolite-protein interactions,
protein-protein signaling (PPI), and gene regulatory networks (GRN) from nine
curated databases and is ready to use directly with the
[cosmosR](https://github.com/saezlab/cosmosR) R package or via the
`metabo.omnipathdb.org` web service.

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
from omnipath_metabo.datasets.cosmos import CosmosEdge

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

After `build_*()` the `network` attribute holds **`Interaction`** namedtuples
(fields: `source`, `target`, `source_type`, `target_type`, `id_type_a`,
`id_type_b`, `interaction_type`, `resource`, `mor`, `locations`, `attrs`).
After `format_*()` it holds **`CosmosEdge`** namedtuples (COSMOS node IDs,
no `id_type_a`/`id_type_b`).  In both cases `pd.DataFrame(bundle.network)`
produces correct column names automatically.

### Transporters

```python
# Human — no orphans, cell-surface only
transporters = cosmos.build_transporters(
    recon3d={'include_orphans': False},
    gem={'include_orphans': False},
    cell_surface_only=True,
)
df_t = pd.DataFrame(transporters.network)

print(df_t.groupby('resource').size())
print(df_t[df_t['resource'] == 'TCDB'])
print(df_t[df_t['resource'] == 'SLC'])
print(df_t[df_t['resource'] == 'Recon3D'])
print(df_t[df_t['resource'] == 'GEM_transporter:Human-GEM'])
print(df_t[df_t['resource'] == 'MRCLinksDB'])

# Inspect GEM transporter column types (Interaction fields, pre-format)
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
df_tm = pd.DataFrame(transporters_mouse.network)
df_gem_mouse = df_tm[df_tm['resource'] == 'GEM_transporter:Mouse-GEM']

# Inspect orphan transport reactions (no gene rule)
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
df_r = pd.DataFrame(receptors.network)
print(df_r.groupby('resource').size())
print(df_r[df_r['resource'] == 'STITCH'])
```

```python
# Mouse
receptors_mouse = cosmos.build_receptors(organism=10090, cell_surface_only=True)
df_rm = pd.DataFrame(receptors_mouse.network)
print(df_rm.groupby('resource').size())
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

### Enzyme–metabolite (GEM stoichiometric reactions)

```python
enzyme_met = cosmos.build_enzyme_metabolite(gem={'include_orphans': False})
df_e = pd.DataFrame(enzyme_met.network)
print(df_e.groupby('resource').size())
print(f'Unique reactions: {df_e["attrs"].apply(lambda a: a.get("reaction_id", "")).nunique()}')
```

### Protein–protein interactions (PPI)

Protein-protein signaling edges from the OmniPath database, covering
kinase-substrate, ligand-receptor, and other signaling interactions.

```python
from omnipath_metabo.datasets.cosmos._build import build_ppi

ppi = build_ppi()
df_ppi = pd.DataFrame(ppi.network)
print(df_ppi.groupby('resource').size())
print(f'PPI edges: {len(df_ppi):,}')

# Filter to stimulatory interactions only
stim = df_ppi[df_ppi['mor'] == 1]
print(f'Stimulatory: {len(stim):,}')
```

```python
# Mouse PPI (OmniPath supports multiple organisms)
ppi_mouse = build_ppi(organism=10090)
df_ppi_m = pd.DataFrame(ppi_mouse.network)
print(f'Mouse PPI edges: {len(df_ppi_m):,}')
```

### Gene regulatory network (GRN)

Transcription factor-target gene edges from CollecTRI (default) or DoRothEA.

```python
from omnipath_metabo.datasets.cosmos._build import build_grn

grn = build_grn()
df_grn = pd.DataFrame(grn.network)
print(df_grn.groupby('resource').size())
print(f'GRN edges: {len(df_grn):,}')

# All GRN edges have interaction_type='gene_regulation'
assert (df_grn['interaction_type'] == 'gene_regulation').all()
```

```python
# Include DoRothEA with confidence levels A,B,C
grn_ext = build_grn(
    ppi=False,
    resources={'grn': {'datasets': 'collectri,dorothea', 'dorothea_levels': 'A,B,C'}},
)
```

---

## 3. Format for cosmosR

After inspecting, format each bundle separately using the category-specific
format functions.  Node IDs are converted to COSMOS R package conventions:

- Metabolites: `Metab__CHEBI:XXXX_<comp1;comp2>` (all locations joined with `;`)
- Proteins (forward): `Gene<N>__<UniProtAC>` or `Gene<N>__<AC1;AC2>` for ambiguous mappings
- Proteins (reverse): `Gene<N>__<UniProtAC>_rev`
- Orphan reactions (no gene rule): `Gene<N>__orphanReac<reaction_id>`
- Receptor / allosteric proteins: bare `<UniProtAC>` (no `Gene<N>__` prefix)

`N` is a sequential reaction index that resets per category (transporters, receptors, etc.).

Connector edges (`interaction_type='connector'`) are added for transporter and GEM gene
nodes only — linking the bare UniProt accession to the formatted `Gene<N>__` node.
No metabolite connector edges are emitted in any category.

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
    df[['source', 'target', 'mor']].rename(columns={'mor': 'sign'}).to_csv(
        f'cosmos_pkn_{name}.csv', index=False,
    )
```

Alternatively, use `to_dataframes()` to get all bundle components at once:

```python
dfs = fmt_t.to_dataframes()
net_df  = dfs['network']      # CosmosEdge rows
met_df  = dfs['metabolites']  # CosmosMetabolite rows
prot_df = dfs['proteins']     # CosmosProtein rows
rxn_df  = dfs['reactions']    # CosmosReaction rows
```

---

## 4. CLI usage

For scripted or one-off exports, use the `cosmos-pkn` command installed with the package:

```bash
# Full human PKN (3-column CSV: source, target, sign)
cosmos-pkn

# Transporters only, TSV output
cosmos-pkn --subset transporters --output cosmos_transporters.tsv

# Mouse, no STITCH, no connector edges
cosmos-pkn --organism 10090 --no-stitch --no-connector-edges

# All columns (source, target, sign, interaction_type, resource, ...)
cosmos-pkn --all-columns --output cosmos_pkn_full.csv

# Higher STITCH confidence, drop pseudo-enzyme orphan nodes
cosmos-pkn --score-threshold 900 --no-orphans
```

Available subsets: `all` (default), `transporters`, `receptors`, `allosteric`, `enzyme_metabolite`, `ppi`, `grn`.

To pre-build Parquet cache files for the server:

```bash
cosmos-pkn build-cache --organism 9606 10090 --category transporters receptors ppi grn
```

---

## 5. Customise the build

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

**Override config via YAML file:**

```python
cfg = cosmos.config('my_config.yaml', organism=10090)
transporters = cosmos.build_transporters(**cfg)
```

---

## 6. Provenance lookup

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

## 7. Multi-organism builds

The `build()` function and all `build_*()` category builders accept an `organism`
parameter (NCBI taxonomy ID). Resources that natively support the organism
query it directly; human-only resources (SLC, Recon3D) are translated via
orthology.

```python
# Full mouse PKN (all categories)
from omnipath_metabo.datasets.cosmos import build

mouse_pkn = build(organism=10090)
df_mouse = pd.DataFrame(mouse_pkn.network)
print(f'Mouse PKN: {len(df_mouse):,} edges')
print(df_mouse.groupby('resource').size())
```

```python
# Check which resources need orthology translation for an organism
from omnipath_metabo.datasets.cosmos._organisms import (
    needs_orthology,
    organism_resources,
    default_gem,
)

# Mouse: most resources are direct, SLC + Recon3D need orthology
print(organism_resources(10090))
print(f'Needs orthology: {needs_orthology(10090)}')
print(f'Default GEM: {default_gem(10090)}')
```

Supported organisms with GEMs: human (9606), mouse (10090), rat (10116),
zebrafish (7955), fruit fly (7227), worm (6239), yeast (4932), *E. coli* (562).
Resources marked "all" (TCDB, BRENDA, STITCH) accept any organism.

---

## 8. Web service

Pre-built PKNs are available from the `metabo.omnipathdb.org` web service.
The server caches each category per organism as Parquet files and assembles
them on request.

```
GET https://metabo.omnipathdb.org/cosmos/pkn?organism=9606&categories=all
GET https://metabo.omnipathdb.org/cosmos/pkn?organism=10090&categories=transporters,ppi
GET https://metabo.omnipathdb.org/cosmos/pkn?organism=9606&format=parquet
```

Additional endpoints:

```
GET /cosmos/categories    # list available categories
GET /cosmos/organisms     # list organisms with cached data
GET /cosmos/resources     # list resources per category
GET /cosmos/status        # cache status and sizes
```

---

## 9. Python client (`omnipath-client`)

The `omnipath-client` package provides a high-level interface to the web service.
Install with `pip install omnipath-client`.

```python
import omnipath_client as oc

# Fetch human transporters as a DataFrame (polars by default, pandas fallback)
df = oc.cosmos.get_pkn(organism='human', categories='transporters')

# Full mouse PKN
df_mouse = oc.cosmos.get_pkn(organism=10090)

# Specific categories
df_subset = oc.cosmos.get_pkn(
    organism=9606,
    categories=['transporters', 'receptors', 'ppi'],
)

# Filter by resource
df_signor = oc.cosmos.get_pkn(
    organism=9606,
    categories='ppi',
    resources='OmniPath:omnipath',
)

# Raw dict (JSON response)
data = oc.cosmos.get_pkn(organism=9606, format='dict')

# Explore available data
print(oc.cosmos.categories())    # ['transporters', 'receptors', ...]
print(oc.cosmos.organisms())     # [9606, 10090, ...]
print(oc.cosmos.resources())     # {'transporters': ['GEM:Human-GEM', ...], ...}
```

The `organism` parameter accepts any form recognised by the OmniPath taxonomy
service: NCBI ID (`9606`), common name (`'human'`), Latin name
(`'Homo sapiens'`), Ensembl code (`'hsapiens'`), or KEGG code (`'hsa'`).

---

## 10. AnnNet graph conversion

The PKN can be converted to an [AnnNet](https://github.com/saezlab/annnet)
graph object for network analysis. This works with both the local build
output and the web service client.

```python
# Via the web service client (recommended for quick access)
import omnipath_client as oc

g = oc.cosmos.get_pkn(organism='human', format='annnet')

# Or convert an existing DataFrame
df = oc.cosmos.get_pkn(organism='human')
g = oc.cosmos.to_annnet(df)
```

```python
# From a local build
from omnipath_metabo.datasets.cosmos import build, format_transporters

bundle = build()
fmt = format_transporters(bundle)
df_local = pd.DataFrame(fmt.network)

# Convert to AnnNet (requires the annnet package)
from omnipath_client.cosmos import to_annnet

g = to_annnet(df_local)
```

The AnnNet graph stores entity types (protein, small_molecule) as vertex
attributes and interaction metadata (interaction_type, resource, MOR as weight)
as edge attributes.

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
| OmniPath PPI | Protein-protein signaling | Multi-species |
| OmniPath GRN | Gene regulation (CollecTRI/DoRothEA) | Multi-species |
