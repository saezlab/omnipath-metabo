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

## 1. Build the full PKN

The default build collects all seven resources for human (NCBI taxonomy ID
9606), translates metabolite IDs to ChEBI and protein IDs to UniProt, and
applies the expert-curation blacklist.

```python
import logging
import warnings
import pandas as pd
from omnipath_metabo.datasets import cosmos
from omnipath_metabo.datasets.cosmos._record import Interaction

warnings.filterwarnings('ignore', module='paramiko')
warnings.filterwarnings('ignore', module='rdata')

logging.basicConfig(level=logging.WARNING, format='%(message)s')
  
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# transporter
transporters = cosmos.build_transporters(
    recon3d={'include_orphans':False},
    gem={'include_orphans': False})
df = pd.DataFrame(transporters.network, columns=Interaction._fields)
df.groupby(['resource']).size()
df[df['resource']=='TCDB']
df[df['resource']=='SLC']
df[df['resource']=='Recon3D']
df[df['resource']=='GEM_transporter:Human-GEM']
df_r = df[df['resource'] == 'GEM_transporter:Human-GEM']    
for col in ['source_type', 'target_type', 'id_type_a',            
    'id_type_b','interaction_type']:  
    print(df_r[col].value_counts())









bundle = cosmos.build()

print(f'Network edges : {len(bundle.network):,}')
print(f'Metabolites   : {len(bundle.metabolites):,}')
print(f'Proteins      : {len(bundle.proteins):,}')
print(f'GEM reactions : {len(bundle.reactions):,}')
```

The returned `CosmosBundle` has four components:

| Attribute | Content |
|---|---|
| `bundle.network` | List of `Interaction` namedtuples — the edges |
| `bundle.metabolites` | ChEBI to original-ID provenance per metabolite |
| `bundle.proteins` | UniProt to original-ID provenance per protein |
| `bundle.reactions` | GEM reaction metadata (genes, metabolites, subsystem) |

**Note:** The first run downloads and caches data. STITCH and the GEMs are
large and may take a few minutes. Subsequent runs use the local cache and are
much faster.

---

## 2. Inspect the network as a DataFrame

```python
import pandas as pd
from omnipath_metabo.datasets.cosmos._record import Interaction

df = pd.DataFrame(bundle.network, columns=Interaction._fields)
print(df.head())
```

Column reference:

| Column | Description |
|---|---|
| `source` | Source entity ID (ChEBI or UniProt after translation) |
| `target` | Target entity ID |
| `source_type` | `'small_molecule'` or `'protein'` |
| `target_type` | `'small_molecule'` or `'protein'` |
| `id_type_a` | ID namespace of source after translation |
| `id_type_b` | ID namespace of target after translation |
| `interaction_type` | `'transport'`, `'ligand_receptor'`, `'allosteric_regulation'`, ... |
| `resource` | Source database, e.g. `'TCDB'`, `'STITCH'`, `'GEM:Human-GEM'`, `'Recon3D'` |
| `mor` | Mode of regulation: `1` (activation), `-1` (inhibition), `0` (unknown) |
| `locations` | Tuple of compartment codes, e.g. `('e', 'c')` |
| `attrs` | Dict with extra metadata: `reaction_id`, `reverse`, `transport_from`, ... |

Summarise by resource and interaction type:

```python
print(df.groupby(['resource', 'interaction_type']).size().to_string())
```

---

## 3. Category-specific subsets

Four convenience functions build subsets without loading irrelevant resources:

```python
# Transporters — TCDB, SLC, GEM_transporter, Recon3D, STITCH-transporter
transporters = cosmos.build_transporters()
df_t = pd.DataFrame(transporters.network, columns=Interaction._fields)
print(f'Transporter edges : {len(df_t):,}')

# Receptors — MRCLinksDB + STITCH-receptor
receptors = cosmos.build_receptors()
df_r = pd.DataFrame(receptors.network, columns=Interaction._fields)
print(f'Receptor edges    : {len(df_r):,}')

# Allosteric regulation — BRENDA + STITCH-other
allosteric = cosmos.build_allosteric()
df_a = pd.DataFrame(allosteric.network, columns=Interaction._fields)
print(f'Allosteric edges  : {len(df_a):,}')

# Enzyme-metabolite — Human-GEM stoichiometric reactions only
enzyme_met = cosmos.build_enzyme_metabolite()
df_e = pd.DataFrame(enzyme_met.network, columns=Interaction._fields)
print(f'Enzyme-metab edges: {len(df_e):,}')
```

---

## 4. Format for cosmosR

`format_pkn()` applies the node-ID prefixes and suffixes expected by the
COSMOS R package:

- Metabolites: `Metab__CHEBI:XXXX_<compartment>`
- Proteins (forward): `Gene<N>__<UniProtAC>`
- Proteins (reverse): `Gene<N>__<UniProtAC>_rev`

```python
from omnipath_metabo.datasets.cosmos import format_pkn

formatted = format_pkn(bundle)
df_fmt = pd.DataFrame(formatted.network)
print(df_fmt[['source', 'target', 'mor', 'interaction_type', 'resource']].head(10))

# Save — ready to load in cosmosR::preprocess_COSMOS_*
df_fmt[['source', 'target', 'mor']].to_csv('cosmos_pkn.csv', index=False)
```

---

## 5. Customise the build

**Change organism (mouse):**

```python
bundle_mouse = cosmos.build(organism=10090)
```

**Lower the STITCH confidence threshold:**

```python
bundle_loose = cosmos.build(stitch={'score_threshold': 500})
```

**Disable the GEMs for a faster build:**

```python
bundle_fast = cosmos.build(gem=False, recon3d=False)
```

**TCDB and SLC only, without ID translation:**

```python
bundle_minimal = cosmos.build(
    brenda=False,
    mrclinksdb=False,
    stitch=False,
    gem=False,
    recon3d=False,
    translate_ids=False,
)
```

---

## 6. Provenance lookup

The bundle tracks the mapping between translated canonical IDs and the original
identifiers from each source database.

```python
# Metabolites: ChEBI <- original source ID
df_met = pd.DataFrame(bundle.metabolites)
# columns: chebi | original_id | id_type | resource | name
print(df_met.head())

# Proteins: UniProt <- original ID (ENSP, gene symbol, ...)
df_prot = pd.DataFrame(bundle.proteins)
# columns: uniprot | original_id | id_type | resource | gene_symbol
print(df_prot.head())

# GEM reactions
df_rxn = pd.DataFrame(bundle.reactions)
# columns: reaction_id | gem | subsystem | genes | metabolites
print(df_rxn.head())
```

---

## Default resources

| Resource | Category | Species |
|---|---|---|
| TCDB | Transporters | Multi-species |
| SLC | Transporters | Human only |
| STITCH | Transporters, receptors, allosteric | Multi-species |
| MRCLinksDB | Receptors | Human, mouse |
| BRENDA | Allosteric regulation | Multi-species |
| Human-GEM | Enzyme-metabolite, GEM transporters | Human |
| Recon3D | GEM transporters | Human |
