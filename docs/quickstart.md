# Quickstart

This guide covers the three main workflows: building a PKN locally,
formatting it for cosmosR, and serving it via the web service.

## Setup

```python
from omnipath_metabo.datasets import cosmos
```

## Build a PKN

Build the full human COSMOS PKN (all 6 categories):

```python
pkn = cosmos.build()
```

Or build individual categories:

```python
transporters = cosmos.build_transporters()
receptors = cosmos.build_receptors()
allosteric = cosmos.build_allosteric()
enzyme_met = cosmos.build_enzyme_metabolite()
```

The PPI and GRN layers are also available:

```python
from omnipath_metabo.datasets.cosmos._build import build_ppi, build_grn

ppi = build_ppi()
grn = build_grn()
```

## Inspect the result

Each builder returns a `CosmosBundle` with a `network` attribute.
Convert to a DataFrame for inspection:

```python
import pandas as pd

df = pd.DataFrame(transporters.network)
print(df.groupby('resource').size())
```

## Format for cosmosR

Format the network with COSMOS node ID conventions:

```python
from omnipath_metabo.datasets.cosmos import (
    format_transporters,
    format_receptors,
    format_allosteric,
    format_enzyme_metabolite,
)

fmt = format_transporters(transporters)
df = pd.DataFrame(fmt.network)
df[['source', 'target', 'mor']].rename(columns={'mor': 'sign'}).to_csv(
    'cosmos_pkn_transporters.csv', index=False,
)
```

## Multi-organism

Build for mouse or rat by passing the NCBI taxonomy ID:

```python
mouse_pkn = cosmos.build(organism=10090)
```

## CLI export

Export from the command line:

```bash
# Full human PKN
cosmos-pkn

# Mouse transporters only
cosmos-pkn --organism 10090 --subset transporters

# All columns, TSV output
cosmos-pkn --all-columns --output cosmos_pkn.tsv
```

## Serve via the web service

Start a local instance of the metabo web service:

```bash
pip install "omnipath-metabo[server]"
omnipath-metabo serve
```

Then query it:

```
GET http://localhost:8000/cosmos/pkn?organism=9606&categories=all
```

Or use the Python client:

```python
import omnipath_client as oc

df = oc.cosmos.get_pkn('human')
```
