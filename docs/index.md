# omnipath-metabo

**Integrated prior-knowledge for metabolomics and multi-omics analysis.**

omnipath-metabo builds and serves the COSMOS prior-knowledge network (PKN)
-- a multi-layer directed signed network for causal reasoning across
metabolomics, proteomics, and transcriptomics data with
[cosmosR](https://github.com/saezlab/cosmosR).

## What is omnipath-metabo?

omnipath-metabo integrates curated databases of metabolite-protein
interactions, protein-protein signaling, and gene regulation into a
unified network. It covers:

- **Metabolite transport** -- membrane transporters from TCDB, SLC,
  genome-scale metabolic models (GEMs), and MRCLinksDB
- **Receptor binding** -- metabolite-receptor interactions from STITCH
  and MRCLinksDB
- **Allosteric regulation** -- enzyme modulation by metabolites from
  BRENDA and STITCH
- **Enzyme-metabolite reactions** -- stoichiometric edges from GEMs
- **Protein-protein signaling** -- kinase-substrate and other signaling
  from OmniPath
- **Gene regulation** -- transcription factor targets from CollecTRI
  and DoRothEA

## Key features

- **7 resources** -- TCDB, SLC, GEMs, MRCLinksDB, STITCH, BRENDA, OmniPath
- **6 network layers** -- transporters, receptors, allosteric,
  enzyme-metabolite, PPI, GRN
- **Multi-organism** -- human, mouse, rat, and more via orthology
  translation
- **Pre-built PKNs** -- ready-to-use networks served via the web service
- **Web service** -- REST API at
  [metabo.omnipathdb.org](https://metabo.omnipathdb.org)
- **AnnNet support** -- convert to
  [AnnNet](https://github.com/saezlab/annnet) graph objects
- **CLI** -- `cosmos-pkn` command for scripted exports

## Quick example

### Build the transporter layer locally

```python
from omnipath_metabo.datasets import cosmos

transporters = cosmos.build_transporters()
```

### Fetch the full PKN via the client

```python
import omnipath_client as oc

df = oc.cosmos.get_pkn('human')
```

## Learn more

- **[Quickstart](quickstart.md)** -- build, format, and serve a PKN
- **[COSMOS PKN vignette](vignettes/cosmos-pkn.md)** -- detailed walkthrough
  of all categories, organisms, and options
- **[API Reference](reference/index.md)** -- full function documentation
- **[Installation](installation.md)** -- setup instructions

## Services

| Service | URL | What it provides |
|---------|-----|-----------------|
| OmniPath Metabo | [metabo.omnipathdb.org](https://metabo.omnipathdb.org) | COSMOS PKN, metabolite-protein interactions |
| OmniPath Utils | [utils.omnipathdb.org](https://utils.omnipathdb.org) | ID translation, taxonomy, orthology |
| OmniPath Database | [dev.omnipathdb.org](https://dev.omnipathdb.org) | Interactions, annotations, complexes, ontology |
