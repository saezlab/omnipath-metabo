# COSMOS PKN — ID Translation Gaps

Observed during `build_transporters()` run (2026-03-14, human / 9606).
These rows are silently dropped after translation fails. Fixing them increases
PKN coverage and should be prioritised before a production release.

---

## Current drop counts

### Source ID failures (metabolite → ChEBI or protein → UniProt)

| Resource | id_type | Rows dropped |
|---|---|---|
| STITCH | pubchem | 23,845 |
| GEM:Human-GEM | metatlas | 15,707 |
| GEM_transporter:Human-GEM | metatlas | 3,257 |
| Recon3D | bigg | 5,379 |
| GEM:Human-GEM | ensembl | 863 |
| Recon3D | entrez | 211 |
| GEM_transporter:Human-GEM | ensembl | 326 |
| **Total** | | **49,588** |

### Target ID failures (protein → UniProt or metabolite → ChEBI)

| Resource | id_type | Rows dropped |
|---|---|---|
| GEM:Human-GEM | metatlas | 16,768 |
| STITCH | ensp | 6,126 |
| Recon3D | bigg | 5,379 |
| GEM_transporter:Human-GEM | metatlas | 3,258 |
| GEM:Human-GEM | ensembl | 817 |
| GEM_transporter:Human-GEM | ensembl | 326 |
| Recon3D | entrez | 211 |
| **Total** | | **32,885** |

---

## Root cause analysis per id_type

### `pubchem` → ChEBI (STITCH, 23,845 rows)

UniChem PubChem→ChEBI bulk mapping does not cover all PubChem CIDs.
STITCH includes many low-confidence or non-drug-like compounds that
have no ChEBI entry. Also includes stereospecific CIDs (s-prefix in
STITCH) that may not map 1:1.

**Possible fixes:**
- Use PubChem→ChEBI via PubChemPy compound annotations (richer than UniChem)
- Accept HMDB as fallback ID when ChEBI is unavailable
- Filter STITCH to drug-like / bioactive compounds before translation

### `metatlas` → ChEBI (Human-GEM, ~19,000 + ~6,500 rows)

The MetAtlas TSV `metabolites.tsv` only has `metChEBIID` populated for a
subset of metabolites. Many GEM metabolites have no ChEBI annotation
in the MetAtlas source file.

**Possible fixes:**
- Add MetaNetX bridge (MNX → ChEBI) for metabolites missing direct ChEBI;
  a partial implementation already exists in `_bigg_to_chebi()` — apply
  the same strategy to `metatlas` metabolites
- Use the BiGG metabolite ID embedded in MetAtlas entries as a secondary
  lookup via `_bigg_to_chebi()`
- Cross-reference via HMDB or KEGG IDs also present in the MetAtlas TSV

### `bigg` → ChEBI (Recon3D, 5,379 rows)

`_bigg_to_chebi()` covers ~28% of Recon3D metabolites directly, with a
MetaNetX bridge adding ~2.7 pp more. The remaining ~70% have no ChEBI
cross-reference in the BiGG JSON or MetaNetX.

**Possible fixes:**
- Add KEGG→ChEBI bridge: BiGG JSON includes `kegg.compound` annotations;
  use UniChem or KEGG API to get ChEBI
- Add HMDB→ChEBI bridge: BiGG JSON includes `hmdb` annotations; already
  have `_hmdb_to_chebi()` — apply to BiGG metabolites
- Combine all cross-reference types in a priority chain

### `ensembl` → UniProt (Human-GEM, ~1,200 rows combined)

Some ENSG IDs in Human-GEM fail `ensg → uniprot` BioMart mapping.
Likely causes: retired Ensembl IDs, readthrough/pseudogene loci, or
IDs from an older Ensembl release than the current BioMart.

**Possible fixes:**
- Use Ensembl REST API as fallback for IDs missing from BioMart
- Pin the GEM Ensembl release version and use a matching BioMart archive
- Accept gene symbol as tertiary fallback (MetAtlas TSV contains gene names)

### `entrez` → UniProt (Recon3D, 211 rows)

Some Entrez Gene IDs in Recon3D BiGG fail both the BiGG gene-symbol
path and the pypath `ncbigene → uniprot` BioMart fallback.
Likely retired gene IDs or non-human contamination.

**Possible fixes:**
- Use NCBI Gene API to resolve retired/merged Entrez IDs
- Cross-reference via gene symbol embedded in BiGG JSON

### `ensp` → UniProt (STITCH, 6,126 rows)

Some Ensembl protein IDs in the STITCH actions file fail
`ensp → uniprot` BioMart mapping. Likely retired ENSP IDs from
older STRING/STITCH releases (STITCH 5 uses Ensembl 89).

**Possible fixes:**
- Use Ensembl REST API `lookup/id` endpoint for retired ENSP IDs
- Map via gene symbol: STITCH proteins can be cross-referenced to
  gene symbols via the STRING protein info file

---

## TCDB TC class filtering

Currently `tcdb_classes()` is used directly to mark all TCDB entries as
`'transporter'`. However, TCDB covers more than just canonical membrane
transporters — TC classes 4 (group translocators / kinases), 8 (accessory
factors), and 9 (incompletely characterized) include proteins like SRC
(P12931) and IL-1β (P01584) that are not transporters in the COSMOS sense.

**Fix:** Filter `tcdb_classes()` to only keep TC classes 1, 2, and 3:

```python
for uniprot, (tc, family) in tcdb_classes().items():
    if tc[0] in ('1', '2', '3'):
        result[uniprot] = 'transporter'
```

Also, the original `tcdb_annotations(organism)` validated UniProt ACs via
`mapping.map_name` and `reflists.check` for the target organism. By using
`tcdb_classes()` directly we lost organism filtering. Consider restoring it.

---

## Priority order

1. **metatlas → ChEBI** — largest absolute loss (>22k rows); MetaNetX
   and BiGG bridges are low-effort extensions of existing code
2. **pubchem → ChEBI (STITCH)** — 23k rows but many are low-interest
   compounds; prioritise drug-like compound coverage first
3. **bigg → ChEBI (Recon3D)** — KEGG + HMDB bridges are straightforward
4. **ensp → UniProt (STITCH)** — Ensembl REST fallback needed
5. **ensembl → UniProt (GEM)** — smaller loss; BioMart archive pinning
6. **entrez → UniProt (Recon3D)** — smallest loss; low priority
