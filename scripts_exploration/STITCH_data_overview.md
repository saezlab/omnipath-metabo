# STITCH new_stitch.interactions() Data Overview

## Summary of Output

The `new_stitch.interactions()` function successfully merges STITCH actions and links data, providing comprehensive chemical-protein interaction information.

## Data Structure

### Record Format: `StitchInteractions`

Each interaction record contains:

```python
StitchInteractions(
    source = Entity(id, type, stereospecific, ncbi_tax_id),
    target = Entity(id, type, stereospecific, ncbi_tax_id),
    directed = bool,
    mode = str,
    activation = bool,
    inhibition = bool,
    experimental = int,
    prediction = int,
    database = int,
    textmining = int,
    combined_score = int,
    final_score = int,
)
```

## Sample Data Analysis (100 interactions)

### Interaction Types

**Source/Target Combinations:**
- **Protein → Small Molecule**: 94% (e.g., protein binding to compound)
- **Small Molecule → Protein**: 6% (e.g., compound affecting protein expression)

### Modes of Interaction

| Mode | Count | Description |
|------|-------|-------------|
| `binding` | 86 | Physical binding interaction |
| `pred_bind` | 8 | Predicted binding |
| `expression` | 5 | Affects gene/protein expression |
| `reaction` | 1 | Catalytic/chemical reaction |

### Direction

- **Undirected**: 95% (bidirectional or non-directional binding)
- **Directed**: 5% (directional effect like expression regulation)

### Activation/Inhibition

**In this sample of 100 interactions:**
- Activation: 0
- Inhibition: 0
- Neither: 100

**Note**: The first 100 interactions don't contain activation/inhibition examples. These effects are present in the full dataset but are less common than binding interactions.

### Stereospecificity

**Target (small molecules):**
- Stereospecific: 74%
- Non-stereospecific: 26%

**Source (small molecules):**
- Only 6 small molecules as source
- 50% stereospecific

### Evidence Scores

**Final Score Statistics** (action score):
- Mean: 333
- Median: 258
- Range: 150-951
- Q1: 207, Q3: 381

**Evidence Types** (from links):
- **Experimental**: Most interactions have experimental evidence (0-884)
- **Database**: Score of 150 for expression-type interactions
- **Prediction**: Used for pred_bind mode
- **Textmining**: Mostly 0 in this sample

## Example Interactions

### Example 1: Expression Regulation (Directed)
```
Source: 10461 (small_molecule)
Target: ENSP00000170630 (protein)
Mode: expression
Directed: True
Combined Score: 150 (database evidence)
```
Compound 10461 affects expression of protein ENSP00000170630.

### Example 2: Protein-Compound Binding (Undirected)
```
Source: ENSP00000353915 (protein)
Target: 23627457 (small_molecule, stereospecific)
Mode: binding
Directed: False
Experimental: 191
Combined Score: 191
```
Protein binds to stereospecific compound with experimental evidence.

### Example 3: Predicted Binding
```
Source: ENSP00000267377 (protein)
Target: 23590374 (small_molecule, stereospecific)
Mode: pred_bind
Directed: False
Experimental: 159
Prediction: 170
Combined Score: 282
Final Score: 170
```
Predicted binding with both experimental and computational evidence.

### Example 4: High-Confidence Binding
```
Source: ENSP00000396308 (protein)
Target: 467905 (small_molecule, stereospecific)
Mode: binding
Experimental: 884
Combined Score: 884
Final Score: 884
```
Strong experimental evidence for binding.

## ID Types

### Chemical IDs (PubChem CIDs)
- Examples: `10461`, `2336`, `2361`, `1567`, `467905`
- Can be stereospecific or non-stereospecific
- Range: 4-8 digit numbers

### Protein IDs (Ensembl Protein)
- Format: `ENSP` + 11 digits (e.g., `ENSP00000170630`)
- Always from specified organism (ncbi_tax_id = 9606 for human)

## Key Observations for COSMOS PKN Integration

### 1. Source/Target Orientation
- **Not consistent**: Sometimes chemical is source, sometimes target
- Need to normalize: Chemical → Protein for consistency
- Check `source.type` and `target.type` to identify roles

### 2. Direction Handling
- Most interactions are undirected (binding)
- Directed interactions (expression) are minority
- May want to filter by `directed` flag depending on network semantics

### 3. Mode Filtering
- **Binding interactions** dominate (86%)
- May want to:
  - Include only `mode == 'binding'` for physical interactions
  - Or include `pred_bind` for broader coverage
  - Or include `expression` for regulatory networks

### 4. Score Thresholds
Common thresholds:
- **150**: Minimum (many low-confidence interactions)
- **400**: Medium confidence
- **700**: High confidence
- **900**: Highest confidence

**Recommendation**: Use 400-700 threshold for COSMOS to balance coverage and quality.

### 5. Activation/Inhibition
- Present in full dataset but rare in early records
- When present, encoded as boolean flags
- Useful for regulatory networks
- Can assign `mor` (mode of regulation):
  - activation → mor = 1
  - inhibition → mor = -1
  - neither → mor = 0

### 6. Stereospecificity
- Majority of small molecules are stereospecific
- Flag preserved in Entity objects
- Could use to create distinct metabolite nodes if needed
- Example: `23627457` (stereo=True) vs `23627457` (stereo=False) as different entities

## Comparison with Phase 1 Sources

### Similar to: Transporter (TCDB/SLC)
- Chemical-protein interactions
- Need Gene__ and Metab__ prefixes
- Location annotation would be beneficial (not provided by STITCH)

### Different from: Transporter
- STITCH: Binding and regulatory interactions (not necessarily transport)
- Directionality less clear (binding is typically undirected)
- No inherent subcellular location information
- Broader interaction types (binding, expression, catalysis)

### Similar to: Receptor-Metabolite (MRCLinksDB)
- Chemical-protein interactions with direction
- Regulatory implications

## Recommendations for COSMOS PKN Script

1. **Normalize orientation**: Always put chemical as source, protein as target
   ```python
   if rec.source.type == 'small_molecule':
       chemical_id = rec.source.id
       protein_id = rec.target.id
   else:
       chemical_id = rec.target.id
       protein_id = rec.source.id
   ```

2. **Score filtering**: Use combined_score >= 400 or 700

3. **Mode filtering**: Decide which modes to include
   - Binding only: most conservative
   - Binding + pred_bind: broader coverage
   - All modes: comprehensive but may include less relevant interactions

4. **Prefixing**:
   - Chemical: `Metab__[PubChem_CID]`
   - Protein: `Gene__[ENSP_ID]`

5. **Location**: Consider adding compartment suffix
   - STITCH doesn't provide location
   - Could integrate with UniProt protein location data
   - Or use default compartment (e.g., `_c` for cytosol)

6. **Reverse reactions**: Probably NOT needed
   - Binding interactions are already undirected
   - Expression interactions are directional effects, not transport

7. **Stereospecificity**: Decide if needed
   - Could append `_s` or `_m` suffix to metabolite IDs
   - Or ignore and use generic CID

## Data Availability

- **Download size**: ~GB (gzipped TSV files)
- **Number of interactions**: Tens of thousands (full dataset)
- **Download time**: Few minutes depending on connection
- **Caching**: Not cached by pypath (uses requests directly)

## Next Steps

1. Decide on filtering criteria (score threshold, modes, direction)
2. Decide on location handling (integrate UniProt or use default)
3. Decide on stereospecificity handling
4. Test script with various thresholds to assess data quality/quantity
5. Integrate into COSMOS PKN builder
