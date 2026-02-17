# Stereospecificity in Chemical-Protein Interactions

## What is Stereospecificity?

**Stereospecificity** refers to whether the 3D spatial arrangement of atoms in a molecule (its stereochemistry) affects its biological activity and interactions.

## The Chemistry Basics

### Stereoisomers

**Stereoisomers** are molecules with:
- ✅ Same chemical formula
- ✅ Same atom-to-atom connections
- ❌ Different 3D spatial arrangements

Think of it like your left and right hands:
- Same parts (fingers, palm, thumb)
- Same connections
- But mirror images - can't be superimposed
- A left glove only fits a left hand!

### Types of Stereoisomers

**1. Enantiomers (Mirror Images)**
```
      COOH                COOH
       |                    |
   H--C--NH2          H2N--C--H
       |                    |
      CH3                  CH3
   L-Alanine           D-Alanine
   (mirror images)
```

**2. Diastereomers (Non-mirror stereoisomers)**
- Multiple chiral centers
- Not mirror images but still different 3D arrangements

## Biological Importance

### Example 1: Glucose
- **D-Glucose**: The form used by all living cells for energy
- **L-Glucose**: Mirror image, tastes sweet but cells can't metabolize it
- **Formula**: Both C₆H₁₂O₆
- **Effect**: Completely different biological function!

### Example 2: Drugs

**Thalidomide (tragic historical example)**
- **(R)-Thalidomide**: Effective sedative/anti-nausea
- **(S)-Thalidomide**: Caused severe birth defects
- Given as mixture in 1950s-60s → tragedy
- Led to strict stereoisomer regulations in pharmaceuticals

**Ibuprofen (common painkiller)**
- **(S)-Ibuprofen**: Active anti-inflammatory
- **(R)-Ibuprofen**: Inactive (slowly converts to S-form in body)
- Sold as racemic mixture (both forms)

### Example 3: Taste and Smell

**Carvone**
- **(R)-Carvone**: Smells like spearmint
- **(S)-Carvone**: Smells like caraway
- Same molecule, different smell!

**Limonene**
- **(R)-Limonene**: Smells like oranges
- **(S)-Limonene**: Smells like lemons

## In Protein-Chemical Interactions

### Why Stereospecificity Matters

Proteins have **chiral binding pockets** (like the left-hand glove):
- Active sites are 3D structures
- Only specific stereoisomers fit properly
- Wrong stereoisomer → no binding or different effect

**Example: Enzyme-Substrate Recognition**
```
Enzyme active site (chiral pocket)
        ___
       /   \
      |  🔑 |  ← Only R-enantiomer fits
       \___/

        ___
       /   \
      | 🔑  |  ← S-enantiomer doesn't fit
       \___/
```

### Real Biological Examples

**1. Amino Acids**
- Proteins only use L-amino acids
- D-amino acids won't be incorporated into proteins
- Some D-amino acids used in cell walls (bacteria)

**2. Sugars**
- Cells primarily use D-sugars
- Enzymes won't recognize L-sugars

**3. Neurotransmitters**
- L-DOPA: Parkinson's disease treatment (crosses blood-brain barrier)
- D-DOPA: Inactive for Parkinson's

## In STITCH Database

### STITCH Stereospecificity Flags

When you see in STITCH data:

**`stereospecific = True`**
- The database is tracking a **specific stereoisomer**
- The PubChem CID refers to a particular 3D form
- Different stereoisomers would have different CIDs
- The interaction is specific to this stereoisomer

**`stereospecific = False`**
- The compound has no stereoisomers (achiral/symmetric)
- OR: Database refers to compound generically (all stereoisomers)
- OR: Stereochemistry doesn't affect this interaction
- Multiple stereoisomers share the same identifier

### STITCH ID Notation

In raw STITCH files, you might see:
- **`CIDs0001234`**: `s` = stereospecific variant
- **`CIDm0001234`**: `m` = mixture (racemic, all stereoisomers)
- **`CID0001234`**: No flag (depends on compound nature)

### Example from Our Data

From the sample:
```
chemical_id = 23627457
stereospecific = True
```
This means PubChem CID 23627457 refers to a **specific stereoisomer**, and the protein interaction is specific to this particular 3D form.

## Implications for COSMOS PKN

### Should We Distinguish Stereoisomers?

**Option 1: Keep Stereo Distinct** (stereospecific = True → separate nodes)
```
Metab__23627457_s_c  (stereospecific, cytosol)
Metab__23627457_m_c  (mixture/other form, cytosol)
```
✅ Pros:
- Chemically accurate
- Preserves biological specificity
- Important for drugs/metabolites with known stereo effects

❌ Cons:
- More complex network
- Many metabolites have unknown stereo in practice
- Potential confusion if stereo not consistently tracked

**Option 2: Merge Stereo** (ignore stereospecific flag)
```
Metab__23627457_c  (all forms)
```
✅ Pros:
- Simpler network
- Matches how many databases handle compounds
- Easier integration with other data sources

❌ Cons:
- Loses chemical precision
- May miss real biological differences

### Recommendation for COSMOS

**Suggested approach:**
1. **Store the stereospecific flag** in node attributes
2. **Don't create separate nodes** by default
3. **Allow filtering** to create stereo-specific networks if needed

Example:
```python
# Single node ID
metabolite_id = 'Metab__23627457_c'

# But store attribute
metabolite_attrs = {
    'pubchem_cid': '23627457',
    'stereospecific': True,
    'compartment': 'c'
}
```

This gives flexibility:
- Default: merged network (simpler)
- Advanced: can split by stereospecific flag for detailed analysis

## Further Reading

**Chemistry:**
- Stereoisomers vs structural isomers
- Chiral centers (carbon with 4 different groups)
- R/S nomenclature (Cahn-Ingold-Prelog rules)
- Enantiomers vs diastereomers

**Pharmacology:**
- Racemic drugs vs enantiopure drugs
- FDA regulations on stereoisomers (since 1992)
- Chiral switching (developing enantiopure version of racemic drug)

**Databases:**
- PubChem CID vs PubChem SID (substance vs compound)
- ChEBI stereo annotation
- InChI strings (encode stereochemistry)

## Summary

**Stereospecificity = "Does 3D shape matter?"**

- **YES (stereospecific = True)**: Specific 3D form required for biological activity
- **NO (stereospecific = False)**: Generic compound or 3D shape doesn't matter

**In COSMOS PKN context:**
- STITCH preserves this information
- We can choose to use it or not
- Recommended: store as attribute, don't split nodes by default
- Allows flexibility for future stereo-specific analysis if needed
