# STITCH Client Comparison Report

## Executive Summary

The new STITCH client (`pypath/inputs/new_stitch/`) represents a significant architectural and functional improvement over the old implementation (`pypath/inputs/stitch.py`). The new version adopts a modular sub-package structure, provides richer data models, implements Cloudflare workarounds, and offers merged interaction data combining actions and links.

**Recommendation: Use new_stitch for COSMOS PKN integration.**

---

## Structural Comparison

### Old Implementation (`stitch.py`)

**Architecture:**
- Single file (~185 lines)
- Two independent functions
- Simple namedtuple records

**File Structure:**
```
inputs/
└── stitch.py
```

### New Implementation (`new_stitch/`)

**Architecture:**
- Modular sub-package
- Separation of concerns across 6 files
- Hierarchical data models with nested entities

**File Structure:**
```
inputs/new_stitch/
├── __init__.py          # Public API exports
├── _raw.py              # Low-level table download
├── _records.py          # Data structure definitions
├── _actions.py          # Actions processing
├── _links.py            # Links processing
└── _interactions.py     # Merged interactions
```

---

## Functional Comparison

### Available Functions

| Old Implementation | New Implementation | Notes |
|--------------------|-------------------|-------|
| `stitch_actions_interactions()` | `actions()` | Enhanced parsing |
| `stitch_links_interactions()` | `links()` | Simplified, no binding score |
| N/A | `interactions()` | **NEW**: Merged actions + links |
| N/A | `tables()` | **NEW**: Generic table reader |

---

## Key Differences

### 1. Data Models

#### Old Implementation

**StitchActionsInteraction:**
```python
namedtuple(
    'StitchActionsInteraction',
    ['partner_a', 'partner_b', 'mechanism', 'action', 'score']
)
```
- Flat structure
- IDs as strings (e.g., '1234')
- Action as string (e.g., 'activation', 'inhibition')

**StitchLinksInteraction:**
```python
namedtuple(
    'StitchLinksInteraction',
    ['partner_a', 'partner_b', 'experimental', 'prediction',
     'database', 'textmining', 'combined_score', 'physical_combined_score']
)
```
- Includes optional physical_combined_score from binding actions
- No entity metadata

#### New Implementation

**Entity (Nested Type):**
```python
namedtuple(
    'Entity',
    ['id', 'type', 'stereospecific', 'ncbi_tax_id']
)
```
- Rich metadata per partner
- Type discrimination: 'small_molecule' vs 'protein'
- Stereospecificity flag
- Taxonomy information

**StitchAction:**
```python
namedtuple(
    'StitchAction',
    ['source', 'target', 'directed', 'mode',
     'activation', 'inhibition', 'score']
)
```
- Source/target as Entity objects
- Boolean `directed` flag
- Separate `activation` and `inhibition` booleans (not string)
- Mode field (mechanism)

**StitchLinks:**
```python
namedtuple(
    'StitchLinks',
    ['chemical_id', 'protein_id', 'experimental', 'prediction',
     'database', 'textmining', 'combined_score',
     'ncbi_tax_id', 'stereospecific']
)
```
- Explicit chemical/protein roles
- Taxonomy and stereo flags included
- **No physical_combined_score** (removed feature)

**StitchInteractions (NEW):**
```python
namedtuple(
    'StitchInteractions',
    ['source', 'target', 'directed', 'mode',
     'activation', 'inhibition', 'experimental', 'prediction',
     'database', 'textmining', 'combined_score', 'final_score']
)
```
- Merges action metadata with link evidence scores
- Provides complete interaction picture in single record
- `final_score` from actions, `combined_score` from links

---

### 2. ID Parsing

#### Old Implementation

```python
sep = re.compile(r'[sm\.]')
a = sep.split(l[0])[1]  # Simple split
b = sep.split(l[1])[1]
```
- Basic split on delimiters
- Loses stereospecificity information
- Loses taxonomy information
- No distinction between CID/ENSP prefixes

#### New Implementation

**Actions:**
```python
REID = re.compile(r'((?:\d+)?)\.?(CID|ENSP)([ms]?)(0*)(\d+)')
# Extracts: (tax, ens, stereo, zeros, id)
```

**Links:**
```python
RECHEMID = re.compile(r'CID([ms]?)0*(\d+)')   # (stereo, id)
REPROTID = re.compile(r'(\d+)\.?(ENSP0*\d+)') # (tax, ensid)
```

**Improvements:**
- Preserves stereospecificity ('s' vs 'm' vs '')
- Extracts taxonomy IDs
- Handles leading zeros properly
- Distinguishes CID (PubChem) from ENSP (Ensembl Protein)

---

### 3. Download Mechanism

#### Old Implementation

```python
import pypath.share.curl as curl

c = curl.Curl(url, silent=False, large=True, slow=True)
for l in c.result:
    # process line
```
- Uses pypath's Curl wrapper
- Can be blocked by Cloudflare
- Integrated with pypath caching system

#### New Implementation

```python
import requests
import gzip
from io import BytesIO

response = requests.get(url, stream=True, timeout=300, allow_redirects=True)
with gzip.open(BytesIO(response.content), 'rt') as f:
    reader = csv.DictReader(f, delimiter='\t')
```

**Improvements:**
- **Cloudflare workaround**: Direct `requests` library bypasses blocking
- CSV DictReader for cleaner column access
- Explicit gzip handling
- Stream processing with configurable max_lines

**Trade-off:**
- Bypasses pypath's caching infrastructure (may need separate caching strategy)

---

### 4. Direction Handling

#### Old Implementation

```python
if l[4] == 'f':  # Direction flag
    a, b = b, a  # Swap partners
```
- Simple swap based on 'f' flag
- All interactions treated as directed after swap
- No distinction between bidirectional and unidirectional

#### New Implementation

```python
directed = action["a_is_acting"].lower() == 't'

# In actions() function:
if set(out[0][:2]) == set(out[1][:2]):  # Same partners
    if not any(act.directed for act in out):
        yield out[0]  # Undirected, yield once
    else:
        for act in out:
            if act.directed:
                yield act  # Directed, yield specific directions
```

**Improvements:**
- Explicit `directed` boolean flag
- Smart deduplication of bidirectional interactions
- Preserves directionality information in output
- Handles both directed and undirected cases correctly

---

### 5. Action/Mechanism Parsing

#### Old Implementation

```python
StitchActionsInteraction(
    mechanism = l[2],      # Raw string
    action = l[3] or None, # Raw string or None
    score = int(l[5]),
)
```
- Action as string: 'activation', 'inhibition', or empty
- Single mechanism field

#### New Implementation

```python
def parse_action(action: str) -> tuple[bool, bool]:
    if action == 'activation':
        return True, False
    if action == 'inhibition':
        return False, True
    else:
        return False, False

StitchAction(
    mode = action["mode"],        # Mechanism
    activation = activation,       # Boolean
    inhibition = inhibition,       # Boolean
    score = int(action["score"]),
)
```

**Improvements:**
- Separate boolean flags for activation/inhibition
- Easier to filter/query by effect type
- Supports programmatic logic (no string comparison needed)

---

### 6. Physical Binding Score Feature

#### Old Implementation

```python
if physical_interaction_score:
    phy_links = dict(
        ((s.partner_a, s.partner_b), s.score)
        for s in stitch_actions_interactions()
        if s.mechanism == 'binding'
    )
    # ... later lookup and include in output
    physical_combined_score = phy_score
```
- Optional feature to enrich links with binding scores
- Requires double download (actions + links)

#### New Implementation

**Not implemented** - physical_combined_score removed.

**Rationale:**
- The new `interactions()` function provides complete merged data
- Users can filter actions by `mode == 'binding'` if needed
- Cleaner separation of actions vs links

**Impact for COSMOS:**
- If physical binding scores are critical, may need custom implementation
- Can be reconstructed by filtering merged interactions

---

## Merged Interactions Feature (NEW)

The new implementation introduces `interactions()`, which combines actions and links:

```python
def interactions(max_lines=None, ncbi_tax_id=9606):
    """Merges actions and links tables into unified interactions"""

    # 1. Load actions
    the_actions = actions(max_lines, ncbi_tax_id)

    # 2. Create lookup from links
    link_lookup = {
        (link.chemical_id, link.protein_id, link.ncbi_tax_id, link.stereospecific): link
        for link in links(max_lines, ncbi_tax_id)
    }

    # 3. Join and yield
    for action in the_actions:
        key = get_action_identifier(action)
        if key in link_lookup:
            yield StitchInteractions(...)
```

**Benefits:**
- Single function call for complete data
- Direction, mechanism, effect from actions
- Evidence scores from links
- Pre-joined, ready for downstream analysis

**For COSMOS PKN:**
This is the **recommended function** to use, as it provides:
- Chemical-protein interactions with direction
- Activation/inhibition effects
- Evidence quality scores
- Taxonomy and stereospecificity metadata

---

## Comparison Summary Table

| Feature | Old | New | Winner |
|---------|-----|-----|--------|
| **Architecture** | Single file | Modular sub-package | **New** |
| **Code organization** | Mixed concerns | Separation of concerns | **New** |
| **Data models** | Flat namedtuples | Hierarchical with Entity | **New** |
| **ID parsing** | Basic split | Comprehensive regex | **New** |
| **Stereospecificity** | Lost | Preserved | **New** |
| **Taxonomy info** | Lost | Preserved | **New** |
| **Download** | Curl (Cloudflare issues) | Requests (workaround) | **New** |
| **Direction handling** | Simple swap | Smart deduplication | **New** |
| **Effect encoding** | String | Boolean flags | **New** |
| **Merged interactions** | Not available | Available | **New** |
| **Physical binding score** | Optional | Not implemented | **Old** |
| **Caching integration** | Yes (pypath cache) | No | **Old** |

---

## Recommendation for COSMOS PKN

### Use `new_stitch`

**Primary reasons:**

1. **Richer metadata**: Entity objects with type, taxonomy, stereospecificity
2. **Merged interactions**: Single function call for complete chemical-protein data
3. **Better directionality**: Proper handling of directed vs undirected interactions
4. **Activation/inhibition flags**: Easy filtering for specific effect types
5. **Cloudflare compatibility**: Avoids download blocking issues
6. **Future-proof**: Modular design allows easier extension

### Migration Notes

**What changes:**

- **Function names**:
  - `stitch_actions_interactions()` → `actions()`
  - `stitch_links_interactions()` → `interactions()` (use merged version)

- **Data access**:
  ```python
  # Old
  from pypath.inputs.stitch import stitch_links_interactions
  for rec in stitch_links_interactions(ncbi_tax_id=9606):
      partner_a = rec.partner_a  # string ID
      partner_b = rec.partner_b  # string ID

  # New
  from pypath.inputs.new_stitch import interactions
  for rec in interactions(ncbi_tax_id=9606):
      chemical = rec.source  # Entity object
      protein = rec.target   # Entity object
      chem_id = rec.source.id
      prot_id = rec.target.id
      is_activation = rec.activation
      is_directed = rec.directed
  ```

- **Loss of physical_combined_score**:
  - If needed, can filter `actions()` for `mode == 'binding'`
  - Or implement custom merger

**Caching considerations:**

- New implementation bypasses pypath cache
- For COSMOS PKN scripts, consider:
  - Local CSV caching of results
  - Infrequent re-downloads (STITCH updates quarterly)
  - Or wrap `tables()` with custom cache decorator

---

## Implementation for COSMOS PKN Script

Based on Phase 1 patterns (transporter.py, receptor_metabolite.py), the STITCH script should:

1. **Use `interactions()` function** from new_stitch for complete data
2. **Extract chemical-protein pairs** with direction and effect
3. **Add Gene__ and Metab__ prefixes** following COSMOS conventions
4. **Handle location/compartment** (STITCH doesn't provide subcellular location directly - may need to use protein location from UniProt as in other scripts)
5. **Filter by score threshold** (e.g., combined_score > 700 for high confidence)
6. **Create forward/reverse reactions** if needed (similar to transporter logic)
7. **Generate node mappings** for genes and metabolites

**Key differences from other sources:**

- STITCH provides **compound-protein interactions** (not transporter substrate)
- May or may not require reverse reactions (depends on COSMOS network design)
- Location information not inherent - would need to be added from protein annotations
- Stereospecificity flag could be valuable for metabolite disambiguation

---

## Conclusion

The **new_stitch** implementation is strongly recommended for COSMOS PKN integration. It provides:

- ✅ Richer, more structured data
- ✅ Better parsing and metadata extraction
- ✅ Merged interactions combining actions and links
- ✅ Cloudflare-proof downloads
- ✅ Cleaner, more maintainable code architecture

**Next steps:**

1. Write STITCH COSMOS PKN build script using `new_stitch.interactions()`
2. Determine if location annotation is needed (integrate with UniProt locations)
3. Define score thresholds for COSMOS quality requirements
4. Decide on forward/reverse reaction logic based on COSMOS network semantics
