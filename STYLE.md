# Style Guide

Coding style guidelines for the omnipath-metabo package.

## Function Naming

Avoid redundant prefixes in function names. Common prefixes to avoid:

- `get_` - omit when the function simply returns data
- `load_` - omit when the function reads from a file/source
- `parse_` - omit when the function transforms input data
- `create_` - omit when the function constructs an object
- `add_` - consider omitting when the action is implied by context

### Examples

```python
# Avoid
def get_uniprot_locations(): ...
def load_tcdb_locations(): ...
def parse_equation(): ...
def create_node_mappings(): ...

# Prefer
def uniprot_locations(): ...
def tcdb_locations(): ...
def equation(): ...
def node_mappings(): ...
```

Exceptions: Keep prefixes when they disambiguate or add meaning:

```python
# Keep: 'strip' describes the transformation
def strip_prefix(): ...

# Keep: distinguishes from related function
def add_reverse(): ...  # if there's also remove_reverse()
```

## Imports

Use lazy imports for heavy dependencies (like pypath) inside functions
rather than at module level. This allows importing lighter parts of the
package without triggering the full dependency chain.

```python
# Avoid at module level
from pypath.inputs.uniprot import uniprot_locations

# Prefer inside functions
def uniprot_locations(organism: int = 9606) -> dict:
    from pypath.inputs.uniprot import uniprot_locations as _uniprot_locations
    return _uniprot_locations(organism=organism)
```

## Data Structures

Return the most appropriate data structure for the use case:

- Use `dict` when data will be used for lookups
- Use `list[tuple]` for pairs/relationships
- Reserve `DataFrame` for tabular data that benefits from pandas operations

```python
# Return dict for lookup tables
def tcdb_locations() -> dict[str, str]:
    df = pd.read_csv(data_path('tcdb_locations.csv'))
    return {loc: abbrev for loc, abbrev in zip(df['location'], df['abbreviation'])}

# Return list of tuples for routes/pairs
def tcdb_routes() -> list[tuple[str, str]]:
    df = pd.read_csv(data_path('tcdb_routes.csv'))
    return [(row['from'], row['to']) for _, row in df.iterrows()]
```
