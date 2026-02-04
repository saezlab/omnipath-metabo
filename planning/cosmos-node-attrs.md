Random edges from the ready COSMOS PKN in its original format:

| From | To |
| --- | --- |
| Metab__CHEBI:9251_l | Gene1__O95210 |
| Gene1__O95210 | Metab__CHEBI:9251_c |
| Metab__CHEBI:9251_c | Gene1__O95210_rev |
| Gene1__O95210_rev | Metab__CHEBI:9251_l |
| O95210 | Gene1__O95210 |
| O95210 | Gene1__O95210_rev |
| CHEBI:9251 | Metab__CHEBI:9251_l |


- `_rev`: reversed direction of reaction (avoid graph loops)
- `Metab__`: metabolite
- `Gene1__`: gene, number is unique identifier to have the same gene as a
    separate node in different reactions
- `_c`: cytoplasm (compartment), others are `_n` nucleus, `_m` mitochondria,
    etc.
