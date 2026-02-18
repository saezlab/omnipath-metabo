# OmniPath Metabo Development Plan

This document outlines the development strategy for the `omnipath-metabo` package,
a Python library for metabolite, compound, reaction, and interaction-related
prior-knowledge database access, processing, and web service functionality.

## Project Overview

The package is part of the OmniPath ecosystem for systems biology and biomedicine
data processing. It will serve as a server-side component with optional client-side
usage capability. Licensed under BSD 3-Clause for industry compatibility.

### Key Objectives

- Compile specialized prior-knowledge networks (e.g., COSMOS PKN) from multiple
  data sources
- Process and transform metabolite/compound/reaction data for network
  optimization and multi-omics integration
- Eventually integrate with the OmniPath PostgreSQL database and web APIs

---

## Long-Term Development Phases

### Phase 1: Foundation and COSMOS PKN Integration (Current)

Establish package structure and incorporate existing COSMOS PKN build scripts.

### Phase 2: Resource Client Expansion

Develop and integrate clients for additional data sources:
- Reactome (reactions and pathways)
- STITCH (compound-protein interactions)
- SwissLipids (lipid hierarchy)
- Genome-scale metabolism models (GEMs)
- Additional sources as identified

### Phase 3: Processing Algorithms

Implement cheminformatics and network processing algorithms:
- RDKit integration for structure processing
- Network assembly and transformation utilities
- ID mapping and cross-referencing tools

### Phase 4: Database Integration

Connect to the OmniPath PostgreSQL database:
- Schema definitions for metabolite/compound entities
- Query builders for complex data retrieval
- Caching strategies for performance

### Phase 5: Web API Development

Build web service layer:
- REST API endpoints for PKN queries
- Data upload endpoints for custom processing
- Combined query result operations

---

## Phase 1: Detailed Plan

### 1.1 Package Instantiation

Create fresh package from the saezlab/python-project cookiecutter template:

- Package name: `omnipath-metabo`
- Module name: `omnipath_metabo`
- Use modern tooling: uv, ruff, pytest, hatchling, mkdocs-material
- BSD 3-Clause license
- GitHub workflows for CI/CD

### 1.2 Module Structure

```
omnipath_metabo/
├── __init__.py
├── _metadata.py
├── _session.py              # Logging and session management
├── _constants.py            # Package-wide constants
│
├── datasets/                # Curated, branded prior-knowledge datasets
│   ├── __init__.py
│   │
│   └── cosmos/              # COSMOS PKN - all build logic contained here
│       ├── __init__.py
│       ├── _build.py        # Main PKN builder, combines all sources
│       ├── _base.py         # Base classes for COSMOS components
│       │
│       ├── sources/         # Data source processors (temporary pypath integration)
│       │   ├── __init__.py
│       │   ├── tcdb.py      # TCDB transporter processing
│       │   ├── slc.py       # SLC-table transporter processing
│       │   ├── brenda.py    # BRENDA allosteric regulation
│       │   ├── mrclinksdb.py # MRCLinksDB receptor-metabolite
│       │   └── rhea.py      # Rhea reactions
│       │
│       ├── network.py       # Edge/node creation, direction handling
│       ├── location.py      # Subcellular location parsing and mapping
│       │
│       └── data/            # COSMOS-specific data files
│           ├── __init__.py
│           └── location_abbreviations.csv
│
└── data/                    # Package-wide data files (if any)
    └── __init__.py
```

The `datasets/` module will host other branded prior-knowledge compilations in the
future (e.g., tissue-specific networks, disease-focused datasets). Each dataset
gets its own submodule containing all necessary build logic and temporary solutions.


### 1.3 Incorporate COSMOS PKN Scripts

Refactor scripts from `scripts_exploration/` into the `datasets/cosmos/` submodule:

| Original Script | Target Module | Key Components |
|-----------------|---------------|----------------|
| `transporter.py` | `cosmos/sources/tcdb.py`, `cosmos/sources/slc.py` | TCDB processing, SLC interactions |
| `allosteric_regulation.py` | `cosmos/sources/brenda.py` | BRENDA allosteric data extraction |
| `receptor_metabolite.py` | `cosmos/sources/mrclinksdb.py` | MRCLinksDB interaction processing |
| `lipidnet_rhea.ipynb` | `cosmos/sources/rhea.py` | Rhea equation parsing, edge creation |

Common functionality extracted to:
- `cosmos/location.py`: Location parsing and abbreviation mapping
- `cosmos/network.py`: Node prefixing, edge creation, direction handling

### 1.4 Refactoring Guidelines

#### Code Organization

- Keep all COSMOS-specific code contained within `datasets/cosmos/`
- Extract common patterns within COSMOS into shared modules (`location.py`, `network.py`)
- Separate data fetching from data transformation
- Use generators for memory-efficient iteration
- Create named tuples or dataclasses for structured records

#### Specific Improvements

1. **Location Parsing** (`cosmos/location.py`)
   - Consolidate duplicated `parse_uniprot_locations()` functions
   - Create reusable location abbreviation mapping loader
   - Handle location edge cases consistently

2. **Node Prefixing** (`cosmos/network.py`)
   - Standardize `Gene__` and `Metab__` prefixing
   - Create utility functions for prefix addition/removal
   - Support configurable prefix schemes

3. **Reaction Direction** (`cosmos/network.py`)
   - Extract reverse reaction logic from multiple scripts
   - Support bidirectional, LR, and RL directions
   - Handle edge pair creation consistently

4. **Species Filtering** (`cosmos/_base.py`)
   - Create base class with NCBI taxonomy ID handling
   - Reuse pypath reflists for species protein sets (temporary)
   - Support multi-species queries where applicable

#### Code Quality Standards

Follow `planning/coding-style.md`:

- PascalCase for classes, snake_case for functions/variables
- Resource names as single words (e.g., `TcdbProcessor`, `tcdb_substrates`)
- Two blank lines between functions/classes
- Blank line after opening and before closing blocks
- Type hints for all public functions
- Napoleon (Google) style docstrings
- Single quotes for strings

### 1.5 Public API Design

```python
# High-level COSMOS PKN builder
from omnipath_metabo.datasets import cosmos

pkn = cosmos.build(
    sources = ['tcdb', 'slc', 'brenda', 'mrclinksdb', 'rhea'],
    ncbi_tax_id = 9606,
    include_reverse = True,
)

# Access individual COSMOS sources (if needed)
from omnipath_metabo.datasets.cosmos import sources

tcdb_edges = sources.tcdb_interactions(ncbi_tax_id = 9606)
brenda_reg = sources.brenda_regulations(organisms = ['human'])

# Future branded datasets will follow the same pattern
# from omnipath_metabo.datasets import lipidnet
# from omnipath_metabo.datasets import metabolic_atlas
```

### 1.6 Testing Plan

#### Unit Tests

Create tests in `tests/` directory:

```
tests/
├── conftest.py              # Shared fixtures
└── test_datasets/
    └── test_cosmos/
        ├── conftest.py      # COSMOS-specific fixtures
        ├── test_build.py    # COSMOS PKN builder tests
        ├── test_location.py # Location parsing tests
        ├── test_network.py  # Edge/node creation tests
        └── test_sources/
            ├── test_tcdb.py     # TCDB processor tests
            ├── test_slc.py      # SLC processor tests
            ├── test_brenda.py   # BRENDA tests
            ├── test_mrclinksdb.py # MRCLinksDB tests
            └── test_rhea.py     # Rhea tests
```

#### Test Categories

1. **Parsing Tests**: Verify location string parsing, equation parsing
2. **Transformation Tests**: Test node prefixing, edge creation, direction handling
3. **Integration Tests**: Test source data processing with mock data
4. **End-to-End Tests**: Test full PKN build with sample data

#### Test Data

- Create small fixture files with representative data samples
- Use pytest fixtures for common test data
- Mock network requests for reproducible tests

### 1.7 Dependencies

Core dependencies:
- `pandas`: Data manipulation
- `pypath-omnipath`: Resource access (temporary, for pypath.inputs)
- `pypath-common`: Common utilities

Optional dependencies:
- `rdkit`: Cheminformatics (Phase 3)

Development dependencies:
- `pytest`, `pytest-cov`: Testing
- `ruff`: Linting and formatting
- `mkdocs-material`: Documentation

### 1.8 Implementation Order

1. **Initialize package** from cookiecutter template
2. **Create base infrastructure**: session, constants
3. **Create `datasets/cosmos/` submodule structure**
4. **Implement COSMOS utilities**: `location.py`, `network.py`
5. **Port transporter sources**: `sources/tcdb.py`, `sources/slc.py`
6. **Port allosteric source**: `sources/brenda.py`
7. **Port receptor source**: `sources/mrclinksdb.py`
8. **Port reaction source**: `sources/rhea.py`
9. **Create COSMOS PKN builder**: `_build.py` integrating all sources
10. **Write unit tests** for COSMOS modules
11. **Documentation**: API reference, tutorials

---

## Coding Conventions Summary

| Aspect | Convention |
|--------|------------|
| Class names | `PascalCase`, resources as single word (`TcdbProcessor`) |
| Functions/variables | `snake_case`, resources as single word (`tcdb_interactions`) |
| Blank lines | 2 between classes/functions, 1 around blocks |
| Indentation | 4 spaces, hanging indent for multi-line arguments |
| Quotes | Single quotes preferred |
| Docstrings | Napoleon (Google) style with type hints |
| Imports | Grouped: stdlib, typing, third-party, package |
| Module naming | Branded datasets as single lowercase word (`cosmos`, `lipidnet`) |

---

## Architecture Repository (saezverse)

A central coordination repository for cross-package architecture and design
decisions is available at `/home/denes/arch` (remote: saezlab/arch).

**When to use it:**
- Before making architectural choices that affect multiple packages, check
  `human/decisions/` for relevant ADRs
- When planning cross-repo work (e.g., pypath + omnipath-metabo), create a
  planning document in `ai/` following the `YYYY-MM-topic.md` naming convention
- Check `human/plans/` for roadmaps and specifications relevant to the current task
- Check `human/packages/` for descriptions of how packages relate to each other
- Check `human/guidelines/` for coding style conventions

**Structure:**
- `human/` — authoritative, human-reviewed documents (may edit for formatting/clarity)
- `ai/` — AI-generated working documents (create freely with YAML front matter)
- `skills/` — reusable instruction sets for AI assistants
- See `AGENTS.md` in that repo for full instructions

**Key documents for this package:**
- `human/plans/omnipath-metabo-cosmos.md` — COSMOS PKN build plan
- `human/plans/omnipath-metabo.md` — overall package roadmap
- `human/packages/omnipath-metabo.md` — package description
- `human/plans/cosmos-gem-integration.md` — GEM integration cross-repo plan

---

## Notes

- This plan prioritizes functionality over database integration
- Temporary solutions (pypath.inputs dependency, ID translation, location access)
  are contained within `datasets/cosmos/` to avoid spreading them across the package
- When proper resource clients, ID mapping, and database integration become available,
  only the `cosmos/` submodule needs refactoring
- Future branded datasets in `datasets/` can use either the temporary approach
  or the proper infrastructure depending on availability
- Web API development (Phase 5) will follow established OmniPath patterns

## Phase 2: Resource Client Expansion
### 2.1 STITCH
We have two versions of STITCH clients, one is the old one stored in pypath pypath-tool/pypath/inputs/stitch.py, and the other one is updated and stored in pypath-tool/pypath/inputs/new_stitch/. Compare their difference first, tell me which parts of functions are updated, if it is necessary to use new_stitch. Generate your investigation in a seperate markdown file named "STITCH_report.md". 
Then let's create a new STITCH COSMOS PKN build script by refering to scripts in Phase 1 using new_stitch. First I want to see outcome of new_stitch.interactions(). Save your investigation in STITCH_data_overview.md, and help me create a csv file with 5000 lines. 
At the moment, we don't add surfix, but in future we will add it when we know more about data. keep Gene__ and Metab__ prefixes same as before. Please normalize orientation: Always put chemical as source, protein as target. And don't create reverse reactions. Put score as a parameter and set 700 as default, use final_score column. Put mode as another parameter and set only keep activation and inhibition as default. Don't consider stereospecificity at present. Set species as a parameter and 9606 as default. For mor column, keep the same as original table: activation = 1, inhibition = -1, unknown = 0.



