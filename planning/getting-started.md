# Python package for metabolite and compound structure, reaction, and interaction related prior-knowledge database access, processing and web service

We are developing a software ecosystem to support systems biology and
biomedicine data processing and analysis. It consits of several components,
which work either in client or server side. We are in the process of building a
large combined database that includes metabolite, compound, reaction, and
interaction related prior-knowledge. The database, called OmniPath, uses
Postgres and multiple web APIs operate on it. Some of the APIs might accept
data upload or do processing of queried data, such as combining multiple query
results or user data or using RDKit for cheminformatics.

After a few exploratory iterations, now we are ready to create a fresh package
which we intend to be the design that goes into deployment soon. The uneven
status (ready vs. missing) of a large variety of related components dictates
the order of implementation of various parts of this package. In the first
phase, we won't connect to the database yet, and won't have a web API yet.
Instead, we'll start by incorporating certain processing algorithms, and before
the database and web API will be integrated, we'll connect these temporarily
directly to our resource (database) access modules and implement a minimal web
API. The newly created package will be called `omnipath-metabo`, and we intend
it to operate primarily on server side, but reserving the possiblity  to use it
also client side. For this reason, to make it attractive to industry as well,
we will use a BSD 3-Clause license.

Outside of this package, we've developed a set of Python scripts to build most
parts of the COSMOS PKN. These we should incorporate, organize, and refactor in
a new submodule within the metabo package. This main submodule could be one
that compiles more specialized prior knowledge datasets that go beyond simple
SQL queries and require postprocessing. We'll have more of these kind of custom
prior-knowledge, one of these is the COSMOS PKN. The COSMOS PKN build scripts
use the resources access modules available in the `pypath` package. You can
browse the code of this package under `../pypath/pypath`, specifically the
resources clients are located under the `../pypath/pypath/inputs` directory.

The COSMOS PKN is a specialized prior-knowledge network that we compile from
several sources. It is used for network optimization (mechanistic modeling,
causal reasoning) in a multi-omics integration context, where the network
covers signaling, gene regulation, kinase-substrate interacitons, metabolism,
and any metabolite-protein or protein-metabolite interactions, such as
transporters and allosteric regulation.  One preliminary script which includes
TCDB and SLC-table in the COSMOS PKN is available here at
`scripts_exploration/transporter.py`, while another draft script covers
allosteric regulation interactions from BRENDA is available at
`scripts_exploration/allosteric_regulation.py`, and receptor-metabolite
interactions from MRCLinksDB is available at
`scripts_exploration/receptors_metabolites.py`, and human lipid metabolism in a
graph format is available in `scripts_exploration/lipid_rhea.ipynb`.  Some
further databases that will be part of this network are Reactome, STITCH,
SwissLipids (for lipid hierarchy), and the genome scale metabolism models
(GEMs). For some of these already clients are implemented in the mentioned
`pypath.inputs`module, some are completely missing at the moment.

We start all of our packages from a cruft/cookiecutter project template which
contains modern tooling and CI (uv, ruff, pytest, GitHub workflows, mkdocs
material). You can find this template in
https://github.com/saezlab/python-project or ../pypath-new/_project_template/.

In general, please consider the `planning/coding-style.md` for every Python
code in this package.

The first steps of establishing this package:

1. Instantiate the template
2. Create a draft module/directory structure
3. Move here COSMOS PKN build scripts
4. Make refactoring plan: organization of functions, classes, modules, public
   API design, code quality, robustness, maintainability, etc.
5. Create testing plan (for unit tests)
