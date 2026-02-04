# Initial plan for COSMOS PKN build implementation in omnipath-metabo

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

Processing and integrations of several resources have been implemented in the
`OmnipathR` package; on the long term we intend to develop only the Python
solution in the current project, hence the old R code you can consider as a
incomplete legacy code, which still carries important information about
resources processing, integration and the COSMOS PKN itself. You can find the
relevant code in `../omnipathr/`, specifically in the `cosmos.R`, `stitch.R`,
`recon3d.R`, `gem.R`, `chalmers-gem.R`. An even earlier, original
implementation of the COSMOS PKN build you can find under the
`../omnipathr/cosmos` directory. These are R scripts with lot of repetitive
code and suboptimal patterns, however they serve as an ultimate reference how
the network should look like. This includes the addition of suffixes or
attributes to enzymes and metabolites, aiming to 1) condider reactions and
interactions in a compartment specific way and 2) to separate the paths in the
causal network that impact the same enzymes or metabolites but in different
reactions, and to avoid creating loops from bi-directional reactions. A few
sample edges from the ready COSMOS PKN in its original format is available in
the `planning/cosmos-node-attrs.md`. However, in the new implmentation we aim
for a more robust and clean design, which means I would add the various
prefixes and suffixes to node IDs only at the very end and only if it's
absolutely necessary, and instead, I would keep each attribute in separate
columns, as it's done in the `OmnipathR` implementation
