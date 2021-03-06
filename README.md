# PHYLDOG
PHYLDOG is a program allowing to simultaneously build gene and species trees when gene families have undergone duplications and losses. 
Trees and parameters are estimated in the maximum likelihood framework.

For instructions about installation, please see the INSTALL file.

### To use it, type:

```sh
mpirun -np NUM_PROCESSORS phyldog param=GeneralOptions.opt
```

Where the file GeneralOptions.opt is presented below.

Just typing:

```sh
./phyldog
```
should show a partial list of options that need to be given to the program.


## HOW TO RUN PHYLDOG

phyldog takes as input a series of files. One file contains general options that are used by all processors on which the job runs. Then there is one file per gene family, where gene family-specific options are given. 
This file containing the general options can be given as an argument to phyldog, as follows:


```
mpirun -np NUM_PROCESSORS phyldog param=GeneralOptions.opt
```

where NUM_PROCESSORS is the number of processors to be used by phyldog, and GeneralOptions.opt is the file containing the general options.

The file contains a list of options as follows:
first_parameter   = value1
second_parameter  = value2

This follows the syntax used for programs of the bppsuite series, as explained on the following webpage:
http://home.gna.org/bppsuite/doc/index.html#Top



## THE FILE CONTAINING GENERAL OPTIONS 


GeneralOptions.opt may contain the following list of options:

PATH= /home/user/path_to_data_files/ #path to the directory where the input files are, and the output files are left.

init.species.tree=user #whether the input species tree is given in a file, or should be "random"

species.tree.file=$(PATH)InputSpeciesTree.tree #gives the path to the input species tree, in newick format.

species.names.file=$(PATH)SpeciesNames.txt # gives the species names to be considered for species tree reconstruction. One species name per line, without spaces.

starting.tree.file=$(PATH)start #the starting species tree is saved in the file start, in the directory $(PATH)

output.tree.file=$(PATH)output.sptree #the end species tree is saved in the file output.sptree, in the directory $(PATH)

output.temporary.tree.file=$(PATH)CurrentSpeciesTree.tree #the species tree obtained is output to this file if the job has not finished but was closed because of time constraints.

genelist.file=$(PATH)listeGene #Contains a list of gene family-specific options

output.duplications.tree.file=$(PATH)SpeciesTreeDuplications.tree # Species tree where branch lengths represent total numbers of duplications

output.losses.tree.file=$(PATH)SpeciesTreeLosses.tree # Species tree where branch lengths represent total numbers of losses

output.numbered.tree.file=$(PATH)SpeciesTreeLosses.tree # Species tree with nodes numbered

optimization.topology=no # whether the species tree topology should be optimized (yes) or not (no)

branch.expected.numbers.optimization = average_then_branchwise # whether the branch-wise parameters of duplications or losses should be optimized and branchwise (branchwise) or optimized and averaged over all branches (average) or not optimized (no). The option "average_then_branchwise" is a good compromise between speed and accuracy.

genome.coverage.file= $(PATH)GenomeCoverage # File giving the expected completeness of the genomes under study, in percents. 

spr.limit=5 # For SPR moves on the species tree, gives the maximum distance between the position of the pruned subtree and its regrafting position.

time.limit=23 # Time limit for the job: beyond 23 hours, the job stops. Should be useful if the job is limited to less than 24 hours

current.step=0 # This option is useful to restart a job that has been stopped due to time.limit. 

output.file.suffix=_extension # An extension that will be added to all output files.

alternate.topology.likelihoods=$(PATH)alternateLks.txt # A file where the likelihoods of alternate topologies encountered during the NNI search on the species tree are saved.

species.duplication.tree.file=previousDuplicationTree.nwk # A file containing a species tree with branch lengths representing duplication parameters. These duplication parameters are then used as starting values for the algorithm. You don't have to give this option if you don't have good estimates for these parameters. 

species.loss.tree.file=previousDuplicationTree.nwk # A file containing a species tree with branch lengths representing loss parameters. These loss parameters are then used as starting values for the algorithm. You don't have to give this option if you don't have good estimates for these parameters.

This GeneralOptions.opt file contains options specific to the search for the best species tree. However, the options included in this file are also read by client processors in charge of gene families. Therefore, it is possible to include options that apply to the gene tree search for all gene families.




## THE FILE LISTING GENE FAMILY-SPECIFIC OPTION FILES


Gene-family specific options need to be given in additional files. The list of these files is given in "genelist.file".
"genelist.file" should look like this:

family_1.option

family_2.option
...

or 


family_1.option:10

family_2.option:30

...

The second way to list the option files contains an additional element of information, which is the "complexity" of a gene family. This "complexity" is used at the beginning of the algorithm to distribute equally the loads on the different computers. So far, it is still unclear what this complexity should be (some function of the number of sequences and the number of sites). I generally use the number of sequences in the gene family for lack of a better metric.



## GENE FAMILY-SPECIFIC OPTION FILES



Inside these option files, gene-family-specific options should be set. These options again follow the bppsuite syntax. A large number of these options are documented in the bppsuite help:
http://home.gna.org/bppsuite/doc/index.html#Top

An example of these options is given below (let's assume it is the contents of the file "family_1.option"):

######## First, data files ########

PATH= home/user/path_to_data_files/ #path to the directory where the input files are, and the output files are left.

DATA= family_1 # Variable used to give the name of the data files.

alphabet=DNA # Could also be "RNA", "protein", or Codon. Please see the bppsuite help for more details.

taxaseq.file=$(PATH)$(DATA).link # File giving the link between species and sequence names (more on this below)

input.sequence.file=$(PATH)$(DATA).fasta # file giving the input sequence alignment for the gene family

input.sequence.format=Fasta # Format of the sequence alignment. Could be Fasta, Phylip, Clustal, Mase, Nexus... Please see the bppsuite help for more details.

output.reconciled.tree.file=$(PATH)$(DATA)_Reconciled.tree # File where to store the output improved and reconciled gene tree, in NHX format. Duplication and speciation nodes are annotated, with the tag "Ev=D" or "Ev=S" respectively.

output.duplications.tree.file=$(PATH)$(DATA)_Duplications.tree # File where the species tree topology is saved, annotated with numbers of duplications for this gene family.

output.losses.tree.file=$(PATH)$(DATA)_Losses.tree # File where the species tree topology is saved, annotated with numbers of losses for this gene family.

output.numbered.tree.file=$(PATH)$(DATA)_Numbered.tree # File where the species tree topology is saved, annotated with node indices.

output.events.file=$(PATH)$(DATA)_Events.txt # File where events of duplication and loss are written, along with the species ID information. The format is one event per line, with an event described as: event(SpeciesID, "FamilyName", duplication|loss).

output.orthologs.file=$(PATH)$(DATA)_Orthologs.txt # File where orthologs and paralogs are written. The format is one orthology/paralogy relationship per line, with first the family name, the type of relationship (Orthology or paralogy) then a list of genes, "<===>" and the other series of genes that are in relationship to the first ones.

input.sequence.sites_to_use=all # tells whether we should use all sites in the alignment or not. Could be "all", "nogap", or "complete". Please see the bppsuite help for more details.

input.sequence.max_gap_allowed=100% # Maximum number of gaps tolerated for including a site in the analysis.

init.gene.tree=user # Starting gene tree. Could be "user", "bionj" or "phyml". "user" requires that a user-input tree is given with the "gene.tree.file" option, whereas the options "bionj" and "phyml" have phyldog use these algorithms to create starting gene trees.

gene.tree.file=$(PATH)$(DATA).tree # File containing the input starting gene tree in newick format. Useful if "init.gene.tree=user".

output.starting.gene.tree.file=$(PATH)$(DATA)_starting.tree # File where the starting gene tree is saved.

######## Then, algorithm options ########

rearrangement.gene.tree = nni # Type of rearrangement: "nni" or "spr". "nni" is much faster but less exhaustive than "spr". If the species tree topology is fixed, we advise spr, which provides better gene trees. Otherwise, we advise the use of nnis.

SPR.limit.gene.tree = 4 # For SPR moves on the gene tree, gives the maximum distance between the position of the pruned subtree and its regrafting position.

######## Then, model options ########

model=GTR(a=1.17322, b=0.27717, c=0.279888, d=0.41831, e=0.344783, initFreqs=observed, initFreqs.observedPseudoCount=1) # options of the model used. Should match the alphabet. Please see the bppsuite help for more details.

rate_distribution=Invariant(dist=Gamma(n=4,alpha=1.0), p=0.1) # Rate heterogeneity option. Here we assume a gamma law with 4 categories and a category of invariants to model rate heterogeneity among sites.

optimization.ignore_parameter=InvariantMixed.dist_Gamma.alpha, InvariantMixed.p, GTR.a, GTR.b, GTR.c, GTR.d, GTR.e, GTR.theta, GTR.theta1, GTR.theta2 # We choose not to optimize these 10 parameters in order to save computing time, as we have provided reasonable input values. However, in cases where good input values are not available, it may be wise to leave this field empty and optimize these parameters.


######## Finally, optimization options ########

optimization.topology=yes # We choose to optimize the topology.

optimization.tolerance=0.01 # We have a large optimization tolerance to speed up the computations.

#### The options below are also tuned to speed up the computations, and may be left untouched. More details on what they mean may be found in the bppsuite help.

optimization.method_DB.nstep=0

optimization.topology.numfirst=false

optimization.topology.tolerance.before=100

optimization.topology.tolerance.during=100

optimization.max_number_f_eval=1000000

optimization.final=none

optimization.verbose=0

optimization.message_handler=none

optimization.profiler=none

optimization.reparametrization=no

## THE FILE GIVING THE LINK BETWEEN SPECIES NAME AND SEQUENCE NAME



The "taxaseq.file" file contains the link between species name and gene name.
For instance: 

Oryctolagus_cuniculus:ENSOCUP00000017695

Dipodomys_ordii:ENSDORP00000011323

Sus_scrofa:ENSSSCP00000001753

Pongo_pygmaeus:ENSPPYP00000018557;ENSPPYP00000018560

...

This means that sequences ENSPPYP00000018557 and ENSPPYP00000018560 correspond to the species Pongo_pygmaeus for instance.

The species names should correspond to the names used in "species.names.file", an option that should found in GeneralOptions.opt, or alternatively should correspond to the names in the input species tree ("inout.species.tree").




## SUMMARY: FILES NEEDED BY PHYLDOG


### Sum-up on the files really needed to run PHYLDOG:
- a file GeneralOptions.opt
- a file giving the input species tree, or alternatively a list of species names (if the option for the input species tree is "random", this file is necessary). Otherwise this file can be used to limit the list of species to include in the study.
- a file giving a list of gene family-specific option files 
- one file per gene family describing the options for this gene family
- one alignment file per gene family
- one file giving the link between species name and sequence name, one per gene family

### Files that can be provided but that are not absolutely necessary:
- a file giving genome coverage per species (otherwise phyldog assumes genomes are 100% complete)
- a file giving the duplication parameters for the species tree (otherwise phyldog will estimate these)
- a file giving the loss parameters for the species tree (otherwise phyldog will estimate these)
- one tree file per gene family (otherwise phyldog can estimate gene tree using bionj or phyml-like algorithms)
