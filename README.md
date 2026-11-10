# TANOS: TAxon jackknife for NOdal Stability

## Table of Contents
   [I. Introduction](#i-introduction)<br>
  [II. Installation Instructions](#ii-installation-instructions)<br>
 [III. Usage Instructions](#iii-usage-instructions)<br>
  [IV. License](#iv-license)<br>
  [V. Contributors](#v-contributors)<br>
 [VI. Citation](#vi-citation)<br>
[VII. Contact](#vii-contact)


## I. Introduction
In phylogenetics, assessing the stability of a node is paramount to having
confidence in the quality of the topology. Several strategies exist to evaluate
a node, with bootstrap support being particularly popular. Most methods rely on
how consistently characters from the input data matrix are in a given state to
assess the node. We seek to evaluate the stability of the node by determining how
resilient nodes are to the removal of taxa. This package can calculate the taxa
resiliency for a provided tree. To determine this, one must remove a taxon from
the alignment (jackknife sampling) and build a new tree, comparing it to the
original tree. To avoid anomalies, a new tree must be built multiple times (say
50). This process must be repeated for each taxon. Thus, if 100 taxa are in a
given tree and 50 replicates are performed for the removal of each taxon then
5,000 new trees are built. These new trees can be provided with the original tree
to our package for comparison to determine how robust a given clade is to varying
sampling. The intent is to assign a value to each node, similar to bootstrapping.
A value of 1 indicates 100% of sampled trees includes the grouping in that clade
(sans the one taxon intentionally removed). Valid values fit in the range [0,1].
Note that this value has no meaning for leaf nodes, parents of leaf nodes, and
the root node.

Please note that this package assumes the following has already been performed:

1. The original alignment of all taxa
2. A tree built with all taxa
3. "jackknifing" the alignment (making one copy per taxon in the tree, minus that taxon)
4. New trees built using the "jackknifed" alignments with a pre-determined level of replication

Assuming the interested user of this package wishes to generate nodal stability
scores for a tree, steps #1-2 will already be completed. Step #3 is trivially
accomplished with basic scripting. See the scripts directory of the repository
for an example. Step #4 is conceptually simple, but more difficult to manage as
the number of taxa and/or replicates increases. An example script to manage such
a task on a computing cluster under [SLURM](https://slurm.schedmd.com) control
using [IQ-TREE](http://www.iqtree.org) is also provided. Note that these are
simply example scripts, not fully-tested programs; modifications will be required
for your system and situation, especially if you intend to use an alternate
directory structure, tree inference software, or grid control engine. Please note
that if you installed TANOS with pip, you will not see these scripts as they are
not included in the PyPI manifest; you will need to download the repo from GitHub
to have access to them. Please see the
[README](https://github.com/pickettbd/tanos/blob/master/scripts/README.md) file
in the scripts directory for additional details.

This program must have access to all the relevant trees (_i.e._, the original/main
tree and the "jackknifed" trees) and a mapping of which trees are associated
with which taxa. A default directory structure and file naming scheme is assumed.
These assumptions can be modified to some degree, or bypassed entirely, with
certain command-line arguments. Please consult the section [III. Usage
Instructions](#iii-usage-instructions) for details.

We are in the process of publishing a paper describing this work. Once
available, we will provide relevant details here:

See our paper in [_journal_](#readme) for further information:<br>
<span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>[https://sub-domain.domain.tld/some/path/to/resource](#readme)


## II. Installation Instructions
This package is written in [Python](https://www.python.org). You must have a version of Python (v3.6+) that supports [f-strings](https://docs.python.org/3/reference/lexical_analysis.html#f-strings). This package also depends on the following Python modules: sys, re, pkgutil, argparse, and pathlib, which are all included in the Python Standard Library. Installation may be accomplished using pip like this:

<span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>`pip install tanos`

Alternately, you may install it manually using setuptools like this:

<span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>`python setup.py install`

If running `setup.py`, you may wish to add `--user` or `--prefix` for an
installation in non-default locations.

## III. Usage Instructions
The workflow's default input is `data/orig/supermatrix_dna.phy`, its output directory is `results`, and it runs 49 jackknife replicates per taxon. Set the input and output parameters when launching it:

```bash
nextflow run main.nf \
  -profile local \
  --input_alignment data/orig/supermatrix_dna.phy \
  --outdir results \
  --replicates 49
```

Use `-profile slurm` or `-profile aws` for the configured execution profiles. The pipeline runs these stages:

1. Select a model with IQ-TREE (`-m TESTONLY`).
2. Infer a main tree using the selected model.
3. Create one FASTA alignment per removed taxon.
4. Infer the requested number of IQ-TREE jackknife trees for each alignment.
5. Combine the jackknife tree files and run IQ-TREE's `-sup` support calculation on the main tree.

Outputs are published beneath the selected output directory in `modelTest/`, `mainTree/`, `jackknife/aln/`, taxon-specific directories, and `jackknifeSupport/`. Nextflow also retains its usual run metadata in its work directory and `.nextflow` files.

The alignment preparation script reads sequential PHYLIP input and checks that the declared taxon count, sequence lengths, and unique taxon names are consistent. The Nextflow jackknife-support stage uses IQ-TREE's support calculation; it does not invoke the Python `tanos` command or retain the mapping between each replicate tree and its removed taxon. To calculate TANOS resiliency from a generated tree set, arrange replicate trees by removed taxon (or create the required two-column file-of-filenames) and run the Python scorer separately.

**Output Formats:**
- **Newick**: the main tree with each node's `taxa-resiliency` metadata in comments by default, or scores substituted for branch lengths or internal labels.
- **JSON** and **pretty JSON**: a tree structure with labels, branch lengths, and node metadata, including the resiliency score.
- **Mermaid**: optional graph syntax; internal labels can be replaced by scores with `--replace-label`.


## IV. License
Please see the [LICENSE](https://github.com/pickettbd/tanos/blob/master/LICENSE).


## V. Contributors
Brandon Pickett (Github: [@pickettbd](https://github.com/pickettbd)) is the sole author of this software.


## VI. Citation
Until a version of the paper is available, we ask that you cite the Github
repository. Thereafter, you will be able to include a citation following
this general form:

> Powell GS, Pickett BD, Ridge PG, Martin GJ, Whiting MF, Bybee SM. Paper title. _journal_. Year. Volume(Issue):pages.


## VII. Contact
For questions, comments, concerns, etc., please open an issue on Github.

**Note:** Before asking for help, please ensure you (a) are using the latest
official release and (b) consult first the section [III. Usage Instructions](#iii-usage-instructions).

