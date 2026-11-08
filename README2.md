# TANOS: TAxon jackknife for NOdal Stability

TANOS evaluates how stable internal clades are when taxa are removed from a phylogenetic alignment. This repository contains two related but distinct workflows: a Python program that calculates TANOS taxon resiliency scores, and a Nextflow pipeline that calculates IQ-TREE support values from inferred jackknife trees. They use jackknife trees in different calculations, so their outputs are not interchangeable.

## How the score works

For each eligible internal clade in the main tree, the Python scorer considers each taxon in that clade in turn. It checks the trees inferred after removing that taxon for the grouping formed by the clade's other taxa. The resiliency score is the fraction of those relevant replicate trees that contain the grouping, on a scale from 0 to 1. The root is assigned 0; leaves and nodes without grandchildren do not have a meaningful score (the latter are assigned 1 by the implementation).

The Python scorer expects the main tree and all jackknifed trees to have already been generated. The Nextflow workflow automates alignment jackknifing and IQ-TREE inference, then asks IQ-TREE to calculate support on the main tree.

### How the two results differ

| Python `tanos` scorer | Nextflow pipeline support stage |
| --- | --- |
| Calculates the TANOS `taxa-resiliency` score for each eligible clade. | Runs IQ-TREE 2 with jackknife trees supplied through `-z` and the main tree as the support target (`-sup`). |
| Associates replicates with the identity of the taxon removed, then checks the corresponding reduced clade. | Concatenates the jackknife tree files into one replicate-tree collection for IQ-TREE. The process does not pass removed-taxon identity into the support calculation. |
| Reports the fraction of relevant replicates supporting the reduced clade, from 0 to 1. | Reports support values produced by IQ-TREE; their exact interpretation and output format are those of the installed IQ-TREE version and command. |

Therefore, the Nextflow pipeline does calculate support values in its final stage, but it does not calculate the Python scorer's TANOS resiliency metric. Use the Python command when the desired result is the taxon-specific resiliency score. Use the Nextflow result when the desired result is IQ-TREE support from the supplied jackknife tree set.

## Requirements

### Python scorer

- Python 3.6 or newer
- No third-party Python runtime packages are required; the scorer uses Python standard-library modules.

### Nextflow workflow

- Nextflow
- IQ-TREE 2 available as `iqtree2` on the execution environment's `PATH`
- Python 3 for the alignment jackknifing process
- An input alignment in sequential PHYLIP format for the included alignment preparation script

The Nextflow profiles are `local`, `slurm`, and `aws`. The Slurm profile selects queue `m9`; adjust that for your cluster. The AWS profile selects Nextflow's `awsbatch` executor and requires the corresponding AWS Batch configuration in your environment.

## Install the Python command

From the repository root, install the package with pip:

```bash
python -m pip install .
```

Or install from a local checkout with setuptools:

```bash
python setup.py install
```

This installs the `tanos` command. The package version is recorded in `VERSION` (currently `0.0.1-beta7`).

## Run the Python scorer

```bash
tanos --help
```

By default it reads `data/mainTree/tree.nwk` and searches `data/jackknife/tree/` for directories named after each taxon. Each taxon directory should contain Newick trees named `tree-1.nwk`, `tree-2.nwk`, and so on. The tree leaf labels and directory names are case-sensitive and must match. Every taxon must have the same number of replicate trees.

For example, with a main tree at `results/main.nwk` and jackknife trees under `results/jackknife/`:

```text
results/
├── main.nwk
└── jackknife/
    ├── taxon_A/
    │   ├── tree-1.nwk
    │   └── tree-2.nwk
    └── taxon_B/
        ├── tree-1.nwk
        └── tree-2.nwk
```

Run with explicit paths and output filenames:

```bash
tanos \
  --main-tree results/main.nwk \
  --jackknife-tree results/jackknife \
  --output-nwk results/scored.nwk \
  --output-json results/scored.json \
  --output-json-pretty results/scored-pretty.json
```

The scorer also accepts a tab-separated tree file-of-filenames instead of the directory layout. Supply exactly two columns per line, with no header: the taxon label followed by the path to one Newick replicate tree.

```text
taxon_A\tresults/jackknife/taxon_A/tree-1.nwk
taxon_A\tresults/jackknife/taxon_A/tree-2.nwk
taxon_B\tresults/jackknife/taxon_B/tree-1.nwk
```

```bash
tanos --main-tree results/main.nwk \
  --jackknife-tree-fofn results/jackknife-trees.tsv
```

By default, the command writes `out.nwk`, `out.json`, and `out_pretty.json` in the current directory. Use an empty string for an output option to skip that output, for example `--output-mmd ""` (Mermaid output is disabled by default). Newick scores are added as comments by default. `--replace-branch-len` writes scores in place of branch lengths; `--replace-label` writes them in place of internal node labels. Other options include `--tree-ext` to change the extension expected by directory scanning, and `--version`, `--license`, and `--cite` for package information.

## Run the Nextflow workflow

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

## Output formats

- **Newick**: the main tree with each node's `taxa-resiliency` metadata in comments by default, or scores substituted for branch lengths or internal labels.
- **JSON** and **pretty JSON**: a tree structure with labels, branch lengths, and node metadata, including the resiliency score.
- **Mermaid**: optional graph syntax; internal labels can be replaced by scores with `--replace-label`.

## Project layout

```text
src/tanos/       Python CLI and Newick tree implementation
processes/       Nextflow processes and alignment jackknifing script
main.nf          Nextflow workflow definition
nextflow.config  Default parameters and local, Slurm, and AWS Batch profiles
test/            Informal Python tree parsing and traversal test scripts/data
```

## Citation and contact

Until the associated paper details are available, cite the [TANOS GitHub repository](https://github.com/pickettbd/tanos). The citation template in the project is:

> Powell GS, Pickett BD, Ridge PG, Martin GJ, Whiting MF, Bybee SM. Paper title. *journal*. Year. Volume(Issue):pages.

For questions or issues, open an issue in the [GitHub repository](https://github.com/pickettbd/tanos/issues). TANOS is distributed under the MIT License; see [LICENSE](LICENSE).
