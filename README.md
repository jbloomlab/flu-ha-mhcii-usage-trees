# Phylogenetic trees of influenza hemagglutinin annotated by MHC II binding
Build interactive Nextstrain phylogenetic trees of influenza HA from different subtypes and annotate by MHC II binding.
Trees are built and maintained by the [Bloom lab](https://jbloomlab.org/).


## Running pipeline
Structure of repo:
  - [config.yaml](config.yaml): configuration
  - [environment.yaml](environment.yaml): conda environment
  - [Snakefile](Snakefile): Snakemake file with pipeline
  - [./scripts/](scripts): scripts used by pipeline
  - [./data/](data): input data
  - `./results/`: results created by pipeline
  - [./auspice/](auspice): final JSON trees for visualization via [https://auspice.us/](https://auspice.us/) or [Nextstrain community builds](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html).

To run the pipeline, build and activate the conda environment in [environment.yaml](environment.yaml), and then run:

    snakemake --software-deployment-method conda -j <n_cpus>

## Steps in pipeline

### `get_ncbi_dataset_zip`: Download influenza sequences from GenBank
Download all Influenza A virus sequences and metadata from GenBank using `ncbi-datasets-cli`.
The download includes genomic sequences, CDS nucleotide sequences, gene annotations, and BioSample metadata as a single zip archive.
The taxon ID for the downloaded sequences is specified in [config.yaml](config.yaml) as *datasets_taxid*.
The datestamp file in [results/ncbi_dataset/ncbi_dataset_download_date.txt](results/ncbi_dataset/ncbi_dataset_download_date.txt) records the download date.

### `format_ncbi_dataset`: Process NCBI zip into TSV intermediates
Extract genome metadata and genomic FASTA headers from the NCBI zip file using `dataformat`.
This processes the large zip once so downstream per-subtype rules can work from smaller TSV files.

### `extract_ha_cds`: Extract HA sequences and metadata per subtype
For each tree defined in [config.yaml](config.yaml), extract full-length hemagglutinin (HA) coding sequences and associated metadata.
Subtypes are parsed from genomic FASTA headers and filtered by the `subtype` regex in the config (e.g., `H5N\d{1,2}` for all H5Nx) and the `cds_length_range` specifying the sequence length in the config; note also we drop any sequences with ambiguous nucleotides, missing dates, or that are not valid coding sequences.
Outputs per-tree FASTA and metadata TSV files in `results/trees/{tree}/`.

### `subsample`: Subsample sequences per subtype
Subsample the HA sequences for each tree using [augur subsample](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/subsample.html).
The subsampling strategy (grouping, max sequences, date filtering, include/exclude lists) is configured per tree in [config.yaml](config.yaml) under `augur_subsample`.
Outputs subsampled FASTA and metadata TSV files in `results/trees/{tree}/`.

### `align`: Align sequences
Align the subsampled HA CDS sequences for each tree using [augur align](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/align.html) (which wraps `mafft`).
A per-tree reference sequence (configured as `reference_sequence` in [config.yaml](config.yaml)) guides the alignment; the reference is removed from the output.
Outputs the aligned FASTA in `results/trees/{tree}/`.

### `tree`: Infer phylogenetic tree
Infer a maximum-likelihood phylogenetic tree from the alignment using [augur tree](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/tree.html) (which wraps IQ-TREE).
Uses a fixed seed for reproducibility and collapses zero-length branches.
Outputs a raw Newick tree in `results/trees/{tree}/` (before any temporal refinement).

### `refine`: Refine tree with temporal information
Refine the raw tree using [augur refine](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/refine.html) (which wraps TreeTime) to build a time-resolved phylogeny.
Outputs a refined Newick tree and a node-data JSON with branch lengths and inferred dates in `results/trees/{tree}/`.

### `ancestral`: Infer ancestral sequences
Infer ancestral nucleotide sequences using [augur ancestral](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/ancestral.html).
Mutations are reported relative to the same reference sequence used for alignment.
Outputs a node-data JSON with nucleotide mutations per branch in `results/trees/{tree}/`.

### `translate`: Translate to amino acid mutations
Translate nucleotide mutations to amino acid mutations using [augur translate](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/translate.html).
Uses a GFF3 annotation file to define gene regions (SigPep, HA1, HA2) for translation.
Outputs a node-data JSON with amino acid mutations per branch in `results/trees/{tree}/`.

### `score_mutation_effects`: Score nodes by mutation phenotype effects
Score each tree node based on the cumulative effect of its amino acid mutations on measured phenotypes.
Per-mutation effects are provided in a TSV file (configured as `mutation_effects` in [config.yaml](config.yaml)).
The script validates that proteins in the effects file match the genes in the tree's amino acid mutations.
For each phenotype (configured under `phenotypes`), computes the total effect, the max magnitude effect, and the identity of the max magnitude mutation.
Outputs a node-data JSON with phenotype scores in `results/trees/{tree}/`.

### `generate_phenotype_auspice_config`: Generate color scales for phenotype scores
Generate an auspice config JSON that defines continuous viridis color scales for each phenotype's total effect and max magnitude effect scores, and a categorical coloring for the max magnitude mutation.
The color scale is configured per tree as `phenotype_color_scale` in [config.yaml](config.yaml).
Outputs an auspice config JSON in `results/trees/{tree}/`.

### `export`: Export auspice JSONs
Export interactive auspice v2 JSONs using [augur export](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/export_v2.html).
Each tree uses an auspice config file (configured as `auspice_config` in [config.yaml](config.yaml)) that defines colorings, filters, display defaults, and metadata; a generated phenotype auspice config with color scales for mutation effect scores; and a `title` for the tree visualization.
The output files are placed in `auspice/` with names like `{auspice_prefix}_{tree}.json` (the prefix is set in [config.yaml](config.yaml), typically matching the repo name for [Nextstrain community builds](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html)).
These are the final pipeline outputs targeted by `rule all`.
