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
A datestamp file in [results/ncbi_dataset](results/ncbi_dataset) records the download date.
