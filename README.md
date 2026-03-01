# Phylogenetic trees of influenza hemagglutinin labeled by MHC II binding
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
