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
Outputs per-tree FASTA and metadata TSV files in `results/{tree}/`.
