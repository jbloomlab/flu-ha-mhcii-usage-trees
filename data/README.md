# Input data

This subdirectory has input data to the pipeline.

## Auspice configuration

[./auspice_config/](auspice_config) contains files relevant to configuration of the exported Auspice JSONs.

## Mutation effect data
The input (experimental) data used for the mutation effects is in [./mutation_effect_data/](mutation_effect_data):

  - [mutation_effect_data/H7_data.csv](mutation_effect_data/H7_data.csv): file downloaded from [https://github.com/dms-vep/Flu-H7-Anhui13-MHCII-binding/blob/master/results/summaries/tufted_duck_MHCII_binding.csv](https://github.com/dms-vep/Flu-H7-Anhui13-MHCII-binding/blob/master/results/summaries/tufted_duck_MHCII_binding.csv)

  - [mutation_effect_data/H5_data.csv](mutation_effect_data/H5_data.csv): file downloaded from [https://github.com/dms-vep/Flu-H5N1-American-Wigeon-2021-HA-MHCII-DMS/blob/main/results/summaries/tufted_duck_MHCII_binding.csv](https://github.com/dms-vep/Flu-H5N1-American-Wigeon-2021-HA-MHCII-DMS/blob/main/results/summaries/tufted_duck_MHCII_binding.csv)


## Strain-token host map

[./strain_token_host_map.tsv](strain_token_host_map.tsv) is a shared hand-curated table that the `infer_missing_host` rule uses to fill `host` / `host_tax_id` for rows where NCBI's structured fields are empty. For each such row the pipeline parses the 2nd `/`-delimited token from the strain name (e.g. `chicken` in `A/chicken/Scotland/1959`), lowercases and trims it, and looks it up here.

Columns:

  - *token*: the lowercased, whitespace-stripped 2nd `/`-delimited token from the strain name.
  - *host_tax_id*: NCBI Taxonomy ID to assign. **Leave blank to explicitly skip** the token (e.g. locations like `california`, environmental samples like `environment`, processed products like `pet food`, or ambiguous terms like `turkey`). Blank rows serve as explicit documentation that the token was considered.
  - *host_scientific_name*: NCBI Taxonomy scientific name for `host_tax_id` (used to populate the `host` column so it matches NCBI's convention). Must be blank iff `host_tax_id` is blank.
  - *notes*: free-text notes (e.g. `location`, `environmental sample`, `typo for chicken`).

Tax IDs should be chosen at a granularity sufficient to produce the correct `host_general` classification in `annotate_host_taxonomy` (species-level for confident identifications, genus/family/order for generic terms, class `Aves` 8782 for very generic avian terms). Tokens not present in this file are treated as unknown: the row is left empty-host and the token is logged as a warning so the file can be extended.

## Tree input data

[./trees/](trees) has the input data for each tree specified under *trees* in [../config.yaml](../config.yaml).
There is a subdirectory for each tree in *trees* (eg, [./trees/H5/](trees/H5) for the *H5* tree).
Each subdirectory has the following files:

  - `reference_sequence.fa`: FASTA with coding sequence of HA chosen as the alignment reference for this subtype. You need to manually decide what HA to use for the reference sequence for your subtype, download it, and include it here. Manually edit the FASTA header to have "reference_sequence" as the first word in the header.

  - `annotation.gff`: a GFF file annotating the proteins in `reference_sequence.fa`. Should annotate *SigPep* (signal peptide), *HA1* (HA1), and *HA2* (HA2). See [this paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4229193/) if you are trying to determine where the signal peptide ends and HA1 (ectodomain) begins.

  - `manual_add_sequences.fa` and `manual_add_metadata.tsv`: FASTA and TSV file giving CDS sequences to manually add to the tree; these are always included. Use this if there are sequences lacking Genbank accessions that you want to include. The `manual_add_sequences.fa` should have the accession as the header, and the TSV file should give the metadata with columns including. If these sequences do not have accessions (eg, they are not on Genbank) you can make accession any string identifier with no spaces:
    + *accession*
    + *strain*
    + *subtype*
    + *date*
    + *host*
    + *host_tax_id*
    + *location*
    + *region*
    + *passage_history*
    + *length*

  - `accessions_to_exclude.txt`: a list of Genbank accessions to exclude from the final tree, listed one-per-line. You should manually add to this file if there are outlier or incorrect sequences you want to exclude from the tree. If there are no sequences to exclude, make this an empty file.

  - `strain_annotations.tsv`: a TSV file with one row per annotated strain, representing directly measured properties to map onto the tree tips. Must include an `accession` column (in any position) whose values match tip accessions in the final tree; every other column becomes a per-tip annotation on the exported Auspice JSON. Numeric columns are colored using the shared `phenotype_color_scale` from `../config.yaml`; non-numeric columns become categorical colorings and are added to the Auspice filters. Make this an empty file (header-only, or fully empty) if there are no strain annotations for a subtype. How the pipeline handles accessions in this file that are not tips in the final tree is controlled by `missing_strain_annotations_action` in `../config.yaml`.
  
  - `accessions_to_include.txt`: a list of Genbank accessions to always include in the final tree (exempt from subsampling filters and the molecular clock filter), listed one-per-line. Lines beginning with `#` (or any trailing `# …` portion of a line) are treated as comments and ignored, so you can annotate accessions with the strain name. You should manually add to this file if there are sequences you want to guarantee are in the tree. If there are none, make this an empty file. The pipeline merges this list with all accessions in `manual_add_metadata.tsv` into `results/trees/{tree}/accessions_to_include_combined.txt`, which is what is actually passed to `augur subsample --include` and `augur refine --keep-ids`. The `validate_includes` rule checks that every accession in this list and in `manual_add_metadata.tsv` appears in the final Auspice JSON; the pipeline errors if any are missing.

  - `cds_length_range.yaml`: YAML file specifying the allowed CDS nucleotide length range for filtering sequences. Has a single key `cds_length_range` with a two-element list `[min, max]`. Either bound can be `null` to impose no limit on that side (e.g., `[1692, null]` for no upper bound). To set this length range for a new subtype, first set contents to `cds_length_range: [null, null]` and then look at the output in `results/trees/{tree}/extract_cds_stats.txt` which prints distribution of length ranges to set reasonable limits.

  - `protein_sites.tsv`: TSV file with the following columns:
    + *reference_sequential_site*: 1, 2, numbering of the protein encoded by `reference_sequence.fa`
    + *reference_aa*: amino acid in the protein encoded by `reference_sequence.fa` at this site
    + *protein*: the region of HA, should be *SigPep*, *HA1*, or *HA2* like in `annotation.gff`.
    + *protein_site*: 1, 2, ... numbering of the *protein*
    + *ectodomain_site*: 1, 2, ... numbering starting with the HA ectodomain (start of HA1), with sites prior to the ectodomain (signal peptide) numbered as ..., -2, -1.
    + *H3_site*: site number in H3 numbering, which is equivalent to 1, 2, ... numbering of the ectodomain of the H3 HA.

  - `mutation_effects.tsv`: TSV file that gives the mutation effects. It has the following columns:
    + *protein*: the protein (*SigPep*, *HA1*, or *HA2*) as in `protein_sites.tsv`
    + *protein_site*: the protein site as in `protein_sites.tsv`
    + *reference_aa*: the reference amino acid at the site as in `protein_sites.tsv`
    + *H3_site*: the H3 site number as in `protein_sites.tsv`
    + *wildtype*: the wildtype amino-acid at the site in the experiment used to measure the mutation effect
    + *mutant*: the mutant amino-acid at the site in the experiment used to measure the mutation effect
    + *phenotype*: the phenotype being measured; these phenotypes must include all phenotypes provided under the *phenotypes* key for the *trees* in [../config.yaml](../config.yaml).
    + *effect*: number giving the effect of the mutation on the phenotype

To add a new tree, manually add `reference_sequence.fa`, `annotation.gff`, `accessions_to_exclude.txt`, and `cds_length_range.yaml`. Then run the script described in next subsection ([build_protein_sites_and_mutation_effects.py](build_protein_sites_and_mutation_effects.py)) to build the `protein_sites.tsv` and `mutation_effects.tsv` files.

## Script to build `protein_sites.tsv` and `mutation_effects.tsv`
The `protein_sites.tsv` and `mutation_effects.tsv` files for the tree input data are built by the script [build_protein_sites_and_mutation_effects.py](build_protein_sites_and_mutation_effects.py) from the input data in [./mutation_effect_data/](mutation_effect_data) and the configuration in [build_protein_sites_and_mutation_effects_config.yaml](build_protein_sites_and_mutation_effects_config.yaml) by running:

    python build_protein_sites_and_mutation_effects.py

This script does the following:

 1. Reads each reference sequence in [./trees/](trees) and translates the CDS into a protein with BioPython using `cds=True` so an error is raised if it does not translate as a valid CDS. The proteins are named by the subdirectory names in [./trees/](trees), eg `H5`, etc.

 2. Reads the YAML configuration in [build_protein_sites_and_mutation_effects_config.yaml](build_protein_sites_and_mutation_effects_config.yaml). This file has a key *mutation_effect_data* with entries that each specify a *csv* file and a set of *phenotypes* (mapping CSV column names to display names). Each CSV must exist and have columns matching the phenotype keys, plus `sequential_site`, `wildtype`, and `mutant` columns. The final display names mapped by *phenotypes* must be unique across all entries in *mutation_effect_data*.

 3. For each CSV in *mutation_effect_data*, reconstructs a protein sequence using the *sequential_site* and *wildtype* columns; if there are any sites where *sequential_site* is missing in sequential numbering, an ambiguous amino acid `X` is inserted. These proteins are named according to the keys in *mutation_effect_data* (eg, `H5_data`); an error is raised if these names overlap with the reference sequence names.

 4. Uses `mafft` to make a multiple sequence alignment of all proteins (both reference sequences and reconstructed experiment proteins), and saves this to `build_protein_sites_and_mutation_effects_alignment.fa`.

 5. Reads the `annotation.gff` file for each tree to map sequential 1, 2, ... numbering of the reference protein to *protein* (*SigPep*, *HA1*, *HA2*) and *protein_site* (1, 2, ... numbering within each protein).

 6. Uses the alignment to compute H3 numbering for each tree. H3 site numbering is the ectodomain numbering of the H3 HA (HA1 starts at site 1, signal peptide sites are negative). If a tree's reference has an insertion relative to H3 at a given alignment position, the H3 site is numbered with letter suffixes like `125a`, `125b`, etc.

 7. Builds `protein_sites.tsv` for every tree, writing *reference_sequential_site*, *reference_aa*, *protein*, *protein_site*, *ectodomain_site*, and *H3_site*.

 8. Builds `mutation_effects.tsv` for every tree. Each mutation effect dataset is mapped to every tree through the alignment (not just the same subtype). For each mutation in each dataset:
    - The experiment's *sequential_site* is mapped to the tree's *reference_sequential_site* via the alignment (columns where both the experiment protein and tree reference are non-gap).
    - Experiment sites that are insertions relative to a tree's reference (gap in reference at that alignment column) are skipped for that tree.
    - Rows where *wildtype* equals *mutant* are excluded.
    - Rows where the *effect* value is missing (NaN) for a given phenotype are excluded for that phenotype.
    - The output columns are *protein*, *protein_site*, *reference_aa*, *H3_site*, *wildtype*, *mutant*, *phenotype* (the display name from the config), and *effect*.

 9. Writes `build_protein_sites_and_mutation_effects_protein_identities.txt` with pairwise percent identity at aligned (both non-gap) sites for all pairs of proteins in the alignment, grouped by each protein for easy reading.
