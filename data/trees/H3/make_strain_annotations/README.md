# Build strain annotations

[virus_titer_QC_pass.csv](virus_titer_QC_pass.csv) is the CSV file taken from [here](https://github.com/dms-vep/Flu-H5N1-American-Wigeon-2021-HA-MHCII-DMS/blob/main/analysis_notebooks/H1_H2_H3_MHCII_entry_titers/results/virus_titer_QC_pass.csv) that has Bernadeta's measurement of how well different strains can infect different 293 cell variants.

The marimo notebook [make_log10_titers_relative_to.py](make_log10_titers_relative_to.py) takes as input [virus_titer_QC_pass.csv](virus_titer_QC_pass.csv) and creates [make_log10_titers_relative_to.csv](make_log10_titers_relative_to.csv), which has the log10 titers of the strains relative to normalization strains.

The script [map_strains_to_accessions.py](map_strains_to_accessions.py) then takes as input [make_log10_titers_relative_to.csv](make_log10_titers_relative_to.csv) and maps strain names to accessions used by the pipeline, and writes the titer data for all strains to can be matched to accessions to [../strain_annotations.tsv](../strain_annotations.tsv). The file [map_strains_to_accessions_summary.txt](map_strains_to_accessions_summary.txt) summarizes any strains that cannot have their accessions determined.
