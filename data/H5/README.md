# Input data for H5 tree

- [accessions_to_exclude.txt](accessions_to_exclude.txt): accessions manually specified for exclusion from trees as they appear to be erroneous sequences.

- [reference_sequence.fasta](reference_sequence.fasta): the reference sequence, HA cds from A/American Wigeon/South Carolina/22-000345-001/2021(H5N1), [Genbank OQ958044.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ958044.1). For the tree, sequences are aligned to this.

- [protein_sites.tsv](protein_sites.tsv): file giving the mapping of protein sites between 1, 2, ... numbering of the protein sequence encoded by [reference_sequence.fasta](reference_sequence.fasta) and other numbering schemes. Derived from numbering schemes described [here](https://dms-vep.org/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/numbering.html). Columns are:
  - *sequential_site*: 1, 2, ... numbering of protein encoded by [reference_sequence.fasta](reference_sequence.fasta)
  - *H3_site*: equivalent number in H3 numbering
  - *H5_site*: equivalent number in mature H5 numbering (ectodomain)
  - *protein*: domain of HA: *SigPep*, *HA1*, *HA2*
  - *protein_site*: sequential 1, 2, ... numbering of *protein*

- [annotation.gff](annotation.gff): GFF3 annotation defining the gene regions (SigPep, HA1, HA2) with nucleotide coordinates in the reference sequence. Generated from [protein_sites.tsv](protein_sites.tsv) by running:
  ```
  python protein_sites_to_gff.py protein_sites.tsv reference_sequence.fa annotation.gff
  ```
- [mutation_effects.tsv](mutation_effects.tsv): File giving effects of individual mutations on properties of interest. It has columns *sequential_site*, *protein*, and *protein_site* that indicate the site number, and then additional columns giving the wildtype and mutant residue at each site and their effect on various phenotypes. This file is generated from [tufted_duck_MHCII_binding.csv](tufted_duck_MHCII_binding.csv) (which was downloaded from [here](https://github.com/dms-vep/Flu-H5N1-American-Wigeon-2021-HA-MHCII-DMS/blob/main/results/summaries/tufted_duck_MHCII_binding.csv)) by running:
  ```
  python get_mutation_effects.py protein_sites.tsv tufted_duck_MHCII_binding.csv mutation_effects.tsv
  ```

