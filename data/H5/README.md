# Input data for H5 tree

- [accessions_to_exclude.txt](accessions_to_exclude.txt): accessions manually specified for exclusion from trees as they appear to be erroneous sequences.

- [reference_sequence.fasta](reference_sequence.fasta): the reference sequence, HA cds from A/American Wigeon/South Carolina/22-000345-001/2021(H5N1), [Genbank OQ958044.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ958044.1). For the tree, sequences are aligned to this.

- [protein_sites.tsv](protein_sites.tsv): file giving the mapping of protein sites between 1, 2, ... numbering of the protein sequence encoded by [reference_sequence.fasta](reference_sequence.fasta) and other numbering schemes:
  - *sequential_site*: 1, 2, ... numbering of protein encoded by [reference_sequence.fasta](reference_sequence.fasta)
  - *H3_site*: equivalent number in H3 numbering
  - *H5_site*: equivalent number in mature H5 numbering (ectodomain)
  - *protein*: domain of HA: *SigPep*, *HA1*, *HA2*
  - *protein_site*: sequential 1, 2, ... numbering of *protein*

- [annotation.gff](annotation.gff): GFF3 annotation defining the gene regions (SigPep, HA1, HA2) with nucleotide coordinates in the reference sequence. Generated from [protein_sites.tsv](protein_sites.tsv) by running:
  ```
  python protein_sites_to_gff.py protein_sites.tsv reference_sequence.fa annotation.gff
  ```

