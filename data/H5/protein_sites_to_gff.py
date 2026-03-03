"""Convert protein_sites.tsv to GFF3 annotation for augur ancestral --annotation.

Usage:
    python protein_sites_to_gff.py protein_sites.tsv reference_sequence.fa annotation.gff
"""

import argparse

import pandas as pd
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protein_sites_tsv", help="Input protein_sites.tsv")
    parser.add_argument("reference_fasta", help="Reference CDS FASTA")
    parser.add_argument("output_gff", help="Output GFF3 file")
    args = parser.parse_args()

    ref = SeqIO.read(args.reference_fasta, "fasta")
    sites = pd.read_csv(args.protein_sites_tsv, sep="\t")

    # validate sequential_site is 1, 2, 3, ...
    expected = list(range(1, len(sites) + 1))
    if sites["sequential_site"].tolist() != expected:
        raise ValueError("sequential_site is not sequential 1, 2, 3, ...")

    # validate total sites matches reference length
    if len(sites) * 3 != len(ref.seq):
        raise ValueError(
            f"{len(sites)} sites * 3 = {len(sites) * 3} nt, "
            f"but reference is {len(ref.seq)} nt"
        )

    # validate protein_site is sequential 1, 2, 3, ... within each protein
    for protein, group in sites.groupby("protein", sort=False):
        expected_sites = list(range(1, len(group) + 1))
        if group["protein_site"].tolist() != expected_sites:
            raise ValueError(
                f"protein_site for {protein} is not sequential 1, 2, 3, ..."
            )

    # build GFF: one gene feature per protein
    proteins = (
        sites.groupby("protein", sort=False)["sequential_site"]
        .agg(["first", "last"])
        .rename(columns={"first": "start_aa", "last": "end_aa"})
    )
    proteins["nt_start"] = (proteins["start_aa"] - 1) * 3 + 1
    proteins["nt_end"] = proteins["end_aa"] * 3

    with open(args.output_gff, "w") as f:
        f.write("##gff-version 3\n")
        f.write(f"##sequence-region {ref.id} 1 {len(ref.seq)}\n")
        for protein_name, row in proteins.iterrows():
            f.write(
                f"{ref.id}\t.\tgene\t{row['nt_start']}\t{row['nt_end']}\t.\t+\t.\t"
                f"gene_name={protein_name}\n"
            )

    print(f"Wrote {len(proteins)} gene features to {args.output_gff}:")
    for protein_name, row in proteins.iterrows():
        n_aa = row["end_aa"] - row["start_aa"] + 1
        print(f"  {protein_name}: {n_aa} aa, nt {row['nt_start']}-{row['nt_end']}")


if __name__ == "__main__":
    main()
