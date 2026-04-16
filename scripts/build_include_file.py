"""Combine user accessions_to_include.txt with manual-add accessions."""

import sys

import pandas as pd


def parse_accessions_file(path):
    """Read an accessions_to_include.txt file, stripping `#` comments and blanks."""
    accs = set()
    with open(path) as f:
        for line in f:
            token = line.split("#", 1)[0].strip()
            if token:
                accs.add(token)
    return accs


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    accessions_to_include = snakemake.input.accessions_to_include  # noqa: F821
    manual_add_metadata = snakemake.input.manual_add_metadata  # noqa: F821
    out_file = snakemake.output.include  # noqa: F821

    genbank_accs = parse_accessions_file(accessions_to_include)
    print(f"Genbank accessions_to_include: {len(genbank_accs)}")

    manual_df = pd.read_csv(manual_add_metadata, sep="\t", dtype=str)
    manual_accs = set(manual_df["accession"])
    print(f"Manual-add accessions: {len(manual_accs)}")

    combined = sorted(genbank_accs | manual_accs)
    with open(out_file, "w") as f:
        f.write("\n".join(combined) + ("\n" if combined else ""))
    print(f"Wrote {len(combined)} combined accessions to {out_file}")


if __name__ == "__main__":
    main()
