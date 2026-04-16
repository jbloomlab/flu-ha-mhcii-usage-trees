"""Merge manual-add CDS sequences and metadata into Genbank-extracted data."""

import gzip
import re
import sys

import pandas as pd
from Bio import SeqIO

from extract_ha_cds import METADATA_COLUMNS

# accession pattern required by the Auspice JSON schema (matches augur export v2)
ACCESSION_RE = re.compile(r"^[0-9A-Za-z\-_.]+$")


def read_fasta(path, open_fn=open):
    """Read FASTA returning {accession: sequence}, using `open_fn` to open the file."""
    with open_fn(path, "rt") as f:
        return {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}


def validate_manual_add(manual_meta, manual_seqs):
    """Validate manual-add metadata and sequences. Returns the set of accessions."""
    if list(manual_meta.columns) != METADATA_COLUMNS:
        raise ValueError(
            f"Manual-add metadata columns {list(manual_meta.columns)} do not "
            f"match expected {METADATA_COLUMNS}"
        )

    accs = manual_meta["accession"].tolist()
    if len(accs) != len(set(accs)):
        dups = [a for a in set(accs) if accs.count(a) > 1]
        raise ValueError(f"Duplicate accessions in manual-add metadata: {dups}")

    bad_chars = [a for a in accs if not ACCESSION_RE.fullmatch(str(a))]
    if bad_chars:
        raise ValueError(
            f"Manual-add accessions must match {ACCESSION_RE.pattern} "
            f"(required by Auspice JSON schema). Offenders: {bad_chars}"
        )

    meta_accs = set(accs)
    seq_accs = set(manual_seqs)
    if meta_accs != seq_accs:
        raise ValueError(
            f"Manual-add metadata and FASTA accessions do not match.\n"
            f"  In metadata only: {sorted(meta_accs - seq_accs)}\n"
            f"  In FASTA only:    {sorted(seq_accs - meta_accs)}"
        )

    return meta_accs


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    genbank_metadata = snakemake.input.genbank_metadata  # noqa: F821
    genbank_sequences = snakemake.input.genbank_sequences  # noqa: F821
    manual_metadata = snakemake.input.manual_add_metadata  # noqa: F821
    manual_sequences = snakemake.input.manual_add_sequences  # noqa: F821
    out_metadata = snakemake.output.metadata  # noqa: F821
    out_sequences = snakemake.output.sequences  # noqa: F821

    manual_meta = pd.read_csv(manual_metadata, sep="\t", dtype=str)
    manual_seqs = read_fasta(manual_sequences)
    manual_accs = validate_manual_add(manual_meta, manual_seqs)
    print(f"Manual-add: {len(manual_accs)} accessions")

    genbank_meta = pd.read_csv(genbank_metadata, sep="\t", dtype=str)
    if list(genbank_meta.columns) != METADATA_COLUMNS:
        raise ValueError(
            f"Genbank metadata columns {list(genbank_meta.columns)} do not "
            f"match expected {METADATA_COLUMNS}"
        )
    genbank_seqs = read_fasta(genbank_sequences, open_fn=gzip.open)
    if set(genbank_meta["accession"]) != set(genbank_seqs):
        raise ValueError("Genbank metadata and FASTA accessions do not match")
    print(f"Genbank: {len(genbank_seqs)} accessions")

    overlap = manual_accs & set(genbank_meta["accession"])
    if overlap:
        print(
            f"Dropping {len(overlap)} Genbank accession(s) that overlap with "
            f"manual-add: {sorted(overlap)}"
        )
        genbank_meta = genbank_meta[~genbank_meta["accession"].isin(overlap)]
        for a in overlap:
            del genbank_seqs[a]

    merged_meta = pd.concat([genbank_meta, manual_meta], ignore_index=True)
    merged_meta.to_csv(out_metadata, sep="\t", index=False)
    print(f"Wrote merged metadata with {len(merged_meta)} rows to {out_metadata}")

    with gzip.open(out_sequences, "wt") as f:
        for acc in merged_meta["accession"]:
            seq = genbank_seqs[acc] if acc in genbank_seqs else manual_seqs[acc]
            f.write(f">{acc}\n{seq}\n")
    print(f"Wrote merged FASTA with {len(merged_meta)} sequences to {out_sequences}")


if __name__ == "__main__":
    main()
