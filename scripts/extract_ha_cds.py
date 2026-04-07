"""Extract full-length HA CDS sequences and metadata for a specific subtype."""

import gzip
import re
import sys
import zipfile

import numpy as np
import pandas as pd
import yaml
from Bio.Seq import Seq

METADATA_COLUMNS = [
    "accession",
    "strain",
    "subtype",
    "date",
    "host",
    "host_tax_id",
    "location",
    "region",
    "passage_history",
    "length",
]

VALID_NUCS = set("ACTGactg")
DATE_RE = re.compile(r"^\d{4}(-\d{2}(-\d{2})?)?$")
HA_GENE_NAMES = {"HA", "ha", "HA1"}
STRAIN_YEAR_RE = re.compile(r"/(\d{4})(?:\(|$)")


def parse_subtypes_from_headers(headers_file):
    """Parse {accession: subtype} from genomic.fna header lines."""
    subtype_re = re.compile(r"\(A/.+?\(([^()]+)\)\)")
    acc_subtypes = {}
    with open(headers_file) as f:
        for line in f:
            acc = line.lstrip(">").split()[0]
            m = subtype_re.search(line)
            if m:
                subtype = m.group(1).upper()
                if "MIXED" not in subtype:
                    acc_subtypes[acc] = subtype
    return acc_subtypes


def filter_metadata(metadata_file, annotation_genes_file, acc_subtypes, subtype_regex):
    """Read genome metadata, filter to HA accessions matching subtype regex.

    HA accessions are identified by segment == "4" in genome metadata OR
    gene-name in HA_GENE_NAMES in the annotation report.

    Returns (filtered_df, stats_lines).
    """
    df = pd.read_csv(metadata_file, sep="\t", dtype=str)
    annot = pd.read_csv(annotation_genes_file, sep="\t", dtype=str)
    ha_gene_accs = set(annot.loc[annot["Gene Name"].isin(HA_GENE_NAMES), "Accession"])
    is_seg4 = df["Segment"] == "4"
    is_ha_gene = df["Accession"].isin(ha_gene_accs)
    df = df[is_seg4 | is_ha_gene].copy()
    df["subtype"] = df["Accession"].map(acc_subtypes)
    df = df.dropna(subset=["subtype"])

    stats = []
    n_seg4_only = int((is_seg4 & ~is_ha_gene).sum())
    n_gene_only = int((~is_seg4 & is_ha_gene).sum())
    n_both = int((is_seg4 & is_ha_gene).sum())
    stats.append(
        f"HA accessions (segment=4 OR gene-name in {HA_GENE_NAMES}): {len(df)}"
    )
    stats.append(f"  segment=4 only: {n_seg4_only}")
    stats.append(f"  gene-name only: {n_gene_only}")
    stats.append(f"  both: {n_both}")
    hxnx_re = re.compile(r"H\d{1,2}N\d{1,2}")
    n_with_subtype = len(df)
    n_hxnx = int(df["subtype"].str.fullmatch(hxnx_re).sum())
    n_tree = int(df["subtype"].str.fullmatch(subtype_regex).sum())
    n_non_hxnx = n_with_subtype - n_hxnx
    stats.append(f"HA accessions with parseable subtype: {n_with_subtype}")
    stats.append(f"  matching any HxNx pattern: {n_hxnx}")
    stats.append(f"  matching tree subtype regex '{subtype_regex}': {n_tree}")
    stats.append(f"  not matching any HxNx pattern: {n_non_hxnx}")

    df = df[df["subtype"].str.fullmatch(subtype_regex)]
    df = df.rename(
        columns={
            "Accession": "accession",
            "Isolate Lineage": "strain",
            "Isolate Collection date": "date",
            "Host Name": "host",
            "Host Taxonomic ID": "host_tax_id",
            "Geographic Location": "location",
            "Geographic Region": "region",
            "Lab Host": "passage_history",
            "Length": "length",
        }
    )
    return df[METADATA_COLUMNS], stats


def extract_cds_from_zip(zip_path, accepted_accessions):
    """Extract full-length hemagglutinin CDS sequences from zip."""
    sequences = {}
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open("ncbi_dataset/data/cds.fna") as f:
            header = None
            seq_lines = []
            for raw_line in f:
                line = raw_line.decode()
                if line.startswith(">"):
                    if header is not None:
                        _maybe_add_seq(
                            header, seq_lines, accepted_accessions, sequences
                        )
                    header = line.rstrip()
                    seq_lines = []
                else:
                    seq_lines.append(line.strip())
            if header is not None:
                _maybe_add_seq(header, seq_lines, accepted_accessions, sequences)
    return sequences


def _maybe_add_seq(header, seq_lines, accepted_accessions, sequences):
    """Add sequence if it's a full-length HA CDS for an accepted accession."""
    acc = header.lstrip(">").split(":")[0]
    if acc not in accepted_accessions:
        return
    # full CDS entries have name "hemagglutinin" without [polyprotein=...]
    parts = header.split(" ", 1)
    if len(parts) < 2:
        return
    defline = parts[1]
    if "[polyprotein=" in defline:
        return
    if not defline.startswith("hemagglutinin"):
        return
    seq = "".join(seq_lines)
    # keep longest if multiple CDS per accession
    if acc not in sequences or len(seq) > len(sequences[acc]):
        sequences[acc] = seq


def filter_sequences(sequences):
    """Filter sequences with ambiguous nucleotides or invalid CDS translation.

    Returns (filtered_sequences, stats_lines).
    """
    stats = []

    # filter ambiguous nucleotides
    clean = {acc: seq for acc, seq in sequences.items() if set(seq) <= VALID_NUCS}
    stats.append(f"Dropped for ambiguous nucleotides: {len(sequences) - len(clean)}")

    # filter by CDS validity: starts with ATG, length multiple of 3, no internal stops
    valid = {}
    for acc, seq in clean.items():
        if not seq.upper().startswith("ATG"):
            continue
        if len(seq) % 3 != 0:
            continue
        protein = str(Seq(seq).translate())
        # allow optional terminal stop but no internal stops
        if "*" in protein[:-1]:
            continue
        valid[acc] = seq
    stats.append(f"Dropped for invalid CDS: {len(clean) - len(valid)}")

    # length distribution
    lengths = np.array([len(seq) for seq in valid.values()])
    stats.append(f"CDS length distribution (n={len(valid)}):")
    stats.append(f"  min: {lengths.min()}, max: {lengths.max()}")
    stats.append(f"  0.1th percentile: {int(np.percentile(lengths, 0.1))}")
    stats.append(f"  1st percentile: {int(np.percentile(lengths, 1))}")
    stats.append(f"  10th percentile: {int(np.percentile(lengths, 10))}")
    stats.append(f"  median: {int(np.median(lengths))}")
    stats.append(f"  90th percentile: {int(np.percentile(lengths, 90))}")
    stats.append(f"  99th percentile: {int(np.percentile(lengths, 99))}")
    stats.append(f"  99.9th percentile: {int(np.percentile(lengths, 99.9))}")

    return valid, stats


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    metadata_file = snakemake.input.genome_metadata  # noqa: F821
    headers_file = snakemake.input.genomic_headers  # noqa: F821
    annotation_genes_file = snakemake.input.annotation_genes  # noqa: F821
    zip_path = snakemake.input.zipfile  # noqa: F821
    subtype_regex = snakemake.params.subtype_regex  # noqa: F821
    with open(snakemake.input.cds_length_range) as f:  # noqa: F821
        cds_length_range = yaml.safe_load(f)["cds_length_range"]
    if len(cds_length_range) != 2:
        raise ValueError(f"cds_length_range must have 2 elements: {cds_length_range}")
    min_len, max_len = cds_length_range
    out_metadata = snakemake.output.metadata  # noqa: F821
    out_fasta = snakemake.output.fasta  # noqa: F821
    out_stats = snakemake.output.stats  # noqa: F821

    acc_subtypes = parse_subtypes_from_headers(headers_file)
    metadata, stats = filter_metadata(
        metadata_file, annotation_genes_file, acc_subtypes, subtype_regex
    )
    accepted_accessions = set(metadata["accession"])
    stats.append(
        f"Metadata accessions matching subtype regex: {len(accepted_accessions)}"
    )

    sequences = extract_cds_from_zip(zip_path, accepted_accessions)
    if not sequences:
        raise RuntimeError(
            f"No HA CDS sequences found for subtype regex '{subtype_regex}'"
        )

    found = accepted_accessions & set(sequences)
    not_found = accepted_accessions - set(sequences)
    stats.append(f"Accessions with full-length HA CDS found: {len(found)}")
    stats.append(f"Accessions without full-length HA CDS: {len(not_found)}")

    # fill missing dates from strain name year as fallback, then validate
    metadata = metadata[metadata["accession"].isin(sequences)]
    no_date = metadata["date"].isna() | (metadata["date"] == "")
    n_no_date = int(no_date.sum())
    parsed_years = metadata.loc[no_date, "strain"].str.extract(
        STRAIN_YEAR_RE, expand=False
    )
    metadata.loc[no_date, "date"] = parsed_years
    n_recovered = int(parsed_years.notna().sum())
    stats.append(
        f"Recovered date from strain name: {n_recovered} of {n_no_date}"
        " missing-date accessions"
    )
    bad_dates = metadata[~metadata["date"].str.match(DATE_RE, na=True)]
    if len(bad_dates):
        raise ValueError(
            f"Dates not matching YYYY[-MM[-DD]]: "
            f"{bad_dates['date'].unique().tolist()}"
        )
    no_date = metadata["date"].isna() | (metadata["date"] == "")
    stats.append(f"Dropped for missing collection date: {no_date.sum()}")
    metadata = metadata[~no_date]
    sequences = {acc: sequences[acc] for acc in metadata["accession"]}

    sequences, seq_stats = filter_sequences(sequences)
    stats.extend(seq_stats)
    if not sequences:
        raise RuntimeError(
            f"No valid HA CDS sequences after filtering for '{subtype_regex}'"
        )

    # filter by CDS length range (either bound can be null for no limit)
    n_before = len(sequences)
    too_short = (
        {acc for acc, seq in sequences.items() if len(seq) < min_len}
        if min_len is not None
        else set()
    )
    too_long = (
        {acc for acc, seq in sequences.items() if len(seq) > max_len}
        if max_len is not None
        else set()
    )
    sequences = {
        acc: seq
        for acc, seq in sequences.items()
        if acc not in too_short and acc not in too_long
    }
    stats.append(f"CDS length range filter [{min_len}, {max_len}]:")
    stats.append(f"  dropped too short (<{min_len}): {len(too_short)}")
    stats.append(f"  dropped too long (>{max_len}): {len(too_long)}")
    stats.append(f"  kept: {len(sequences)} of {n_before}")
    if not sequences:
        raise RuntimeError(f"No sequences in length range [{min_len}, {max_len}]")

    # only keep metadata rows that have a CDS sequence
    metadata = metadata[metadata["accession"].isin(sequences)]
    metadata = metadata.sort_values(["date", "strain"]).reset_index(drop=True)
    metadata.to_csv(out_metadata, sep="\t", index=False)

    with gzip.GzipFile(out_fasta, "wb", mtime=0) as gz:
        for acc in metadata["accession"]:
            gz.write(f">{acc}\n{sequences[acc]}\n".encode())

    stats.append(f"Written {len(sequences)} HA CDS sequences for '{subtype_regex}'")

    for line in stats:
        print(line)

    with open(out_stats, "w") as f:
        f.write(
            f"All HA CDSs extracted from NCBI metadata for subtype regex"
            f" '{subtype_regex}' (before subsampling)\n\n"
        )
        for line in stats:
            f.write(line + "\n")


if __name__ == "__main__":
    main()
