"""Collapse low-frequency host orders to 'other' in subsampled metadata."""

import sys

import pandas as pd


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    metadata_file = snakemake.input.metadata  # noqa: F821
    out_metadata = snakemake.output.metadata  # noqa: F821
    min_freq = snakemake.params.collapse_low_freq_host_order  # noqa: F821

    df = pd.read_csv(metadata_file, sep="\t", dtype=str)

    non_null = df["host_order"].notna()
    counts = df.loc[non_null, "host_order"].value_counts()
    freqs = counts / counts.sum()

    print(f"host_order frequencies before collapsing:\n{freqs.to_string()}\n")

    low_freq_orders = sorted(order for order, freq in freqs.items() if freq < min_freq)
    print(
        f"Collapsing {len(low_freq_orders)} orders with freq < {min_freq} to 'other':"
    )
    for order in low_freq_orders:
        print(f"  {order}: {freqs[order]:.4f} ({counts[order]} sequences)")

    df.loc[non_null & df["host_order"].isin(low_freq_orders), "host_order"] = "other"

    print(
        f"\nhost_order counts after collapsing:\n{df['host_order'].value_counts(dropna=False)}"
    )

    df.to_csv(out_metadata, sep="\t", index=False)
    print(f"Written metadata with {len(df)} rows")


if __name__ == "__main__":
    main()
