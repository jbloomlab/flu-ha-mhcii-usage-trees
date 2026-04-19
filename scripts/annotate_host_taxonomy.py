"""Annotate metadata with host taxonomy classification using taxonkit."""

import subprocess
import sys

import pandas as pd

HUMAN_TAX_ID = "9606"
AVES_CLASS = "Aves"
MAMMALIA_CLASS = "Mammalia"
CARNIVORA_ORDER = "Carnivora"
PERISSODACTYLA_ORDER = "Perissodactyla"
SUS_GENUS = "Sus"
BOS_GENUS = "Bos"


def get_taxonomy_for_ids(tax_ids, taxonomy_dir):
    """Get lineage, class, order, and genus for each tax ID using taxonkit.

    Returns dict of {tax_id: (lineage, class_name, order_name, genus_name)}.
    """
    input_text = "\n".join(tax_ids)
    result = subprocess.run(
        ["taxonkit", "lineage", "--data-dir", taxonomy_dir],
        input=input_text,
        capture_output=True,
        text=True,
        check=True,
    )
    result2 = subprocess.run(
        ["taxonkit", "reformat", "--data-dir", taxonomy_dir, "-f", "{c}\t{o}\t{g}"],
        input=result.stdout,
        capture_output=True,
        text=True,
        check=True,
    )
    taxonomy = {}
    for line in result2.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        tax_id = parts[0]
        lineage = parts[1] if len(parts) >= 2 else ""
        class_name = parts[-3] if len(parts) >= 5 else ""
        order_name = parts[-2] if len(parts) >= 5 else ""
        genus_name = parts[-1] if len(parts) >= 5 else ""
        taxonomy[tax_id] = (lineage, class_name, order_name, genus_name)
    return taxonomy


def classify_host_general(tax_id, taxonomy):
    """Classify a tax ID as human, avian, swine, bovine, equine, carnivore, or other mammal."""
    if pd.isna(tax_id) or str(tax_id).strip() == "":
        return pd.NA
    tax_id_str = str(int(float(tax_id)))
    if tax_id_str not in taxonomy:
        raise ValueError(f"No taxonomy found for tax ID {tax_id_str}")
    lineage, class_name, order_name, genus_name = taxonomy[tax_id_str]
    # deleted tax IDs have empty lineage; treat as unknown
    if not lineage:
        print(
            f"WARNING: tax ID {tax_id_str} has no lineage (deleted?), treating as unknown"
        )
        return pd.NA
    if tax_id_str == HUMAN_TAX_ID:
        return "human"
    # use class from reformat, fall back to checking lineage string
    if class_name == AVES_CLASS or (not class_name and AVES_CLASS in lineage):
        return "avian"
    if class_name == MAMMALIA_CLASS or (not class_name and MAMMALIA_CLASS in lineage):
        if order_name == CARNIVORA_ORDER:
            return "carnivore"
        if order_name == PERISSODACTYLA_ORDER:
            return "equine"
        if genus_name == SUS_GENUS:
            return "swine"
        if genus_name == BOS_GENUS:
            return "bovine"
        return "other mammal"
    raise ValueError(
        f"Tax ID {tax_id_str} has class '{class_name}' and lineage "
        f"'{lineage}', expected '{AVES_CLASS}' or '{MAMMALIA_CLASS}'"
    )


def get_host_order(tax_id, taxonomy):
    """Get the order for a tax ID."""
    if pd.isna(tax_id) or str(tax_id).strip() == "":
        return pd.NA
    tax_id_str = str(int(float(tax_id)))
    lineage, _, order_name, _ = taxonomy[tax_id_str]
    if not lineage:
        return pd.NA
    return order_name if order_name else pd.NA


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    metadata_file = snakemake.input.metadata  # noqa: F821
    taxonomy_dir = snakemake.input.taxonomy_dir  # noqa: F821
    out_metadata = snakemake.output.metadata  # noqa: F821

    df = pd.read_csv(metadata_file, sep="\t", dtype=str)

    # get unique non-null tax IDs
    tax_ids = df["host_tax_id"].dropna()
    tax_ids = tax_ids[tax_ids.str.strip() != ""]
    unique_ids = sorted(set(str(int(float(x))) for x in tax_ids))
    print(f"Looking up taxonomy for {len(unique_ids)} unique host tax IDs")

    taxonomy = get_taxonomy_for_ids(unique_ids, taxonomy_dir)

    # rename host -> host_specific
    df = df.rename(columns={"host": "host_specific"})

    # add host_general and host_order
    df["host_general"] = df["host_tax_id"].apply(
        lambda x: classify_host_general(x, taxonomy)
    )
    df["host_order"] = df["host_tax_id"].apply(lambda x: get_host_order(x, taxonomy))

    # reorder columns
    cols = [
        "accession",
        "strain",
        "subtype",
        "date",
        "host_specific",
        "host_tax_id",
        "host_general",
        "host_order",
        "location",
        "region",
        "passage_history",
        "length",
    ]
    df = df[cols]

    print(f"host_general counts:\n{df['host_general'].value_counts(dropna=False)}")
    print(f"host_order counts:\n{df['host_order'].value_counts(dropna=False)}")

    df.to_csv(out_metadata, sep="\t", index=False)
    print(f"Written annotated metadata with {len(df)} rows")


if __name__ == "__main__":
    main()
