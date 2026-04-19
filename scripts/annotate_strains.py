"""Annotate tree tips from a per-strain annotations TSV.

Emits a node-data JSON with per-tip annotations and an auspice config JSON
with a coloring for every annotation column (continuous for numeric columns
using the shared `phenotype_color_scale`, categorical otherwise). Non-numeric
columns are also added to `filters`.
"""

import json
import sys

import pandas as pd
from Bio import Phylo


ID_COL = "accession"


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    tsv_path = snakemake.input.strain_annotations  # noqa: F821
    tree_path = snakemake.input.tree  # noqa: F821
    node_data_path = snakemake.output.node_data  # noqa: F821
    auspice_config_path = snakemake.output.auspice_config  # noqa: F821
    action = snakemake.params.missing_strain_annotations_action  # noqa: F821
    color_scale = list(snakemake.params.phenotype_color_scale)  # noqa: F821

    if action not in ("ignore", "error"):
        raise ValueError(
            f"Invalid missing_strain_annotations_action={action!r}, "
            "must be 'ignore' or 'error'"
        )

    # Empty-file short-circuit: file has no parseable content (0-byte or
    # whitespace-only) or a header with no data rows.
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
        is_empty = len(df) == 0
    except pd.errors.EmptyDataError:
        df = None
        is_empty = True
    if is_empty:
        print(f"{tsv_path} has no data; writing empty node data and auspice config.")
        with open(node_data_path, "w") as f:
            json.dump({"nodes": {}}, f, indent=2)
        with open(auspice_config_path, "w") as f:
            json.dump({"colorings": [], "filters": []}, f, indent=2)
        return

    if ID_COL not in df.columns:
        raise ValueError(
            f"{tsv_path} has no {ID_COL!r} column; got columns {list(df.columns)}"
        )

    accessions = df[ID_COL]
    if accessions.isna().any() or (accessions.str.strip() == "").any():
        raise ValueError(f"{tsv_path} has empty/null values in {ID_COL!r} column")
    if accessions.str.contains(r"\s").any():
        raise ValueError(f"{tsv_path} has whitespace in {ID_COL!r} values")
    dupes = sorted(accessions[accessions.duplicated()].unique())
    if dupes:
        raise ValueError(f"{tsv_path} has duplicate {ID_COL!r} values: {dupes}")

    tree = Phylo.read(tree_path, "newick")
    tips = {t.name for t in tree.get_terminals()}

    missing = sorted(set(accessions) - tips)
    if missing:
        if action == "error":
            raise ValueError(
                f"{len(missing)} accessions in {tsv_path} are not tips in "
                f"{tree_path} and missing_strain_annotations_action='error': "
                f"{missing}"
            )
        print(
            f"WARNING: dropping {len(missing)} accession(s) not present as tips "
            f"in {tree_path}: {missing}"
        )
        df = df[df[ID_COL].isin(tips)].reset_index(drop=True)

    annotation_cols = [c for c in df.columns if c != ID_COL]
    if not annotation_cols:
        raise ValueError(
            f"{tsv_path} has only {ID_COL!r} column and no annotation columns"
        )

    # Detect numeric columns: cast a copy to numeric, numeric iff no non-null
    # value fails to parse.
    numeric_cols = []
    numeric_values = {}
    for col in annotation_cols:
        coerced = pd.to_numeric(df[col], errors="coerce")
        non_null = df[col].notna() & (df[col].astype(str).str.strip() != "")
        if non_null.any() and coerced[non_null].notna().all():
            numeric_cols.append(col)
            numeric_values[col] = coerced

    nodes = {}
    for _, row in df.iterrows():
        acc = row[ID_COL]
        entry = {}
        for col in annotation_cols:
            raw = row[col]
            if pd.isna(raw) or (isinstance(raw, str) and raw.strip() == ""):
                continue
            if col in numeric_cols:
                entry[col] = float(numeric_values[col].loc[row.name])
            else:
                entry[col] = str(raw)
        if entry:
            nodes[acc] = entry

    colorings = []
    filters = []
    for col in annotation_cols:
        if col in numeric_cols:
            vals = numeric_values[col].dropna()
            vmin, vmax = float(vals.min()), float(vals.max())
            num_colors = len(color_scale)
            if vmax == vmin or num_colors == 1:
                scale = [[vmin, color_scale[0]]]
            else:
                scale = [
                    [vmin + (vmax - vmin) * i / (num_colors - 1), color]
                    for i, color in enumerate(color_scale)
                ]
            colorings.append(
                {"key": col, "title": col, "type": "continuous", "scale": scale}
            )
        else:
            colorings.append({"key": col, "title": col, "type": "categorical"})
            filters.append(col)

    with open(node_data_path, "w") as f:
        json.dump({"nodes": nodes}, f, indent=2)
    with open(auspice_config_path, "w") as f:
        json.dump({"colorings": colorings, "filters": filters}, f, indent=2)
    print(
        f"Wrote annotations for {len(nodes)} tips across "
        f"{len(annotation_cols)} column(s) "
        f"({len(numeric_cols)} numeric, {len(annotation_cols) - len(numeric_cols)} "
        f"categorical)."
    )


if __name__ == "__main__":
    main()
