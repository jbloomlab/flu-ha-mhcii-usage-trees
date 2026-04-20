"""Infer host from strain name for rows where NCBI host info is missing.

For each metadata row where both `host` and `host_tax_id` are empty, parse the
2nd `/`-delimited token from the strain name (lowercased and whitespace-stripped)
and look it up in `strain_token_host_map.tsv`. If the map has a non-empty
`host_tax_id` for the token, fill in `host` and `host_tax_id`. Otherwise (token
missing or explicitly flagged skip) leave both fields empty. Unknown tokens are
logged as warnings so the map can be extended over time.
"""

import sys

import pandas as pd

from extract_ha_cds import METADATA_COLUMNS

MAP_COLUMNS = ["token", "host_tax_id", "host_scientific_name", "notes"]


def parse_strain_token(strain):
    """Return the lowercased, stripped 2nd `/`-delimited token, or None."""
    if not isinstance(strain, str):
        return None
    parts = strain.split("/")
    if len(parts) < 4:
        return None
    return parts[1].strip().lower()


def load_map(path):
    """Load and validate the strain-token map. Returns {token: (tax_id, sci_name)}."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if list(df.columns) != MAP_COLUMNS:
        raise ValueError(
            f"Map columns {list(df.columns)} do not match expected {MAP_COLUMNS}"
        )
    if df["token"].duplicated().any():
        dups = sorted(df.loc[df["token"].duplicated(), "token"].unique())
        raise ValueError(f"Duplicate tokens in map: {dups}")
    if (df["token"] != df["token"].str.strip().str.lower()).any():
        raise ValueError("Map tokens must be whitespace-stripped and lowercase")

    mapping = {}
    for _, row in df.iterrows():
        tid = row["host_tax_id"]
        name = row["host_scientific_name"]
        if (tid == "") != (name == ""):
            raise ValueError(
                f"Row for token {row['token']!r} must have both host_tax_id and "
                f"host_scientific_name set, or both blank (got "
                f"tax_id={tid!r}, name={name!r})"
            )
        mapping[row["token"]] = (tid, name)
    return mapping


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    in_meta = snakemake.input.metadata  # noqa: F821
    in_map = snakemake.input.strain_token_host_map  # noqa: F821
    out_meta = snakemake.output.metadata  # noqa: F821

    df = pd.read_csv(in_meta, sep="\t", dtype=str)
    if list(df.columns) != METADATA_COLUMNS:
        raise ValueError(
            f"Metadata columns {list(df.columns)} do not match expected "
            f"{METADATA_COLUMNS}"
        )

    mapping = load_map(in_map)
    print(f"Loaded {len(mapping)} tokens from {in_map}")

    missing_mask = df["host"].isna() & df["host_tax_id"].isna()
    n_missing_before = int(missing_mask.sum())
    print(f"Rows missing host info: {n_missing_before} / {len(df)}")

    filled = 0
    skipped_explicit = 0
    unknown_tokens = {}  # token -> count
    unparseable_strains = 0

    for idx in df.index[missing_mask]:
        token = parse_strain_token(df.at[idx, "strain"])
        if token is None:
            unparseable_strains += 1
            continue
        if token not in mapping:
            unknown_tokens[token] = unknown_tokens.get(token, 0) + 1
            continue
        tid, name = mapping[token]
        if tid == "":
            skipped_explicit += 1
            continue
        df.at[idx, "host"] = name
        df.at[idx, "host_tax_id"] = tid
        filled += 1

    remaining_missing = n_missing_before - filled

    print(f"Filled host for {filled} rows")
    print(f"Explicit-skip tokens (blank tax_id in map): {skipped_explicit} rows")
    print(f"Unparseable strain names: {unparseable_strains} rows")
    print(f"Unknown tokens encountered: {len(unknown_tokens)} unique")
    if unknown_tokens:
        print("Unknown tokens (consider adding to strain_token_host_map.tsv):")
        for tok, n in sorted(unknown_tokens.items(), key=lambda x: (-x[1], x[0])):
            print(f"  {n:>4}  {tok!r}")
    print(f"Rows still missing host info: {remaining_missing} / {len(df)}")

    df.to_csv(out_meta, sep="\t", index=False)
    print(f"Wrote {out_meta}")


if __name__ == "__main__":
    main()
