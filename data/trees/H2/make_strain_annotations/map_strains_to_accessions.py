"""Map strain names in ``log10_titers_relative_to.csv`` to GenBank accessions.

Run manually (not part of the Snakemake pipeline) from this directory:

    python map_strains_to_accessions.py

First consults ``data/trees/H5/accessions_to_include.txt`` (parsing the
``accession  # strain`` comments, including ``strain_a = strain_b`` aliases),
then falls back to ``data/trees/H5/manual_add_metadata.tsv``.

Writes the per-strain titer table augmented with an ``accession`` column
(``strain`` column dropped, unmapped rows dropped) to
``data/trees/H5/accession_log10_titers_relative_to.tsv`` and a human-readable
summary to ``map_strains_to_accessions_summary.txt`` in this directory.
"""

import os
import re

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
H5_DIR = os.path.dirname(SCRIPT_DIR)

TITERS_CSV = os.path.join(SCRIPT_DIR, "log10_titers_relative_to.csv")
ACCESSIONS_TO_INCLUDE = os.path.join(H5_DIR, "accessions_to_include.txt")
MANUAL_ADD_METADATA = os.path.join(H5_DIR, "manual_add_metadata.tsv")

OUT_TSV = os.path.join(H5_DIR, "strain_annotations.tsv")
OUT_SUMMARY = os.path.join(SCRIPT_DIR, "map_strains_to_accessions_summary.txt")

STRAINS_TO_EXCLUDE = [
    "A/Netherlands/B1/1968",  # experimental data may have high noSA background for this strain
]


def normalize(strain):
    """Lowercase, convert underscores to spaces, collapse whitespace."""
    return re.sub(r"\s+", " ", strain.replace("_", " ")).strip().lower()


def parse_accessions_to_include(path):
    """Yield (accession, [strain_aliases]) for each non-blank line.

    Each non-comment line must be ``accession  # strain [= alias ...]``.
    """
    with open(path) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if "#" not in line:
                raise ValueError(
                    f"{path} line {lineno}: expected '# strain' comment: {line!r}"
                )
            acc_part, comment = line.split("#", 1)
            accession = acc_part.strip()
            if not accession:
                raise ValueError(f"{path} line {lineno}: empty accession: {line!r}")
            aliases = [s.strip() for s in comment.split("=")]
            aliases = [s for s in aliases if s]
            if not aliases:
                raise ValueError(
                    f"{path} line {lineno}: no strain name in comment: {line!r}"
                )
            yield accession, aliases


def build_mapping():
    """Return (norm_strain -> (accession, source, original_name), collisions)."""
    mapping = {}
    collisions = []

    for accession, aliases in parse_accessions_to_include(ACCESSIONS_TO_INCLUDE):
        for alias in aliases:
            key = normalize(alias)
            if key in mapping:
                if mapping[key][0] != accession:
                    collisions.append(
                        (alias, mapping[key][0], "accessions_to_include", accession)
                    )
                continue
            mapping[key] = (accession, "accessions_to_include", alias)

    manual_df = pd.read_csv(MANUAL_ADD_METADATA, sep="\t", dtype=str)
    for col in ("accession", "strain"):
        if col not in manual_df.columns:
            raise ValueError(f"{MANUAL_ADD_METADATA} missing column {col!r}")
    for _, row in manual_df.iterrows():
        key = normalize(row["strain"])
        if key in mapping:
            if mapping[key][0] != row["accession"]:
                collisions.append(
                    (
                        row["strain"],
                        mapping[key][0],
                        "manual_add_metadata",
                        row["accession"],
                    )
                )
            continue
        mapping[key] = (row["accession"], "manual_add_metadata", row["strain"])

    return mapping, collisions


def main():
    mapping, collisions = build_mapping()

    titers = pd.read_csv(TITERS_CSV)
    if "strain" not in titers.columns:
        raise ValueError(f"{TITERS_CSV} missing 'strain' column")

    titers = titers[~titers["strain"].isin(STRAINS_TO_EXCLUDE)]

    accessions = []
    sources = []
    matched_names = []
    for strain in titers["strain"]:
        hit = mapping.get(normalize(strain))
        if hit is None:
            accessions.append(None)
            sources.append(None)
            matched_names.append(None)
        else:
            accessions.append(hit[0])
            sources.append(hit[1])
            matched_names.append(hit[2])

    titers = titers.assign(
        accession=accessions, _source=sources, _matched=matched_names
    )

    mapped = titers[titers["accession"].notna()].copy()
    unmapped = titers[titers["accession"].isna()].copy()

    out = mapped.drop(columns=["strain", "_source", "_matched"])
    out = out[["accession"] + [c for c in out.columns if c != "accession"]]
    out.to_csv(OUT_TSV, sep="\t", index=False)

    src_counts = mapped["_source"].value_counts().to_dict()
    lines = []
    lines.append(f"Total strains in {os.path.basename(TITERS_CSV)}: {len(titers)}")
    lines.append(f"Mapped: {len(mapped)}")
    for src in ("accessions_to_include", "manual_add_metadata"):
        lines.append(f"  {src}: {src_counts.get(src, 0)}")
    lines.append(f"Unmapped: {len(unmapped)}")
    lines.append("")

    lines.append("Mappings (strain -> accession [source; matched name]):")
    for _, row in mapped.sort_values("strain").iterrows():
        lines.append(
            f"  {row['strain']} -> {row['accession']} "
            f"[{row['_source']}; matched {row['_matched']!r}]"
        )
    lines.append("")

    lines.append(f"Unmapped strains ({len(unmapped)}):")
    for strain in unmapped["strain"].sort_values():
        lines.append(f"  {strain}")
    lines.append("")

    lines.append(f"Cross-source collisions ({len(collisions)}):")
    for strain, existing_acc, new_source, new_acc in collisions:
        lines.append(
            f"  {strain}: kept {existing_acc}, ignored {new_acc} from {new_source}"
        )

    with open(OUT_SUMMARY, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Wrote {len(out)} mapped rows to {OUT_TSV}")
    print(f"Wrote summary to {OUT_SUMMARY}")
    print(f"Unmapped strains: {len(unmapped)}")


if __name__ == "__main__":
    main()
