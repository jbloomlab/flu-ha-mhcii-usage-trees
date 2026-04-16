"""Validate that every accession in the include list is a tip in the Auspice JSON."""

import json
import sys

import pandas as pd

from build_include_file import parse_accessions_file


def collect_leaf_names(node, leaves):
    """Walk the tree, collecting `name` values of all leaf nodes (no children)."""
    if "children" not in node or not node["children"]:
        leaves.add(node["name"])
        return
    for child in node["children"]:
        collect_leaf_names(child, leaves)


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    auspice_json = snakemake.input.auspice_json  # noqa: F821
    accessions_to_include = snakemake.input.accessions_to_include  # noqa: F821
    manual_add_metadata = snakemake.input.manual_add_metadata  # noqa: F821
    out_tsv = snakemake.output.tsv  # noqa: F821

    genbank_accs = parse_accessions_file(accessions_to_include)
    manual_accs = set(
        pd.read_csv(manual_add_metadata, sep="\t", dtype=str)["accession"]
    )
    # manual_add wins on overlap
    expected = {acc: "genbank" for acc in genbank_accs}
    expected.update({acc: "manual_add" for acc in manual_accs})
    print(
        f"Expected accessions: {len(expected)} "
        f"(genbank: {len(genbank_accs)}, manual_add: {len(manual_accs)}, "
        f"overlap: {len(genbank_accs & manual_accs)})"
    )

    with open(auspice_json) as f:
        tree = json.load(f)["tree"]
    leaves = set()
    collect_leaf_names(tree, leaves)
    print(f"Leaves in Auspice JSON: {len(leaves)}")

    rows = [
        {"accession": acc, "source": src, "in_tree": acc in leaves}
        for acc, src in expected.items()
    ]
    df = pd.DataFrame(rows, columns=["accession", "source", "in_tree"])
    df = df.sort_values(["in_tree", "source", "accession"])
    df.to_csv(out_tsv, sep="\t", index=False)

    missing = df[~df["in_tree"]]
    print(f"Missing from tree: {len(missing)}")
    if len(missing):
        for source, group in missing.groupby("source"):
            print(f"  {source}: {sorted(group['accession'])}")
        raise ValueError(
            f"{len(missing)} of {len(expected)} include accessions are missing "
            f"from {auspice_json}; see {out_tsv} and log."
        )
    print("All include accessions are present in the tree.")


if __name__ == "__main__":
    main()
