"""Prune an Auspice JSON to tips whose strain annotations match a filter.

Collapses unary internal nodes after pruning by merging their branch mutations
into the surviving child's branch so no mutation info is lost.
"""

import json
import sys

import pandas as pd


ID_COL = "accession"


def merge_branch_mutations(dst_child, src_node):
    """Prepend `src_node`'s branch mutations onto `dst_child`'s branch mutations."""
    src_muts = src_node.get("branch_attrs", {}).get("mutations", {})
    if not src_muts:
        return
    dst_branch = dst_child.setdefault("branch_attrs", {})
    dst_muts = dst_branch.setdefault("mutations", {})
    for gene, muts in src_muts.items():
        dst_muts[gene] = list(muts) + list(dst_muts.get(gene, []))


def prune(node, keep):
    """Return pruned subtree (possibly collapsed), or None if empty."""
    children = node.get("children")
    if not children:
        return node if node.get("name") in keep else None
    kept = [c for c in (prune(c, keep) for c in children) if c is not None]
    if not kept:
        return None
    if len(kept) == 1:
        child = kept[0]
        merge_branch_mutations(child, node)
        return child
    node["children"] = kept
    return node


def collect_tip_names(node, names):
    children = node.get("children")
    if not children:
        names.add(node.get("name"))
        return
    for c in children:
        collect_tip_names(c, names)


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

    source_json_path = snakemake.input.source_json  # noqa: F821
    strain_annotations_path = snakemake.input.strain_annotations  # noqa: F821
    out_json_path = snakemake.output.auspice_json  # noqa: F821
    keep_where = dict(snakemake.params.keep_where)  # noqa: F821
    title = snakemake.params.title  # noqa: F821

    if not keep_where:
        raise ValueError("keep_where is empty; must specify at least one column/value")

    with open(source_json_path) as f:
        auspice = json.load(f)

    df = pd.read_csv(strain_annotations_path, sep="\t", dtype=str)
    if ID_COL not in df.columns:
        raise ValueError(
            f"{strain_annotations_path} has no {ID_COL!r} column; "
            f"got {list(df.columns)}"
        )
    missing_cols = [c for c in keep_where if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"keep_where columns {missing_cols} not in {strain_annotations_path} "
            f"columns {list(df.columns)}"
        )

    mask = pd.Series(True, index=df.index)
    for col, value in keep_where.items():
        mask &= df[col].astype(str) == str(value)
    keep_set = set(df.loc[mask, ID_COL])
    if not keep_set:
        raise ValueError(f"No rows in {strain_annotations_path} satisfy {keep_where}")
    print(f"Keeping {len(keep_set)} accession(s) matching {keep_where}")

    source_tips = set()
    collect_tip_names(auspice["tree"], source_tips)
    missing_tips = sorted(keep_set - source_tips)
    if missing_tips:
        raise ValueError(
            f"{len(missing_tips)} accession(s) satisfy {keep_where} but are not "
            f"tips in {source_json_path}: {missing_tips}"
        )

    pruned_root = prune(auspice["tree"], keep_set)
    if pruned_root is None:
        raise ValueError("Pruning removed the entire tree")
    auspice["tree"] = pruned_root

    if title is not None:
        auspice.setdefault("meta", {})["title"] = title
        print(f"Set meta.title to {title!r}")

    final_tips = set()
    collect_tip_names(auspice["tree"], final_tips)
    if final_tips != keep_set:
        extra = sorted(final_tips - keep_set)
        dropped = sorted(keep_set - final_tips)
        raise ValueError(
            f"Unexpected pruning result. Extra tips: {extra[:5]}. "
            f"Dropped: {dropped[:5]}."
        )
    print(f"Pruned tree has {len(final_tips)} tips")

    with open(out_json_path, "w") as f:
        json.dump(auspice, f)


if __name__ == "__main__":
    main()
