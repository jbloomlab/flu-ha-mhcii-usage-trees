#!/usr/bin/env python3
"""Score tree nodes based on per-mutation phenotype effects.

Computes cumulative amino acid mutations from root to each node, then scores
each node for each phenotype using measured mutation effects.
"""

import argparse
import json
import sys

import pandas as pd
from Bio import Phylo

sys.setrecursionlimit(10000)


class MutationScorer:
    """Score lists of mutations based on measured effects.

    Parameters
    ----------
    mutation_effects_df : pandas.DataFrame
        Must have columns `protein`, `protein_site`, `wildtype`, `mutant`, `effect`.
    """

    def __init__(self, mutation_effects_df):
        df = mutation_effects_df[
            ["protein", "protein_site", "wildtype", "mutant", "effect"]
        ].dropna()
        assert len(df) == len(
            df[["protein", "protein_site", "mutant"]].drop_duplicates()
        ), "Duplicate (protein, protein_site, mutant) entries in mutation effects"
        assert len(df[["protein", "protein_site"]].drop_duplicates()) == len(
            df[["protein", "protein_site", "wildtype"]].drop_duplicates()
        ), "Inconsistent wildtype for same (protein, protein_site)"
        self.wts = (
            df.drop_duplicates(subset=["protein", "protein_site"])
            .set_index(["protein", "protein_site"])["wildtype"]
            .to_dict()
        )
        self.effects = {}
        for (protein, site), site_df in df.groupby(["protein", "protein_site"]):
            self.effects[(protein, site)] = site_df.set_index("mutant")[
                "effect"
            ].to_dict()
        for (protein, site), wt in self.wts.items():
            assert (
                wt not in self.effects[(protein, site)]
                or self.effects[(protein, site)][wt] == 0
            )
            self.effects[(protein, site)][wt] = 0.0

    def total_effect(self, muts):
        """Total effect for mutations as (protein, wt, site, mutant) tuples."""
        total = 0.0
        for protein, wt, site, m in muts:
            key = (protein, site)
            if (
                key in self.effects
                and wt in self.effects[key]
                and m in self.effects[key]
            ):
                total += self.effects[key][m] - self.effects[key][wt]
        return total

    def max_magnitude_effect(self, muts):
        """Returns (effect, mutation_str) for mutation with max |effect|."""
        max_abs = 0.0
        effect = 0.0
        mutation = ""
        for protein, wt, site, m in muts:
            key = (protein, site)
            if (
                key in self.effects
                and wt in self.effects[key]
                and m in self.effects[key]
            ):
                mut_effect = self.effects[key][m] - self.effects[key][wt]
                if abs(mut_effect) > max_abs:
                    max_abs = abs(mut_effect)
                    effect = mut_effect
                    mutation = f"{protein} {wt}{site}{m}"
        return effect, mutation


def compute_cumulative_mutations(aa_muts_json, tree_newick):
    """Compute cumulative AA mutations from root for each node across all genes.

    Returns dict mapping node name to list of (protein, wt, site, mutant) tuples.
    """
    with open(aa_muts_json) as f:
        aa_data = json.load(f)

    tree = Phylo.read(tree_newick, "newick")

    branch_mutations = {}
    for node, node_data in aa_data["nodes"].items():
        muts = []
        for gene, mut_strs in node_data.get("aa_muts", {}).items():
            for mut_str in mut_strs:
                wt = mut_str[0]
                mutant = mut_str[-1]
                site = int(mut_str[1:-1])
                muts.append((gene, wt, site, mutant))
        branch_mutations[node] = muts

    parent_map = {}
    for clade in tree.find_clades():
        for child in clade.clades:
            assert child.name not in parent_map, child.name
            parent_map[child.name] = clade.name

    cumulative = {}

    def get_cumulative(node_name):
        if node_name in cumulative:
            return cumulative[node_name]
        if node_name not in parent_map:
            assert node_name == tree.root.name
            cumulative[node_name] = []
            return []
        parent_muts = get_cumulative(parent_map[node_name])
        branch_muts = branch_mutations.get(node_name, [])
        muts_dict = {}
        for protein, wt, site, mut in parent_muts:
            muts_dict[(protein, site)] = (protein, wt, site, mut)
        for protein, wt, site, mut in branch_muts:
            key = (protein, site)
            if key in muts_dict:
                orig_wt = muts_dict[key][1]
                if mut == orig_wt:
                    del muts_dict[key]
                else:
                    muts_dict[key] = (protein, orig_wt, site, mut)
            else:
                muts_dict[key] = (protein, wt, site, mut)
        cumulative[node_name] = list(muts_dict.values())
        return cumulative[node_name]

    for node_name in branch_mutations:
        get_cumulative(node_name)

    return cumulative


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree-newick", required=True)
    parser.add_argument("--aa-muts", required=True)
    parser.add_argument("--mutation-effects", required=True)
    parser.add_argument("--phenotypes", required=True, nargs="+")
    parser.add_argument("--phenotype-names", required=True, nargs="+")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if len(args.phenotypes) != len(args.phenotype_names):
        raise ValueError("--phenotypes and --phenotype-names must have same length")

    effects_df = pd.read_csv(args.mutation_effects, sep="\t")
    required_cols = {
        "protein",
        "protein_site",
        "wildtype",
        "mutant",
        "phenotype",
        "effect",
    }
    if not required_cols.issubset(effects_df.columns):
        raise ValueError(
            f"mutation_effects must have columns {required_cols}, "
            f"got {set(effects_df.columns)}"
        )

    available_phenotypes = set(effects_df["phenotype"].unique())
    for p in args.phenotypes:
        if p not in available_phenotypes:
            raise ValueError(
                f"Phenotype {p!r} not in mutation_effects. "
                f"Available: {available_phenotypes}"
            )

    with open(args.aa_muts) as f:
        aa_data = json.load(f)
    aa_muts_genes = set()
    for nd in aa_data["nodes"].values():
        aa_muts_genes.update(nd.get("aa_muts", {}).keys())
    effects_proteins = set(effects_df["protein"].unique())
    if aa_muts_genes != effects_proteins:
        raise ValueError(
            f"Proteins in mutation_effects {effects_proteins} != "
            f"genes in aa_muts {aa_muts_genes}"
        )

    cumulative = compute_cumulative_mutations(args.aa_muts, args.tree_newick)
    print(f"Computed cumulative mutations for {len(cumulative)} nodes")

    node_data = {"nodes": {}}
    for phenotype, name in zip(args.phenotypes, args.phenotype_names):
        pheno_df = effects_df[effects_df["phenotype"] == phenotype]
        scorer = MutationScorer(pheno_df)
        for node_name, mutations in cumulative.items():
            if node_name not in node_data["nodes"]:
                node_data["nodes"][node_name] = {}
            total = scorer.total_effect(mutations)
            max_effect, max_mut = scorer.max_magnitude_effect(mutations)
            node_data["nodes"][node_name][f"{name} total effect"] = total
            node_data["nodes"][node_name][f"{name} max magnitude effect"] = max_effect
            node_data["nodes"][node_name][f"{name} max magnitude mutation"] = max_mut

    print(f"Writing node data to {args.output}")
    with open(args.output, "w") as f:
        json.dump(node_data, f, indent=2)
    print(
        f"Scored {len(node_data['nodes'])} nodes for {len(args.phenotypes)} phenotypes"
    )


if __name__ == "__main__":
    main()
