"""Get mutation effects from tufted duck MHCII binding data.

Usage:
    python get_mutation_effects.py protein_sites.tsv tufted_duck_MHCII_binding.csv mutation_effects.tsv
"""

import sys

import pandas as pd

protein_sites_tsv, mhcii_binding_csv, output_tsv = sys.argv[1:]

protein_sites = pd.read_csv(protein_sites_tsv, sep="\t")

mhcii_binding = pd.read_csv(mhcii_binding_csv)[
    [
        "sequential_site",
        "wildtype",
        "mutant",
        "tufted duck MHCII binding escape",
        "entry difference between noSA tufted duck MHCII and SA23 cells",
    ]
]

mutation_effects = mhcii_binding.melt(
    id_vars=["sequential_site", "wildtype", "mutant"],
    var_name="phenotype",
    value_name="effect",
)

mutation_effects = mutation_effects.merge(
    protein_sites,
    on="sequential_site",
    how="left",
    validate="many_to_one",
)

missing = set(mutation_effects["sequential_site"]) - set(
    protein_sites["sequential_site"]
)
if missing:
    raise ValueError(f"sequential_site values not in protein_sites: {missing}")

mutation_effects = mutation_effects.sort_values(
    ["sequential_site", "mutant", "phenotype"]
)

mutation_effects = mutation_effects[mutation_effects["effect"].notnull()]

col_order = [
    c
    for c in mutation_effects.columns
    if c not in {"wildtype", "mutant", "phenotype", "effect"}
]
mutation_effects = mutation_effects[
    col_order + ["wildtype", "mutant", "phenotype", "effect"]
]

mutation_effects.to_csv(output_tsv, sep="\t", index=False)

print(f"Wrote {len(mutation_effects)} rows to {output_tsv}")
