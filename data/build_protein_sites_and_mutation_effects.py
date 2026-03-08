"""Build ``protein_sites.tsv`` and ``mutation_effects.tsv`` for each tree.

See the README for details on what this script does.
"""

import itertools
import os
import subprocess
import tempfile

import pandas as pd
import yaml
from Bio import AlignIO, SeqIO

TREES_DIR = "trees"
CONFIG_YAML = "build_protein_sites_and_mutation_effects_config.yaml"
ALIGNMENT_OUT = "build_protein_sites_and_mutation_effects_alignment.fa"
IDENTITIES_OUT = "build_protein_sites_and_mutation_effects_protein_identities.txt"
EXPECTED_PROTEINS = ["SigPep", "HA1", "HA2"]


def read_reference_proteins(trees_dir):
    """Read and translate reference sequences for all trees."""
    proteins = {}
    for tree in sorted(os.listdir(trees_dir)):
        tree_dir = os.path.join(trees_dir, tree)
        if not os.path.isdir(tree_dir):
            continue
        fa_path = os.path.join(tree_dir, "reference_sequence.fa")
        if not os.path.isfile(fa_path):
            raise ValueError(f"No reference_sequence.fa in {tree_dir}")
        records = list(SeqIO.parse(fa_path, "fasta"))
        if len(records) != 1:
            raise ValueError(f"Expected 1 record in {fa_path}, got {len(records)}")
        protein = str(records[0].seq.translate(cds=True))
        proteins[tree] = protein
        print(f"  {tree}: {len(protein)} amino acids")
    if not proteins:
        raise ValueError(f"No tree subdirectories found in {trees_dir}")
    return proteins


def read_config(config_path):
    """Read and validate the YAML configuration."""
    with open(config_path) as f:
        config = yaml.safe_load(f)
    if "mutation_effect_data" not in config:
        raise ValueError(f"Config {config_path} missing 'mutation_effect_data' key")
    all_display_names = []
    for name, entry in config["mutation_effect_data"].items():
        for key in ["csv", "phenotypes"]:
            if key not in entry:
                raise ValueError(f"mutation_effect_data[{name}] missing '{key}'")
        if not os.path.isfile(entry["csv"]):
            raise ValueError(f"CSV file not found: {entry['csv']}")
        all_display_names.extend(entry["phenotypes"].values())
    if len(all_display_names) != len(set(all_display_names)):
        raise ValueError(f"Duplicate phenotype display names: {all_display_names}")
    return config


def reconstruct_experiment_proteins(config, ref_names):
    """Reconstruct protein sequences from mutation effect CSVs."""
    proteins = {}
    dataframes = {}
    for name, entry in config["mutation_effect_data"].items():
        if name in ref_names:
            raise ValueError(
                f"mutation_effect_data name '{name}' collides with a tree name"
            )
        df = pd.read_csv(entry["csv"])
        for col in ["sequential_site", "wildtype", "mutant"]:
            if col not in df.columns:
                raise ValueError(f"CSV {entry['csv']} missing column '{col}'")
        for pheno_col in entry["phenotypes"]:
            if pheno_col not in df.columns:
                raise ValueError(
                    f"CSV {entry['csv']} missing phenotype column '{pheno_col}'"
                )
        # get unique (sequential_site, wildtype) pairs
        site_wt = df[["sequential_site", "wildtype"]].drop_duplicates()
        site_wt = site_wt.dropna(subset=["sequential_site"])
        site_wt["sequential_site"] = site_wt["sequential_site"].astype(int)
        # validate consistency: each site has exactly one wildtype
        wt_per_site = site_wt.groupby("sequential_site")["wildtype"].nunique()
        inconsistent = wt_per_site[wt_per_site > 1]
        if len(inconsistent) > 0:
            raise ValueError(
                f"Inconsistent wildtype at sites in {entry['csv']}: "
                f"{inconsistent.index.tolist()}"
            )
        site_to_wt = dict(zip(site_wt["sequential_site"], site_wt["wildtype"]))
        max_site = max(site_to_wt)
        protein_chars = []
        for i in range(1, max_site + 1):
            if i in site_to_wt:
                aa = site_to_wt[i]
                if len(aa) != 1:
                    raise ValueError(
                        f"Wildtype '{aa}' at site {i} in {entry['csv']} "
                        f"is not a single character"
                    )
                protein_chars.append(aa)
            else:
                protein_chars.append("X")
        protein = "".join(protein_chars)
        proteins[name] = protein
        dataframes[name] = df
        n_x = protein.count("X")
        print(f"  {name}: {len(protein)} amino acids ({n_x} gap-filled with X)")
    return proteins, dataframes


def run_mafft(all_proteins, output_path):
    """Run mafft on all proteins and return the alignment."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        for name, seq in all_proteins.items():
            tmp.write(f">{name}\n{seq}\n")
        tmp_path = tmp.name
    try:
        result = subprocess.run(
            ["mafft", "--auto", tmp_path],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"mafft failed:\n{e.stderr}") from e
    finally:
        os.unlink(tmp_path)
    with open(output_path, "w") as f:
        f.write(result.stdout)
    alignment = AlignIO.read(output_path, "fasta")
    # return as dict of name -> aligned sequence string
    aln_dict = {rec.id: str(rec.seq) for rec in alignment}
    for name in all_proteins:
        if name not in aln_dict:
            raise ValueError(f"Protein '{name}' not found in mafft output")
    return aln_dict


def parse_gff(gff_path, protein_length):
    """Parse a GFF file and return amino acid ranges for each protein.

    The GFF covers the full CDS including the stop codon, but the translated
    protein excludes it. ``protein_length`` is the length of the translated
    protein (without stop). The last protein's end is adjusted accordingly.

    Returns dict mapping protein name to (aa_start, aa_end) with 1-based
    inclusive coordinates.
    """
    proteins = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                raise ValueError(f"Unexpected GFF line in {gff_path}: {line}")
            nuc_start = int(parts[3])
            nuc_end = int(parts[4])
            attrs = parts[8]
            gene_name = None
            for attr in attrs.split(";"):
                if attr.startswith("gene_name="):
                    gene_name = attr.split("=", 1)[1]
            if gene_name is None:
                raise ValueError(f"No gene_name in GFF line: {line}")
            if gene_name not in EXPECTED_PROTEINS:
                raise ValueError(
                    f"Unexpected protein '{gene_name}' in {gff_path}, "
                    f"expected {EXPECTED_PROTEINS}"
                )
            # convert nucleotide 1-based coords to amino acid 1-based coords
            if nuc_start % 3 != 1:
                raise ValueError(
                    f"Protein {gene_name} nuc_start {nuc_start} not at codon "
                    f"boundary in {gff_path}"
                )
            if nuc_end % 3 != 0:
                raise ValueError(
                    f"Protein {gene_name} nuc_end {nuc_end} not at codon "
                    f"boundary in {gff_path}"
                )
            aa_start = (nuc_start - 1) // 3 + 1
            aa_end = nuc_end // 3
            proteins[gene_name] = (aa_start, aa_end)
    for prot in EXPECTED_PROTEINS:
        if prot not in proteins:
            raise ValueError(f"Missing protein '{prot}' in {gff_path}")
    # adjust last protein's end to exclude stop codon
    last_prot = EXPECTED_PROTEINS[-1]
    aa_start, aa_end = proteins[last_prot]
    if aa_end != protein_length + 1:
        raise ValueError(
            f"GFF {gff_path} last protein {last_prot} ends at aa {aa_end} "
            f"but expected {protein_length + 1} (protein length {protein_length} "
            f"+ 1 for stop codon)"
        )
    proteins[last_prot] = (aa_start, aa_end - 1)
    return proteins


def build_site_to_protein_map(gff_proteins):
    """Build mapping from sequential site to (protein, protein_site)."""
    site_map = {}
    for prot in EXPECTED_PROTEINS:
        aa_start, aa_end = gff_proteins[prot]
        for seq_site in range(aa_start, aa_end + 1):
            protein_site = seq_site - aa_start + 1
            site_map[seq_site] = (prot, protein_site)
    return site_map


def compute_h3_sites(aln_dict, h3_gff_proteins):
    """Compute H3 site for each alignment column.

    Returns a list (one per alignment column) of H3 site strings.
    H3 ectodomain site 1 = first residue of HA1. Signal peptide sites are
    negative. Insertions relative to H3 get letter suffixes (e.g., 125a).
    """
    h3_seq = aln_dict["H3"]
    aln_len = len(h3_seq)
    ha1_start_aa = h3_gff_proteins["HA1"][0]  # sequential site where HA1 starts

    h3_sites = []
    h3_pos = 0  # tracks sequential position in H3 (1-based, after increment)
    last_h3_ecto_site = None
    insertion_counter = 0

    for col in range(aln_len):
        if h3_seq[col] != "-":
            h3_pos += 1
            # HA1 start -> ecto site 1, no site 0 (goes ..., -2, -1, 1, 2, ...)
            ecto_site = h3_pos - ha1_start_aa + 1
            if ecto_site <= 0:
                ecto_site -= 1
            last_h3_ecto_site = ecto_site
            insertion_counter = 0
            h3_sites.append(str(ecto_site))
        else:
            # insertion relative to H3
            if last_h3_ecto_site is None:
                # before any H3 residue; shouldn't normally happen with a
                # good alignment but handle gracefully
                h3_sites.append(None)
            else:
                suffix = chr(ord("a") + insertion_counter)
                h3_sites.append(f"{last_h3_ecto_site}{suffix}")
                insertion_counter += 1
    return h3_sites


def build_protein_sites(tree, aln_dict, h3_sites, trees_dir, ref_proteins):
    """Build protein_sites.tsv for a single tree."""
    gff_path = os.path.join(trees_dir, tree, "annotation.gff")
    protein_length = len(ref_proteins[tree])
    gff_proteins = parse_gff(gff_path, protein_length)
    site_map = build_site_to_protein_map(gff_proteins)
    ha1_start_aa = gff_proteins["HA1"][0]

    tree_seq = aln_dict[tree]
    aln_len = len(tree_seq)
    total_aa = sum(1 for c in tree_seq if c != "-")

    # validate site_map covers all sequential sites
    if max(site_map) != total_aa:
        raise ValueError(
            f"GFF for {tree} covers up to site {max(site_map)} but protein "
            f"has {total_aa} amino acids"
        )

    rows = []
    seq_pos = 0
    for col in range(aln_len):
        if tree_seq[col] == "-":
            continue
        seq_pos += 1
        aa = tree_seq[col]
        protein, protein_site = site_map[seq_pos]
        ecto_site = seq_pos - ha1_start_aa + 1
        if ecto_site <= 0:
            ecto_site -= 1
        h3_site = h3_sites[col]
        rows.append(
            {
                "reference_sequential_site": seq_pos,
                "reference_aa": aa,
                "protein": protein,
                "protein_site": protein_site,
                "ectodomain_site": ecto_site,
                "H3_site": h3_site,
            }
        )

    df = pd.DataFrame(rows).sort_values("reference_sequential_site")
    out_path = os.path.join(trees_dir, tree, "protein_sites.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  {tree}: {len(df)} sites -> {out_path}")
    return df


def build_mutation_effects(tree, tree_protein_sites, aln_dict, config, dataframes):
    """Build mutation_effects.tsv for a single tree."""
    tree_seq = aln_dict[tree]
    aln_len = len(tree_seq)

    # build mapping: alignment column -> reference sequential site (for tree)
    col_to_tree_site = {}
    seq_pos = 0
    for col in range(aln_len):
        if tree_seq[col] != "-":
            seq_pos += 1
            col_to_tree_site[col] = seq_pos

    # build lookup from reference_sequential_site to protein_sites row
    site_lookup = tree_protein_sites.set_index("reference_sequential_site")

    all_rows = []
    for data_name, entry in config["mutation_effect_data"].items():
        data_seq = aln_dict[data_name]
        df = dataframes[data_name]

        # build mapping: experiment sequential site -> alignment column
        exp_site_to_col = {}
        exp_pos = 0
        for col in range(aln_len):
            if data_seq[col] != "-":
                exp_pos += 1
                exp_site_to_col[exp_pos] = col

        # build mapping: experiment sequential site -> tree sequential site
        exp_to_tree = {}
        for exp_site, aln_col in exp_site_to_col.items():
            if aln_col in col_to_tree_site:
                exp_to_tree[exp_site] = col_to_tree_site[aln_col]

        phenotype_map = entry["phenotypes"]

        for _, row in df.iterrows():
            if pd.isna(row["sequential_site"]):
                continue
            exp_site = int(row["sequential_site"])
            wt = row["wildtype"]
            mut = row["mutant"]
            if wt == mut:
                continue
            if exp_site not in exp_to_tree:
                continue  # insertion relative to this tree's reference
            tree_site = exp_to_tree[exp_site]
            site_info = site_lookup.loc[tree_site]

            for csv_col, display_name in phenotype_map.items():
                effect = row[csv_col]
                if pd.isna(effect):
                    continue
                all_rows.append(
                    {
                        "protein": site_info["protein"],
                        "protein_site": site_info["protein_site"],
                        "reference_aa": site_info["reference_aa"],
                        "H3_site": site_info["H3_site"],
                        "wildtype": wt,
                        "mutant": mut,
                        "phenotype": display_name,
                        "effect": effect,
                    }
                )

    result_df = pd.DataFrame(all_rows)
    protein_order = pd.CategoricalDtype(EXPECTED_PROTEINS, ordered=True)
    result_df["protein"] = result_df["protein"].astype(protein_order)
    result_df = result_df.sort_values(
        ["protein", "protein_site", "mutant", "phenotype"]
    )
    out_path = os.path.join(TREES_DIR, tree, "mutation_effects.tsv")
    result_df.to_csv(out_path, sep="\t", index=False)
    print(f"  {tree}: {len(result_df)} mutation effect rows -> {out_path}")
    return result_df


def write_pairwise_identities(aln_dict, output_path):
    """Write pairwise identity at aligned (both non-gap) sites for all pairs."""
    names = list(aln_dict)
    aln_len = len(next(iter(aln_dict.values())))
    # compute all pairwise identities
    identities = {}
    for name1, name2 in itertools.combinations(names, 2):
        seq1 = aln_dict[name1]
        seq2 = aln_dict[name2]
        aligned = 0
        identical = 0
        for i in range(aln_len):
            if seq1[i] != "-" and seq2[i] != "-":
                aligned += 1
                if seq1[i] == seq2[i]:
                    identical += 1
        identities[(name1, name2)] = (identical, aligned)
        identities[(name2, name1)] = (identical, aligned)
    # write grouped by each protein
    lines = []
    for i, name1 in enumerate(names):
        if i > 0:
            lines.append("")
        lines.append(f"Identities to {name1}")
        for name2 in names:
            if name2 == name1:
                continue
            identical, aligned = identities[(name1, name2)]
            pct = identical / aligned * 100 if aligned > 0 else 0.0
            lines.append(f"  {name2}: {identical}/{aligned} ({pct:.1f}%)")
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    print("Reading reference sequences...")
    ref_proteins = read_reference_proteins(TREES_DIR)

    print("Reading configuration...")
    config = read_config(CONFIG_YAML)

    print("Reconstructing experiment proteins...")
    exp_proteins, dataframes = reconstruct_experiment_proteins(
        config, set(ref_proteins)
    )

    # combine all proteins preserving order: references first, then experiments
    all_proteins = {**ref_proteins, **exp_proteins}
    if "H3" not in ref_proteins:
        raise ValueError("H3 tree required for H3 numbering but not found")

    print(f"Running mafft on {len(all_proteins)} proteins...")
    aln_dict = run_mafft(all_proteins, ALIGNMENT_OUT)
    print(f"  Alignment saved to {ALIGNMENT_OUT}")

    print("Parsing H3 annotation for H3 numbering...")
    h3_gff = parse_gff(
        os.path.join(TREES_DIR, "H3", "annotation.gff"), len(ref_proteins["H3"])
    )
    h3_sites = compute_h3_sites(aln_dict, h3_gff)

    print("Building protein_sites.tsv for each tree...")
    protein_sites_by_tree = {}
    for tree in sorted(ref_proteins):
        protein_sites_by_tree[tree] = build_protein_sites(
            tree, aln_dict, h3_sites, TREES_DIR, ref_proteins
        )

    print("Building mutation_effects.tsv for each tree...")
    for tree in sorted(ref_proteins):
        build_mutation_effects(
            tree, protein_sites_by_tree[tree], aln_dict, config, dataframes
        )

    print("Computing pairwise protein identities...")
    write_pairwise_identities(aln_dict, IDENTITIES_OUT)
    print(f"  Saved to {IDENTITIES_OUT}")

    print("Done.")


if __name__ == "__main__":
    main()
