"""Snakemake file that runs analysis."""

import copy
import os
import shlex

import yaml


configfile: "config.yaml"


# Constrain {tree} and {pruned_tree} so neither wildcard matches the other
# rule's output (rule `export` and rule `prune_tree` both write to
# `auspice/{auspice_prefix}_*.json`).
wildcard_constraints:
    tree="|".join(config["trees"]),
    pruned_tree="|".join(config["pruned_trees"]),


rule all:
    input:
        auspice_jsons=expand(
            os.path.join("auspice", config["auspice_prefix"] + "_{tree}.json"),
            tree=config["trees"],
        ),
        pruned_auspice_jsons=expand(
            os.path.join("auspice", config["auspice_prefix"] + "_{pruned_tree}.json"),
            pruned_tree=config["pruned_trees"],
        ),
        bad_dates=expand(
            "results/trees/{tree}/accessions_with_bad_dates_reset.tsv",
            tree=config["trees"],
        ),
        include_validation=expand(
            "results/trees/{tree}/include_validation.tsv",
            tree=config["trees"],
        ),


rule get_ncbi_dataset_zip:
    output:
        zipfile="results/ncbi_dataset/ncbi_dataset.zip",
        datestamp="results/ncbi_dataset/ncbi_dataset_download_date.txt",
    log:
        "results/logs/get_ncbi_dataset_zip.txt",
    conda:
        "environment.yaml"
    """Get zipped dataset of all influenza sequences."""
    params:
        taxid=config["dataset_taxid"],
    shell:
        """
        datasets download virus genome taxon {params.taxid} \
            --include genome,cds,annotation,biosample \
            --no-progressbar \
            --filename {output.zipfile} \
            &> {log}
        date +%Y-%m-%d > {output.datestamp} 2> {log}
        """


rule format_ncbi_dataset:
    """Extract genome metadata TSV and genomic FASTA headers from the zip."""
    input:
        zipfile=rules.get_ncbi_dataset_zip.output.zipfile,
    output:
        genome_metadata="results/ncbi_dataset/genome_metadata.tsv",
        genomic_headers="results/ncbi_dataset/genomic_headers.txt",
        annotation_genes="results/ncbi_dataset/annotation_genes.tsv",
    log:
        "results/logs/format_ncbi_dataset.txt",
    conda:
        "environment.yaml"
    shell:
        """
        dataformat tsv virus-genome \
            --package {input.zipfile} \
            --fields accession,isolate-lineage,isolate-collection-date,\
isolate-lineage-source,host-name,host-tax-id,geo-location,geo-region,\
lab-host,is-lab-host,completeness,length,segment,release-date,biosample-acc,\
sra-accs,is-vaccine-strain,purpose-of-sampling \
            > {output.genome_metadata} 2> {log}
        unzip -p {input.zipfile} ncbi_dataset/data/genomic.fna \
            | grep '^>' > {output.genomic_headers} 2>> {log}
        dataformat tsv virus-annotation \
            --package {input.zipfile} \
            --fields accession,gene-name \
            > {output.annotation_genes} 2>> {log}
        """


rule extract_ha_cds:
    """Extract full-length HA CDS and metadata for a subtype."""
    input:
        genome_metadata=rules.format_ncbi_dataset.output.genome_metadata,
        genomic_headers=rules.format_ncbi_dataset.output.genomic_headers,
        annotation_genes=rules.format_ncbi_dataset.output.annotation_genes,
        zipfile=rules.get_ncbi_dataset_zip.output.zipfile,
        cds_length_range=lambda wc: config["trees"][wc.tree]["cds_length_range"],
    output:
        metadata="results/trees/{tree}/metadata_all_pre_host.tsv",
        fasta="results/trees/{tree}/cds_all.fasta.gz",
        stats="results/trees/{tree}/extract_cds_stats.txt",
    log:
        "results/logs/extract_ha_cds_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        subtype_regex=lambda wc: config["trees"][wc.tree]["subtype"],
    script:
        "scripts/extract_ha_cds.py"


rule download_taxonomy:
    """Download NCBI taxonomy dump for taxonkit."""
    output:
        taxonomy_dir=directory("results/ncbi_taxonomy"),
    log:
        "results/logs/download_taxonomy.txt",
    conda:
        "environment.yaml"
    shell:
        """
        mkdir -p {output.taxonomy_dir} 2> {log}
        wget -q -O {output.taxonomy_dir}/taxdump.tar.gz \
            http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz \
            2>> {log}
        tar -xzf {output.taxonomy_dir}/taxdump.tar.gz \
            -C {output.taxonomy_dir} \
            names.dmp nodes.dmp delnodes.dmp merged.dmp \
            2>> {log}
        rm {output.taxonomy_dir}/taxdump.tar.gz 2>> {log}
        """


rule merge_manual_add:
    """Merge manual-add CDS sequences and metadata into Genbank-extracted data."""
    input:
        genbank_metadata=rules.extract_ha_cds.output.metadata,
        genbank_sequences=rules.extract_ha_cds.output.fasta,
        manual_add_metadata=lambda wc: config["trees"][wc.tree]["manual_adds"][
            "metadata"
        ],
        manual_add_sequences=lambda wc: config["trees"][wc.tree]["manual_adds"][
            "sequences"
        ],
    output:
        metadata="results/trees/{tree}/metadata_all_pre_host_merged.tsv",
        sequences="results/trees/{tree}/cds_all_merged.fasta.gz",
    log:
        "results/logs/merge_manual_add_{tree}.txt",
    conda:
        "environment.yaml"
    script:
        "scripts/merge_manual_add.py"


rule build_include_file:
    """Combine accessions_to_include.txt with manual-add accessions."""
    input:
        accessions_to_include=lambda wc: config["trees"][wc.tree][
            "accessions_to_include"
        ],
        manual_add_metadata=lambda wc: config["trees"][wc.tree]["manual_adds"][
            "metadata"
        ],
    output:
        include="results/trees/{tree}/accessions_to_include_combined.txt",
    log:
        "results/logs/build_include_file_{tree}.txt",
    conda:
        "environment.yaml"
    script:
        "scripts/build_include_file.py"


rule annotate_host_taxonomy:
    """Annotate metadata with host taxonomy (general class, order)."""
    input:
        metadata=rules.merge_manual_add.output.metadata,
        taxonomy_dir=rules.download_taxonomy.output.taxonomy_dir,
    output:
        metadata="results/trees/{tree}/metadata_all.tsv",
    log:
        "results/logs/annotate_host_taxonomy_{tree}.txt",
    conda:
        "environment.yaml"
    script:
        "scripts/annotate_host_taxonomy.py"


def _subsample_config_with_include(wc):
    """Return augur_subsample config with the combined-include path injected."""
    cfg = copy.deepcopy(config["trees"][wc.tree]["augur_subsample"])
    include_path = rules.build_include_file.output.include.format(tree=wc.tree)
    for key in ["defaults", "samples"]:
        if key not in cfg:
            continue
        entries = [cfg[key]] if key == "defaults" else cfg[key].values()
        for entry in entries:
            entry["include"] = include_path
    return cfg


rule subsample:
    """Subsample the CDSs for a tree using augur subsample."""
    input:
        sequences=rules.merge_manual_add.output.sequences,
        metadata=rules.annotate_host_taxonomy.output.metadata,
        combined_include=rules.build_include_file.output.include,
        # the list below is files listing sequences to exclude
        exclude_files=lambda wc: [
            cfg["exclude"]
            for key1 in ["defaults", "samples"]
            if key1 in config["trees"][wc.tree]["augur_subsample"]
            for cfg in (
                [config["trees"][wc.tree]["augur_subsample"][key1]]
                if key1 == "defaults"
                else config["trees"][wc.tree]["augur_subsample"][key1].values()
            )
            if "exclude" in cfg
        ],
    output:
        sequences="results/trees/{tree}/cds_subsampled.fasta",
        metadata="results/trees/{tree}/metadata_subsampled.tsv",
        config="results/trees/{tree}/subsample_config.yaml",
    log:
        "results/logs/subsample_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        build_config_yaml=lambda wc: yaml.dump(
            _subsample_config_with_include(wc), default_flow_style=False
        ),
        seed=1,
        strain_id="accession",
    shell:
        """
        cat > {output.config} <<'EOF'
{params.build_config_yaml}EOF

        augur subsample \
            --metadata-id-columns {params.strain_id} \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --config {output.config} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --seed {params.seed} \
            &> {log}
        """


rule collapse_host_order:
    """Collapse low-frequency host orders to 'other' in subsampled metadata."""
    input:
        metadata=rules.subsample.output.metadata,
    output:
        metadata="results/trees/{tree}/metadata_subsampled_collapsed.tsv",
    log:
        "results/logs/collapse_host_order_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        collapse_low_freq_host_order=config["collapse_low_freq_host_order"],
    script:
        "scripts/collapse_host_order.py"


rule align:
    """Align the sequences using augur align."""
    input:
        sequences=rules.subsample.output.sequences,
        reference_sequence=lambda wc: config["trees"][wc.tree]["reference_sequence"],
    output:
        alignment="results/trees/{tree}/alignment.fasta",
    log:
        "results/logs/align_{tree}.txt",
    conda:
        "environment.yaml"
    threads: 4
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --output {output.alignment} \
            --nthreads {threads} \
            --reference-sequence {input.reference_sequence} \
            --remove-reference \
            &> {log}
        """


rule tree:
    """Infer tree using augur tree."""
    input:
        alignment=rules.align.output.alignment,
    output:
        tree="results/trees/{tree}/tree_raw.nwk",
    log:
        "results/logs/tree_{tree}.txt",
    conda:
        "environment.yaml"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads} \
            --tree-builder-args='-seed 1 -czb' \
            &> {log}
        """


rule refine:
    """Refine the tree using augur refine."""
    input:
        tree=rules.tree.output.tree,
        alignment=rules.align.output.alignment,
        metadata=rules.collapse_host_order.output.metadata,
        keep_ids=rules.build_include_file.output.include,
    output:
        tree="results/trees/{tree}/tree.nwk",
        node_data="results/trees/{tree}/branch_lengths.json",
        refine_output="results/trees/{tree}/refine_output.txt",
    log:
        "results/logs/refine_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        strain_id="accession",
        addtl_flags=lambda wc: " ".join(
            f"--{k} {v}" for k, v in config["trees"][wc.tree]["augur_refine"].items()
        ),
    shell:
        """
        set -euo pipefail
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --use-fft \
            --keep-ids {input.keep_ids} \
            {params.addtl_flags} \
            2>&1 | tee {output.refine_output} > {log}
        """


rule parse_refine_outliers:
    """Parse augur refine output to extract strains with dates reset as outliers."""
    input:
        refine_output=rules.refine.output.refine_output,
    output:
        tsv="results/trees/{tree}/accessions_with_bad_dates_reset.tsv",
    log:
        "results/logs/parse_refine_outliers_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        bad_dates_action=config["bad_dates_in_keep_accessions_action"],
    script:
        "scripts/parse_refine_outliers.py"


rule ancestral:
    """Infer ancestral nucleotide sequences with augur ancestral."""
    input:
        tree=rules.refine.output.tree,
        alignment=rules.align.output.alignment,
        root_sequence=lambda wc: config["trees"][wc.tree]["reference_sequence"],
    output:
        node_data="results/trees/{tree}/nt_muts.json",
    log:
        "results/logs/ancestral_{tree}.txt",
    conda:
        "environment.yaml"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --root-sequence {input.root_sequence} \
            --output-node-data {output.node_data} \
            --seed 1 \
            &> {log}
        """


rule translate:
    """Translate nucleotide mutations to amino acid mutations."""
    input:
        tree=rules.refine.output.tree,
        node_data=rules.ancestral.output.node_data,
        annotation=lambda wc: config["trees"][wc.tree]["annotation"],
    output:
        node_data="results/trees/{tree}/aa_muts.json",
    log:
        "results/logs/translate_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        genes=lambda wc: " ".join(config["trees"][wc.tree]["genes"]),
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.annotation} \
            --genes {params.genes} \
            --output-node-data {output.node_data} \
            &> {log}
        """


rule score_mutation_effects:
    """Score tree nodes based on per-mutation phenotype effects."""
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        mutation_effects=lambda wc: config["trees"][wc.tree]["mutation_effects"],
    output:
        node_data="results/trees/{tree}/mutation_effects_scores.json",
    log:
        "results/logs/score_mutation_effects_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        phenotypes=lambda wc: " ".join(
            shlex.quote(p) for p in config["trees"][wc.tree]["phenotypes"]
        ),
        phenotype_names=lambda wc: " ".join(
            shlex.quote(config["trees"][wc.tree]["phenotypes"][p])
            for p in config["trees"][wc.tree]["phenotypes"]
        ),
    shell:
        """
        python scripts/score_mutation_effects.py \
            --tree-newick {input.tree} \
            --aa-muts {input.aa_muts} \
            --mutation-effects {input.mutation_effects} \
            --phenotypes {params.phenotypes} \
            --phenotype-names {params.phenotype_names} \
            --output {output.node_data} \
            &> {log}
        """


rule generate_phenotype_auspice_config:
    """Generate auspice config with color scales for phenotype scores."""
    input:
        node_data=rules.score_mutation_effects.output.node_data,
    output:
        auspice_config="results/trees/{tree}/auspice_config_phenotypes.json",
    log:
        "results/logs/generate_phenotype_auspice_config_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        phenotype_names=lambda wc: " ".join(
            shlex.quote(config["trees"][wc.tree]["phenotypes"][p])
            for p in config["trees"][wc.tree]["phenotypes"]
        ),
        continuous_scale=lambda wc: " ".join(
            shlex.quote(c) for c in config["trees"][wc.tree]["phenotype_color_scale"]
        ),
    shell:
        """
        python scripts/generate_phenotype_auspice_config.py \
            --phenotypes {params.phenotype_names} \
            --continuous-scale {params.continuous_scale} \
            --node-data {input.node_data} \
            --output {output.auspice_config} \
            &> {log}
        """


rule annotate_strains:
    """Annotate tips from the per-strain annotations TSV for a tree."""
    input:
        tree=rules.refine.output.tree,
        strain_annotations=lambda wc: config["trees"][wc.tree]["strain_annotations"],
    output:
        node_data="results/trees/{tree}/strain_annotations_node_data.json",
        auspice_config="results/trees/{tree}/auspice_config_strain_annotations.json",
    log:
        "results/logs/annotate_strains_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        missing_strain_annotations_action=config["missing_strain_annotations_action"],
        phenotype_color_scale=lambda wc: config["trees"][wc.tree][
            "phenotype_color_scale"
        ],
    script:
        "scripts/annotate_strains.py"


rule format_description:
    """Format description markdown, replacing {tree} with the tree name."""
    input:
        description=lambda wc: config["trees"][wc.tree]["description"],
    output:
        description="results/trees/{tree}/description.md",
    log:
        "results/logs/format_description_{tree}.txt",
    conda:
        "environment.yaml"
    shell:
        """
        sed 's/{{tree}}/{wildcards.tree}/g' {input.description} \
            > {output.description} 2> {log}
        """


rule export:
    """Export auspice json."""
    input:
        tree=rules.refine.output.tree,
        branch_lengths=rules.refine.output.node_data,
        nt_muts=rules.ancestral.output.node_data,
        aa_muts=rules.translate.output.node_data,
        mutation_effects_scores=rules.score_mutation_effects.output.node_data,
        strain_annotations_node_data=rules.annotate_strains.output.node_data,
        metadata=rules.collapse_host_order.output.metadata,
        auspice_config=lambda wc: config["trees"][wc.tree]["auspice_config"],
        phenotype_auspice_config=rules.generate_phenotype_auspice_config.output.auspice_config,
        strain_annotations_auspice_config=rules.annotate_strains.output.auspice_config,
        description=rules.format_description.output.description,
    output:
        auspice_json=os.path.join("auspice", config["auspice_prefix"] + "_{tree}.json"),
    log:
        "results/logs/export_{tree}.txt",
    conda:
        "environment.yaml"
    params:
        strain_id="accession",
        title=lambda wc: shlex.quote(config["trees"][wc.tree]["title"]),
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.mutation_effects_scores} {input.strain_annotations_node_data} \
            --include-root-sequence-inline \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --output {output.auspice_json} \
            --title {params.title} \
            --auspice-config {input.auspice_config} {input.phenotype_auspice_config} {input.strain_annotations_auspice_config} \
            --description {input.description} \
            &> {log}
        """


rule prune_tree:
    """Prune a source Auspice JSON to tips matching a strain-annotations filter."""
    input:
        source_json=lambda wc: os.path.join(
            "auspice",
            config["auspice_prefix"]
            + "_"
            + config["pruned_trees"][wc.pruned_tree]["source_tree"]
            + ".json",
        ),
        strain_annotations=lambda wc: config["trees"][
            config["pruned_trees"][wc.pruned_tree]["source_tree"]
        ]["strain_annotations"],
    output:
        auspice_json=os.path.join(
            "auspice", config["auspice_prefix"] + "_{pruned_tree}.json"
        ),
    params:
        keep_where=lambda wc: config["pruned_trees"][wc.pruned_tree]["keep_where"],
        title=lambda wc: config["pruned_trees"][wc.pruned_tree].get("title"),
    log:
        "results/logs/prune_tree_{pruned_tree}.txt",
    conda:
        "environment.yaml"
    script:
        "scripts/prune_auspice_json.py"


rule validate_includes:
    """Check every accession specified to include is a tip in the final Auspice JSON."""
    input:
        auspice_json=rules.export.output.auspice_json,
        accessions_to_include=lambda wc: config["trees"][wc.tree][
            "accessions_to_include"
        ],
        manual_add_metadata=lambda wc: config["trees"][wc.tree]["manual_adds"][
            "metadata"
        ],
    output:
        tsv="results/trees/{tree}/include_validation.tsv",
    log:
        "results/logs/validate_includes_{tree}.txt",
    conda:
        "environment.yaml"
    script:
        "scripts/validate_includes.py"
