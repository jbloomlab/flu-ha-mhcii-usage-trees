"""Snakemake file that runs analysis."""

import os
import shlex

import yaml


configfile: "config.yaml"


rule all:
    input:
        auspice_jsons=expand(
            os.path.join("auspice", config["auspice_prefix"] + "_{tree}.json"),
            tree=config["trees"],
        ),
        bad_dates=expand(
            "results/trees/{tree}/accessions_with_bad_dates_reset.tsv",
            tree=config["trees"],
        ),


rule get_ncbi_dataset_zip:
    """Get zipped dataset of all influenza sequences."""
    params:
        taxid=config["dataset_taxid"],
    output:
        zipfile="results/ncbi_dataset/ncbi_dataset.zip",
        datestamp="results/ncbi_dataset/ncbi_dataset_download_date.txt",
    log:
        "results/logs/get_ncbi_dataset_zip.txt",
    conda:
        "environment.yaml"
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
    params:
        subtype_regex=lambda wc: config["trees"][wc.tree]["subtype"],
    log:
        "results/logs/extract_ha_cds_{tree}.txt",
    conda:
        "environment.yaml"
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


rule annotate_host_taxonomy:
    """Annotate metadata with host taxonomy (general class, order)."""
    input:
        metadata=rules.extract_ha_cds.output.metadata,
        taxonomy_dir=rules.download_taxonomy.output.taxonomy_dir,
    output:
        metadata="results/trees/{tree}/metadata_all.tsv",
    log:
        "results/logs/annotate_host_taxonomy_{tree}.txt",
    conda:
        "environment.yaml"
    script:
        "scripts/annotate_host_taxonomy.py"


rule subsample:
    """Subsample the CDSs for a tree using augur subsample."""
    input:
        sequences=rules.extract_ha_cds.output.fasta,
        metadata=rules.annotate_host_taxonomy.output.metadata,
        # the list below is files listing sequences to include / exclude
        include_exclude_files=lambda wc: [
            cfg[key2]
            for key1 in ["defaults", "samples"]
            if key1 in config["trees"][wc.tree]["augur_subsample"]
            for cfg in (
                [config["trees"][wc.tree]["augur_subsample"][key1]]
                if key1 == "defaults"
                else config["trees"][wc.tree]["augur_subsample"][key1].values()
            )
            for key2 in ["exclude", "include"]
            if key2 in cfg
        ],
    output:
        sequences="results/trees/{tree}/cds_subsampled.fasta",
        metadata="results/trees/{tree}/metadata_subsampled.tsv",
        config="results/trees/{tree}/subsample_config.yaml",
    params:
        build_config_yaml=lambda wc: yaml.dump(
            config["trees"][wc.tree]["augur_subsample"], default_flow_style=False
        ),
        seed=1,
        strain_id="accession",
    conda:
        "environment.yaml"
    log:
        "results/logs/subsample_{tree}.txt",
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
    params:
        collapse_low_freq_host_order=config["collapse_low_freq_host_order"],
    log:
        "results/logs/collapse_host_order_{tree}.txt",
    conda:
        "environment.yaml"
    script:
        "scripts/collapse_host_order.py"


rule align:
    """Align the sequences using augur align."""
    input:
        sequences=rules.subsample.output.sequences,
        reference_sequence=lambda wc: config["trees"][wc.tree]["reference_sequence"],
    output:
        alignment="results/trees/{tree}/alignment.fasta",
    threads: 4
    conda:
        "environment.yaml"
    log:
        "results/logs/align_{tree}.txt",
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
    threads: 4
    conda:
        "environment.yaml"
    log:
        "results/logs/tree_{tree}.txt",
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
        lambda wc: (
            [config["trees"][wc.tree]["augur_refine"]["keep-ids"]]
            if "keep-ids" in config["trees"][wc.tree]["augur_refine"]
            else []
        ),
        tree=rules.tree.output.tree,
        alignment=rules.align.output.alignment,
        metadata=rules.collapse_host_order.output.metadata,
    output:
        tree="results/trees/{tree}/tree.nwk",
        node_data="results/trees/{tree}/branch_lengths.json",
        refine_output="results/trees/{tree}/refine_output.txt",
    params:
        strain_id="accession",
        addtl_flags=lambda wc: " ".join(
            f"--{k} {v}" for k, v in config["trees"][wc.tree]["augur_refine"].items()
        ),
    conda:
        "environment.yaml"
    log:
        "results/logs/refine_{tree}.txt",
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
    conda:
        "environment.yaml"
    log:
        "results/logs/ancestral_{tree}.txt",
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
    params:
        genes=lambda wc: " ".join(config["trees"][wc.tree]["genes"]),
    conda:
        "environment.yaml"
    log:
        "results/logs/translate_{tree}.txt",
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
    params:
        phenotypes=lambda wc: " ".join(
            shlex.quote(p) for p in config["trees"][wc.tree]["phenotypes"]
        ),
        phenotype_names=lambda wc: " ".join(
            shlex.quote(config["trees"][wc.tree]["phenotypes"][p])
            for p in config["trees"][wc.tree]["phenotypes"]
        ),
    conda:
        "environment.yaml"
    log:
        "results/logs/score_mutation_effects_{tree}.txt",
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
    params:
        phenotype_names=lambda wc: " ".join(
            shlex.quote(config["trees"][wc.tree]["phenotypes"][p])
            for p in config["trees"][wc.tree]["phenotypes"]
        ),
        continuous_scale=lambda wc: " ".join(
            shlex.quote(c) for c in config["trees"][wc.tree]["phenotype_color_scale"]
        ),
    conda:
        "environment.yaml"
    log:
        "results/logs/generate_phenotype_auspice_config_{tree}.txt",
    shell:
        """
        python scripts/generate_phenotype_auspice_config.py \
            --phenotypes {params.phenotype_names} \
            --continuous-scale {params.continuous_scale} \
            --node-data {input.node_data} \
            --output {output.auspice_config} \
            &> {log}
        """


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
        metadata=rules.collapse_host_order.output.metadata,
        auspice_config=lambda wc: config["trees"][wc.tree]["auspice_config"],
        phenotype_auspice_config=rules.generate_phenotype_auspice_config.output.auspice_config,
        description=rules.format_description.output.description,
    output:
        auspice_json=os.path.join("auspice", config["auspice_prefix"] + "_{tree}.json"),
    params:
        strain_id="accession",
        title=lambda wc: shlex.quote(config["trees"][wc.tree]["title"]),
    conda:
        "environment.yaml"
    log:
        "results/logs/export_{tree}.txt",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.mutation_effects_scores} \
            --include-root-sequence-inline \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --output {output.auspice_json} \
            --title {params.title} \
            --auspice-config {input.auspice_config} {input.phenotype_auspice_config} \
            --description {input.description} \
            &> {log}
        """
