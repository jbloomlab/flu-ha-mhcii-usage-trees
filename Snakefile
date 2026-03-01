"""Snakemake file that runs analysis."""

import yaml


configfile: "config.yaml"


rule all:
    input:
        "results/ncbi_dataset/ncbi_dataset_download_date.txt",
        expand(
            [
                "results/{tree}/ha_metadata_subsampled.tsv",
                "results/{tree}/ha_cds_subsampled.fasta",
            ],
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
    log:
        "results/logs/format_ncbi_dataset.txt",
    conda:
        "environment.yaml"
    shell:
        """
        dataformat tsv virus-genome \
            --package {input.zipfile} \
            --fields accession,isolate-lineage,isolate-collection-date,\
isolate-lineage-source,host-name,host-common-name,geo-location,geo-region,\
lab-host,is-lab-host,completeness,length,segment,release-date,biosample-acc,\
sra-accs,is-vaccine-strain,purpose-of-sampling \
            > {output.genome_metadata} 2> {log}
        unzip -p {input.zipfile} ncbi_dataset/data/genomic.fna \
            | grep '^>' > {output.genomic_headers} 2>> {log}
        """


rule extract_ha_cds:
    """Extract full-length HA CDS and metadata for a subtype."""
    input:
        genome_metadata=rules.format_ncbi_dataset.output.genome_metadata,
        genomic_headers=rules.format_ncbi_dataset.output.genomic_headers,
        zipfile=rules.get_ncbi_dataset_zip.output.zipfile,
    output:
        metadata="results/{tree}/ha_metadata_all.tsv",
        fasta="results/{tree}/ha_cds_all.fasta.gz",
        stats="results/{tree}/extract_ha_cds_stats.txt",
    params:
        subtype_regex=lambda wc: config["trees"][wc.tree]["subtype"],
        cds_length_range=lambda wc: config["trees"][wc.tree]["cds_length_range"],
    log:
        "results/logs/extract_ha_cds_{tree}.txt",
    conda:
        "environment.yaml"
    script:
        "scripts/extract_ha_cds.py"


rule subsample:
    """Subsample the HAs for a tree using augur subsample."""
    input:
        sequences=rules.extract_ha_cds.output.fasta,
        metadata=rules.extract_ha_cds.output.metadata,
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
        sequences="results/{tree}/ha_cds_subsampled.fasta",
        metadata="results/{tree}/ha_metadata_subsampled.tsv",
        config="results/{tree}/subsample_config.yaml",
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
