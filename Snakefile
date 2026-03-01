"""Snakemake file that runs analysis."""


configfile: "config.yaml"


rule all:
    input:
        "results/ncbi_dataset/ncbi_dataset_download_date.txt",
        expand(
            [
                "results/{tree}/ha_metadata_all.tsv",
                "results/{tree}/ha_cds_all.fasta.gz",
                "results/{tree}/extract_ha_cds_stats.txt",
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
        subtype_regex=lambda w: config["trees"][w.tree]["subtype"],
        cds_length_range=lambda w: config["trees"][w.tree]["cds_length_range"],
    log:
        "results/logs/extract_ha_cds_{tree}.txt",
    conda:
        "environment.yaml"
    script:
        "scripts/extract_ha_cds.py"
