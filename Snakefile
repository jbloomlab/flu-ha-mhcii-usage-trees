"""Snakemake file that runs analysis."""


configfile: "config.yaml"


rule all:
    input:
        f"results/ncbi_dataset/{config['datasets_taxid']}_download_date.txt",


rule get_ncbi_dataset_zip:
    """Get zipped dataset of all influenza sequences."""
    params:
        taxid=config["datasets_taxid"],
    output:
        zipfile=f"results/ncbi_dataset/{config['datasets_taxid']}.zip",
        datestamp=f"results/ncbi_dataset/{config['datasets_taxid']}_download_date.txt",
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
