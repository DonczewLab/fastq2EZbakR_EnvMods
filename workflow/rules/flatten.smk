# Flatten annotation with DEXseq script
rule flatten:
    input:
        gtf=config["annotation"],
    output:
        flatgtf="results/raw_flattened/flat_genome.gtf",
    log:
        "logs/get-genome.log",
    envmodules:
        config["modules"]["htseq"],
    threads: 1
    script:
        "../scripts/flatten/dexseq_prepare_annotation.py"


# Add exon ID to flattened annotation
rule add_exon:
    input:
        "results/raw_flattened/flat_genome.gtf",
    output:
        config["flat_annotation"],
    log:
        "logs/add_exon.log",
    params:
        shellscript=workflow.source_path("../scripts/flatten/exon_ID.sh"),
    envmodules:
        config["modules"]["python"],
    threads: 1
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {input} {output}
        """
