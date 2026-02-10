# Create index for RSEM
rule RSEM_index:
    input:
        reference_genome=config["genome"],
        gtf=AandQ_ANNOTATION,
    output:
        seq="rsem_index/reference.seq",
        grp="rsem_index/reference.grp",
        ti="rsem_index/reference.ti",
        tfa="rsem_index/reference.transcripts.fa",
        idxfa="rsem_index/reference.idx.fa",
        n2g="rsem_index/reference.n2g.idx.fa",
    params:
        extra="--gtf {} {}".format(str(AandQ_ANNOTATION), str(config["rsem_index_params"])),
        prefix="rsem_index/reference",
    log:
        "logs/RSEM_index/prepare-reference.log",
    threads: 20
    envmodules:
        config["modules"]["rsem"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p rsem_index logs/RSEM_index
        rsem-prepare-reference \
          --num-threads {threads} \
          {params.extra} \
          {input.reference_genome} \
          {params.prefix} \
          1> {log} 2>&1
        """


# Run RSEM to quantify transcript abundances
rule RSEM:
    input:
        bam="results/align/{sample}-Aligned.toTranscriptome.out.bam",
        reference=multiext(
            "rsem_index/reference",
            ".grp",
            ".ti",
            ".transcripts.fa",
            ".seq",
            ".idx.fa",
            ".n2g.idx.fa",
        ),
    output:
        genes_results="results/rsem/{sample}.genes.results",
        isoforms_results="results/rsem/{sample}.isoforms.results",
        bam="results/rsem/{sample}.transcript.bam",
    params:
        # keep these so your script can use them if it references params/config
        paired_end=config["PE"],
        extra=config["rsem_quant_params"],
    log:
        "logs/RSEM/{sample}.log",
    threads: 20
    envmodules:
        config["modules"]["rsem"],
        config["modules"]["samtools"],
        config["modules"]["python"],
    script:
        "../scripts/rsem/rsem-calc.py"            

