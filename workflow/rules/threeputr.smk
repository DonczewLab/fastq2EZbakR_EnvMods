"""
Rules to build a custom 3'UTR annotation from 3'-end data
"""

rule get_informative_read:
    input:
        "results/sf_reads/{sample}.s.bam",
    output:
        "results/informative_read/{sample}_informative.bam",
    log:
        "logs/get_informative_read/{sample}.log",
    envmodules:
        config["modules"]["samtools"],
    threads: 8
    params:
        informative_flag=lambda wildcards: (
            "64" if config["strandedness"] == "reverse" else "128"
        ),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})" logs/get_informative_read
        samtools view -@ {threads} -hb -f {params.informative_flag} {input} -o {output} 1> {log} 2>&1
        """


rule bam_to_3pend_bg:
    input:
        fetch_informative_read,
    output:
        "results/bam2bg/{sample}_informative_{strand}.bg",
    log:
        "logs/bam_to_3pend_bg/{sample}_{strand}.log",
    envmodules:
        config["modules"]["samtools"],
        config["modules"]["bedtools"],
    params:
        strandedness=config.get("strandedness", "reverse"),
        coverage_cutoff=config.get("coverage_cutoff", 10),
    threads: 8
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})" logs/bam_to_3pend_bg

        T={params.coverage_cutoff}
        strand="{wildcards.strand}"
        library="{params.strandedness}"

        if [[ "$strand" == "plus" ]]; then
            if [[ "$library" == "reverse" ]]; then
                strand_symbol="-"
                pos_symbol="-5"
            else
                strand_symbol="+"
                pos_symbol="-3"
            fi
        elif [[ "$strand" == "minus" ]]; then
            if [[ "$library" == "reverse" ]]; then
                strand_symbol="+"
                pos_symbol="-5"
            else
                strand_symbol="-"
                pos_symbol="-3"
            fi
        else
            echo "ERROR: strand must be plus or minus; got: $strand" >&2
            exit 1
        fi

        samtools sort -@ {threads} {input} 2> {log} \
        | genomeCoverageBed -ibam - -bg -strand "$strand_symbol" "$pos_symbol" 2>> {log} \
        | awk -v T="$T" -v OFS='\t' '$4>=T' \
        | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """


rule merge_3pend_bg:
    input:
        expand("results/bam2bg/{sample}_informative_{{strand}}.bg", sample=SAMP_NAMES),
    output:
        temp("results/merge_3pend_bg/merged_3pend_{strand}.bg"),
    params:
        strandedness=config.get("strandedness", "reverse"),
    envmodules:
        config["modules"]["bedtools"],
    log:
        "logs/merge_3pend_bg/{strand}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})" logs/merge_3pend_bg

        strand="{wildcards.strand}"
        library="{params.strandedness}"

        if [[ "$strand" == "plus" ]]; then
            if [[ "$library" == "reverse" ]]; then strand_symbol="-"; else strand_symbol="+"; fi
        elif [[ "$strand" == "minus" ]]; then
            if [[ "$library" == "reverse" ]]; then strand_symbol="+"; else strand_symbol="-"; fi
        else
            echo "ERROR: strand must be plus or minus; got: $strand" >&2
            exit 1
        fi

        bedtools unionbedg -i {input} \
        | awk 'BEGIN{{OFS="\t"}}{{s=0; for(i=4;i<=NF;i++) s+=$i; print $1,$2,$3,s}}' \
        | LC_COLLATE=C sort -k1,1 -k2,2n > {output} 2> {log}
        """


rule cluster_PAS_bedtools:
    input:
        "results/merge_3pend_bg/merged_3pend_{strand}.bg",
    output:
        temp("results/call_PAS/bedtools_clusters_{strand}.bg"),
    log:
        "logs/call_PAS_bedtools/bedtools_cluster_{strand}.log",
    params:
        extra=config.get("bedtools_cluster_distance", 50),
    envmodules:
        config["modules"]["bedtools"],
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})" logs/call_PAS_bedtools
        bedtools cluster -d {params.extra} -i {input} > {output} 2> {log}
        """


rule summarise_PAS_clusters:
    input:
        "results/call_PAS/bedtools_clusters_{strand}.bg",
    output:
        "results/summarise_PAS_clusters/summarise_clusters_{strand}.bg",
    log:
        "logs/summarise_PAS_clusters/{strand}.log",
    envmodules:
        config["modules"]["bedtools"],
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})" logs/summarise_PAS_clusters
        bedtools groupby -g 1,5 -c 2,3,4 -o min,max,sum -i {input} > {output} 2> {log}
        """


rule make_threepUTR_gtf:
    input:
        bed=expand(
            "results/summarise_PAS_clusters/summarise_clusters_{strand}.bg",
            strand=["minus", "plus"],
        ),
        gtf=config.get("annotation"),
        fasta=config.get("genome"),
    output:
        "annotations/threepUTR_annotation.gtf",
    log:
        "logs/make_threepUTR_gtf/make_threepUTR_gtf.log",
    envmodules:
        config["modules"]["r"],
        config["modules"]["bioconductor"],
    params:
        rscript=workflow.source_path("../scripts/threeputrs/make_3pgtf.R"),
        coverage=config.get("cluster_coverage", 20 * NUM_SAMPS),
        fxn=config.get("cluster_fxn", 0.0),
        extension=config.get("extension", 0),
        polyA=config.get("false_polyA_len", 7),
        CPA=config.get("require_CPA_site", False),
        only_annotated=config.get("only_annotated_threeputrs", False),
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})" logs/make_threepUTR_gtf
        chmod +x {params.rscript}
        {params.rscript} \
            --bed_minus {input.bed[0]} \
            --bed_plus {input.bed[1]} \
            --gtf {input.gtf} \
            --fasta {input.fasta} \
            --output {output} \
            --extension {params.extension} \
            --min_coverage {params.coverage} \
            --false_polyA_len {params.polyA} \
            --require_CPA {params.CPA} \
            --only_annotated {params.only_annotated} \
            --min_fxn {params.fxn} 1> {log} 2>&1
        """
