### THESE RULES PERTAIN TO THE ASSIGNMENT OF READS TO FEATURES WITH FEATURECOUNTS

# Assign reads to genes
rule featurecounts_genes:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_genes/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_genes/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_genes_extra"] + FC_GENES_PARAMS,
        outprefix="results/featurecounts_genes/{sample}",
    log:
        "logs/featurecounts_genes/{sample}.log",
    envmodules:
        config["modules"]["subread"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/featurecounts_genes
        featureCounts \
          -T {threads} \
          -s {params.strand} \
          -a {input.annotation} \
          -o {params.outprefix}.featureCounts \
          {params.extra} \
          {input.samples} \
          1> {log} 2>&1

        # REQUIRED downstream: featureCounts -R CORE output (DO NOT overwrite it)
        core="results/featurecounts_genes/{wildcards.sample}.s.bam.featureCounts"
        test -s "${{core}}"
        """


# Assign reads to exons
rule featurecounts_exons:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_exons/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
        temp("results/featurecounts_exons/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_exons_extra"] + FC_EXONS_PARAMS,
        outprefix="results/featurecounts_exons/{sample}",
    log:
        "logs/featurecounts_exons/{sample}.log",
    envmodules:
        config["modules"]["subread"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/featurecounts_exons
        featureCounts \
          -T {threads} \
          -s {params.strand} \
          -a {input.annotation} \
          -o {params.outprefix}.featureCounts \
          {params.extra} \
          {input.samples} \
          1> {log} 2>&1

        # REQUIRED downstream: featureCounts -R CORE output (DO NOT overwrite it)
        core="results/featurecounts_exons/{wildcards.sample}.s.bam.featureCounts"
        test -s "${{core}}"
        """


# Assign reads to transcripts
rule featurecounts_transcripts:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_transcripts/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_transcripts/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_transcripts_extra"] + FC_TRANSCRIPTS_PARAMS,
        outprefix="results/featurecounts_transcripts/{sample}",
    log:
        "logs/featurecounts_transcripts/{sample}.log",
    envmodules:
        config["modules"]["subread"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/featurecounts_transcripts
        featureCounts \
          -T {threads} \
          -s {params.strand} \
          -a {input.annotation} \
          -o {params.outprefix}.featureCounts \
          {params.extra} \
          {input.samples} \
          1> {log} 2>&1

        # REQUIRED downstream: featureCounts -R CORE output (DO NOT overwrite it)
        core="results/featurecounts_transcripts/{wildcards.sample}.s.bam.featureCounts"
        test -s "${{core}}"
        """


# Assign reads to exonic bins
rule featurecounts_exonbins:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["flat_annotation"],
    output:
        multiext(
            "results/featurecounts_exonbins/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_exonbins/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_exonbins_extra"] + FC_EXONBINS_PARAMS,
        outprefix="results/featurecounts_exonbins/{sample}",
    log:
        "logs/featurecounts_exonbins/{sample}.log",
    envmodules:
        config["modules"]["subread"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/featurecounts_exonbins
        featureCounts \
          -T {threads} \
          -s {params.strand} \
          -a {input.annotation} \
          -o {params.outprefix}.featureCounts \
          {params.extra} \
          {input.samples} \
          1> {log} 2>&1

        # REQUIRED downstream: featureCounts -R CORE output (DO NOT overwrite it)
        core="results/featurecounts_exonbins/{wildcards.sample}.s.bam.featureCounts"
        test -s "${{core}}"
        """


# Get the set of isoforms a read maps to from transcriptome bam
# (same script used regardless of rsem; only the input bam path differs)
if config.get("run_rsem", True):

    rule read_to_transcripts:
        input:
            bam="results/rsem/{sample}.transcript.bam",
        output:
            table=temp("results/read_to_transcripts/{sample}.csv"),
        log:
            "logs/read_to_transcripts/{sample}.log",
        threads: 1
        envmodules:
            config["modules"]["python"],
        script:
            "../scripts/features/transcript_assignment.py"

else:

    rule read_to_transcripts:
        input:
            bam="results/align/{sample}-Aligned.toTranscriptome.out.bam",
        output:
            table=temp("results/read_to_transcripts/{sample}.csv"),
        log:
            "logs/read_to_transcripts/{sample}.log",
        threads: 1
        envmodules:
            config["modules"]["python"],
        script:
            "../scripts/features/transcript_assignment.py"


# Get set of junctions a read overlaps
rule read_to_junctions:
    input:
        "results/sf_reads/{sample}.s.bam",
    output:
        temp("results/read_to_junctions/{sample}.csv.gz"),
        temp("results/read_to_junctions/{sample}_check.txt"),
    params:
        shellscript=workflow.source_path("../scripts/features/junction_assignment.sh"),
        pythonscript=workflow.source_path("../scripts/features/junction_assignment.py"),
        awkscript=workflow.source_path("../scripts/rsem_plus/fragment_sam_rsem.awk"),
    log:
        "logs/read_to_junctions/{sample}.log",
    threads: 32
    envmodules:
        # This step typically uses samtools + python + coreutils/awk.
        # awk is usually system; samtools + python are safe to load.
        config["modules"]["samtools"],
        config["modules"]["python"],
    shell:
        r"""
        set -euo pipefail
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        chmod +x {params.awkscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.pythonscript} {params.awkscript}  1> {log} 2>&1
        """


# Make junction annotation that featureCounts can assign reads with respect to
rule junction_annotation:
    input:
        config["annotation"],
    output:
        "junction_annotation/junctions.gtf",
    params:
        rscript=workflow.source_path("../scripts/features/junction_annotation.R"),
        extra=config["junction_annotation_params"],
    log:
        "logs/junction_annotation/junctions.log",
    threads: 1
    envmodules:
        config["modules"]["r"],
        config["modules"]["bioconductor"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p junction_annotation logs/junction_annotation
        chmod +x {params.rscript}
        {params.rscript} -r {input} -o {output} {params.extra} 1> {log} 2>&1
        """


# Experimental: assign reads to exon-exon junctions
rule featurecounts_eej:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation="junction_annotation/junctions.gtf",
    output:
        multiext(
            "results/featurecounts_eej/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_eej/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_eej_extra"] + FC_EEJ_PARAMS,
        outprefix="results/featurecounts_eej/{sample}",
    log:
        "logs/featurecounts_eej/{sample}.log",
    envmodules:
        config["modules"]["subread"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/featurecounts_eej
        featureCounts \
          -T {threads} \
          -s {params.strand} \
          -a {input.annotation} \
          -o {params.outprefix}.featureCounts \
          {params.extra} \
          {input.samples} \
          1> {log} 2>&1

        # REQUIRED downstream: featureCounts -R CORE output (DO NOT overwrite it)
        core="results/featurecounts_eej/{wildcards.sample}.s.bam.featureCounts"
        test -s "${{core}}"
        """


# Experimental: assign reads to exon-intron junctions
rule featurecounts_eij:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation="junction_annotation/junctions.gtf",
    output:
        multiext(
            "results/featurecounts_eij/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_eij/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_eij_extra"] + FC_EIJ_PARAMS,
        outprefix="results/featurecounts_eij/{sample}",
    log:
        "logs/featurecounts_eij/{sample}.log",
    envmodules:
        config["modules"]["subread"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/featurecounts_eij
        featureCounts \
          -T {threads} \
          -s {params.strand} \
          -a {input.annotation} \
          -o {params.outprefix}.featureCounts \
          {params.extra} \
          {input.samples} \
          1> {log} 2>&1

        # REQUIRED downstream: featureCounts -R CORE output (DO NOT overwrite it)
        core="results/featurecounts_eij/{wildcards.sample}.s.bam.featureCounts"
        test -s "${{core}}"
        """


# Assign reads to 3'-UTRs
rule featurecounts_3utr:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=THREEPUTR_ANNOTATION,
    output:
        multiext(
            "results/featurecounts_3utr/{sample}",
            ".featureCounts",
            ".featureCounts.summary"
        ),
        temp("results/featurecounts_3utr/{sample}.s.bam.featureCounts"),
    params:
        strand=FC_STRAND,
        extra=FC_3UTR_PARAMS,
        outprefix="results/featurecounts_3utr/{sample}",
    threads: 20
    log:
        "logs/featurecounts_3utr/{sample}.log",
    envmodules:
        config["modules"]["subread"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/featurecounts_3utr
        featureCounts \
          -T {threads} \
          -s {params.strand} \
          -a {input.annotation} \
          -o {params.outprefix}.featureCounts \
          {params.extra} \
          {input.samples} \
          1> {log} 2>&1

        # REQUIRED downstream: featureCounts -R CORE output (DO NOT overwrite it)
        core="results/featurecounts_3utr/{wildcards.sample}.s.bam.featureCounts"
        test -s "${{core}}"
        """
