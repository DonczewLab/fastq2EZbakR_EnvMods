### Compile all QC into one report
rule multiqc:
    input:
        MULTIQC_INPUT,
    output:
        "results/multiqc/multiqc_report.html",
    params:
        extra=config.get("multiqc_extra", ""),
        outdir="results/multiqc",
        # MultiQC works best when you point it at directories, not thousands of files
        input_dirs=lambda wc, input: sorted({os.path.dirname(fp) for fp in input}),
        report_name="multiqc_report.html",
    log:
        "logs/multiqc/multiqc.log",
    threads: 1
    envmodules:
        config["modules"]["multiqc"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        multiqc {params.extra} --force \
          -o {params.outdir} -n {params.report_name} \
          {params.input_dirs} \
          1> {log} 2>&1
        """
