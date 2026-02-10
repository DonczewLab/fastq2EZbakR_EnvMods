### Little bit of Snakemake smoving here. If downloading
### fastq files, SAMP_NAMES is set equal to SRA accession IDs
### provided in config rather than the sample IDs specified under
### samples in the config. This allows for the {accession} wild
### card to be determined.
if config["PE"]:

    rule download_fastq:
        output:
            "results/download_fastq/{accession}_1.fastq",
            "results/download_fastq/{accession}_2.fastq",
        params:
            extra=config["fasterq_dump_extras"],
        envmodules:
            config["modules"]["ncbi_sra"],
        threads: 12
        log:
            "logs/download_fastq/{accession}.log",
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output[0]})" logs/download_fastq

            tmpdir="$(mktemp -d)"
            trap 'rm -rf "$tmpdir"' EXIT

            outdir="$(dirname {output[0]})"
            fasterq-dump \
              --temp "$tmpdir" \
              --threads {threads} \
              {params.extra} \
              --outdir "$outdir" \
              "{wildcards.accession}" \
              1> {log} 2>&1

            # Wrapper-like compression auto-detect:
            for f in {output}; do
              case "$f" in
                *.gz)
                  base="${f%.gz}"
                  test -s "$base"
                  pigz -p {threads} -f "$base"
                  test -s "$f"
                  ;;
                *.bz2)
                  base="${f%.bz2}"
                  test -s "$base"
                  # pbzip2 supports -pN and optional -m (mem); wrapper may add mem, but we omit.
                  pbzip2 -p{threads} -f "$base"
                  test -s "$f"
                  ;;
                *)
                  test -s "$f"
                  ;;
              esac
            done
            """

else:

    rule download_fastq:
        output:
            "results/download_fastq/{accession}.fastq",
        params:
            extra=config["fasterq_dump_extras"],
        envmodules:
            config["modules"]["ncbi_sra"],
        threads: 12
        log:
            "logs/download_fastq/{accession}.log",
        shell:
            r"""
            set -euo pipefail
            mkdir -p "$(dirname {output})" logs/download_fastq

            tmpdir="$(mktemp -d)"
            trap 'rm -rf "$tmpdir"' EXIT

            outdir="$(dirname {output})"
            fasterq-dump \
              --temp "$tmpdir" \
              --threads {threads} \
              {params.extra} \
              --outdir "$outdir" \
              "{wildcards.accession}" \
              1> {log} 2>&1

            # Wrapper-like compression auto-detect:
            f="{output}"
            case "$f" in
              *.gz)
                base="${f%.gz}"
                test -s "$base"
                pigz -p {threads} -f "$base"
                test -s "$f"
                ;;
              *.bz2)
                base="${f%.bz2}"
                test -s "$base"
                pbzip2 -p{threads} -f "$base"
                test -s "$f"
                ;;
              *)
                test -s "$f"
                ;;
            esac
            """
