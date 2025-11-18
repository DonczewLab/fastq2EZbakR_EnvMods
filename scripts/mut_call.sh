#!/usr/bin/env bash
set -euo pipefail

bam=""
snps=""
out=""
tracks="TC"
threads=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam)        bam="$2"; shift 2;;
    --snps)       snps="$2"; shift 2;;
    --out)        out="$2"; shift 2;;
    --mut-tracks) tracks="$2"; shift 2;;
    --threads)    threads="$2"; shift 2;;
    *) echo "Unknown arg $1" >&2; exit 1;;
  esac
done

mkdir -p "$(dirname "$out")"

# Stub: emit a tiny valid CSV.GZ to satisfy downstream rules
tmp="$(mktemp)"
{
  echo "chrom,start,end,sample,track,count"
  echo "chr1,1000,1001,${bam##*/},${tracks},0"
} > "$tmp"
gzip -c "$tmp" > "$out"
rm -f "$tmp"
