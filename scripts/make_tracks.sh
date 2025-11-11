#!/usr/bin/env bash
set -euo pipefail
bamdir=""
outdir=""
wsl=0
while [[ $# -gt 0 ]]; do
  case $1 in
    --bam-dir) bamdir="$2"; shift 2;;
    --out-dir) outdir="$2"; shift 2;;
    --wsl) wsl="$2"; shift 2;;
    *) echo "Unknown arg $1"; exit 1;;
  esac
done

mkdir -p "$outdir"
# Discover samples by bam basename before first dot
for b in "$bamdir"/*.bam; do
  s=$(basename "$b"); s="${s%%.*}"
  for tr in TC; do
    for strand in pos neg; do
      for idx in 0 1 2; do
        touch "$outdir/${s}.${tr}.${idx}.${strand}.tdf"
      done
    end
  done
done
