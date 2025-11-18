#!/usr/bin/env bash
set -euo pipefail

bamdir=""
outdir=""
wsl=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir) bamdir="$2"; shift 2;;
    --out-dir) outdir="$2"; shift 2;;
    --wsl)     wsl="$2"; shift 2;;
    *) echo "Unknown arg $1" >&2; exit 1;;
  esac
done

mkdir -p "$outdir"

# Stub implementation:
# Discover samples by bam basename before first dot, and just touch the
# expected TDF files so downstream rules see them.
for b in "$bamdir"/*.bam; do
  [ -e "$b" ] || continue
  s="$(basename "$b")"
  s="${s%%.*}"     # strip extension & suffix
  for tr in TC; do
    for strand in pos neg; do
      for idx in 0 1 2; do
        touch "$outdir/${s}.${tr}.${idx}.${strand}.tdf"
      done
    done
  done
done
