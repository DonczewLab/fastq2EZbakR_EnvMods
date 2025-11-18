#!/usr/bin/env python3
import argparse
import os
import gzip
import glob
import json

ap = argparse.ArgumentParser()
ap.add_argument("--counts-dir", required=True)
ap.add_argument("--features", required=True)   # JSON/string, unused in stub
ap.add_argument("--annotation", required=True) # path, unused in stub
ap.add_argument("--out-dir", required=True)
ap.add_argument("--lowram", required=True)     # "0"/"1", unused in stub
args = ap.parse_args()

os.makedirs(args.out_dir, exist_ok=True)

for c in sorted(glob.glob(os.path.join(args.counts_dir, "*_counts.csv.gz"))):
    base = os.path.basename(c)
    out = os.path.join(args.out_dir, base)
    with gzip.open(c, "rb") as fi, gzip.open(out, "wb") as fo:
        fo.write(fi.read())
