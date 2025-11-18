#!/usr/bin/env python3
import argparse
import os
import gzip

ap = argparse.ArgumentParser()
ap.add_argument("--merged-dir", required=True)
ap.add_argument("--finalize-cB", required=True)
ap.add_argument("--finalize-cUP", required=True)
ap.add_argument("--finalize-arrow", required=True)
ap.add_argument("--lowram", required=True)
ap.add_argument("--out-cb", required=True)
ap.add_argument("--out-cup", required=True)
ap.add_argument("--arrow-dir", required=True)
args = ap.parse_args()

os.makedirs(os.path.dirname(args.out_cb), exist_ok=True)
os.makedirs(os.path.dirname(args.out_cup), exist_ok=True)
os.makedirs(args.arrow_dir, exist_ok=True)

if int(args.finalize_cB):
    with gzip.open(args.out_cb, "wb") as fo:
        fo.write(b"feature,count\nstub,0\n")

if int(args.finalize_cUP):
    with gzip.open(args.out_cup, "wb") as fo:
        fo.write(b"feature,count\nstub,0\n")

if int(args.finalize_arrow):
    open(os.path.join(args.arrow_dir, "_DONE"), "w").close()
