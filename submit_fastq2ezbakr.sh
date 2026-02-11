#!/usr/bin/env bash
#SBATCH -A donczew-lab
#SBATCH -p serial
#SBATCH -t 08:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH -J fastq2EZbakR.master
#SBATCH -o logs/slurm-master-%j.out
#SBATCH -e logs/slurm-master-%j.err

set -euo pipefail
mkdir -p logs

module load python/3.10

snakemake \
  -j 20 \
  --use-envmodules \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 600 \
  --printshellcmds \
  --cluster-config config/cluster_config.yaml \
  --cluster "sbatch -A {cluster.account} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} -J {cluster.name} -o {cluster.output} -e {cluster.error}"
