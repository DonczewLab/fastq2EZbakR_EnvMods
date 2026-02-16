# fastq2EZbakR_EnvMods
Environment-module–based execution of fastq2EZbakR

This repository is a reproducible re-implementation of the original  
https://github.com/isaacvock/fastq2EZbakR Snakemake workflow that **replaces Conda/Singularity execution with HPC environment modules**, while preserving the original pipeline structure, scripts, and configuration semantics as closely as possible.

All scientific logic, scripts, and parameterization are inherited from the original workflow unless explicitly stated otherwise.

---

## Key differences from upstream

+ No Conda environments
+ No Singularity containers
+ Uses `snakemake --use-envmodules`
+ Centralized `modules.yaml` for software versions
+ Original scripts remain unmodified
+ Snakemake rules minimally edited to remove `conda:` blocks
+ Wrapper usage preserved where compatible (e.g. `fastp`)

This approach is intended for shared HPC systems where Conda and Singularity are unavailable or discouraged.

---

# Project Description
A Snakemake pipeline designed to process nucleotide recoding RNA-seq data (NR-seq, e.g., TimeLapse-seq, SLAM-seq, TUC-seq, etc.). fastq2EZbakR provides output readily compatible with EZbakR, but is also designed to provide processed NR-seq data in a convenient form. fastq2EZbakR_EnvMods takes the same original workflow but converts to env mods 

### Attribution
This workflow is adapted from [fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR), originally developed by Isaac Vock under the Apache License 2.0. Modifications by Kevin A. Boyd (OMRF) include conversion from Conda to HPC environment modules, refactoring of the Snakemake workflow, and OMRF-specific cluster configuration files.

---

# Sample Configuration Guide

This section explains **exactly how to define samples** for use with this pipeline.  
No other files need to be modified when changing samples, provided paths and parameters are valid.

---

## 1. Sample identifiers

Samples are defined by **sample IDs**, which are used consistently across:

+ directory names
+ output files
+ configuration dictionaries (`pnews`, `polds`, etc.)

Example sample IDs:

```yaml
WT_1
WT_2
WT_ctl
KO_1
KO_2
KO_ctl
```

Sample ID rules:
+ Must be unique
+ Should not contain spaces
+ Filesystem-safe characters only (letters, numbers, underscores recommended)

---

## 2. Required sample-level parameters

Each sample must appear in the following config sections.

pnews
```yaml
pnews:
  WT_1: 0.5
  WT_2: 0.5
  WT_ctl: 0.5
  KO_1: 0.5
  KO_2: 0.5
  KO_ctl: 0.5
```
polds
```yaml
polds:
  WT_1: 0.5
  WT_2: 0.5
  WT_ctl: 0.5
  KO_1: 0.5
  KO_2: 0.5
  KO_ctl: 0.5
```

+ Every sample ID must appear in both dictionaries
+ Values are passed directly to EZbakR modeling steps

---

## 3. Control samples

Control samples must be explicitly listed.

```yaml
control_samples:
  - WT_ctl
  - KO_ctl
```

+ Sample IDs must match exactly
+ Used for EZbakR model specification

---

## 4. Input data modes
The pipeline supports two mutually exclusive input modes.

### A. FASTQ input (default)
Used when starting from raw sequencing data.

```yaml
bam2bakr: False
```

Each sample maps to a directory containing FASTQ files.

```yaml
samples:
  WT_1: data/fastq/WT_1
  WT_2: data/fastq/WT_2
  WT_ctl: data/fastq/WT_ctl
  KO_1: data/fastq/KO_1
  KO_2: data/fastq/KO_2
  KO_ctl: data/fastq/KO_ctl
FASTQ directory expectations
```
Paired-end data (`PE: True`)
```yaml
*_R1*.fastq.gz
*_R2*.fastq.gz
```
Single-end data (`PE: False`)
```yaml
*.fastq.gz
```
+ Both compressed and uncompressed FASTQs are supported.

### B. BAM input (skip preprocessing)
Used when alignments already exist.
```yaml
bam2bakr: True
```
Each sample maps directly to a BAM file.
```yaml
samples:
  WT_1: data/bams/WT_1.bam
  WT_2: data/bams/WT_2.bam
  WT_ctl: data/bams/WT_ctl.bam
  KO_1: data/bams/KO_1.bam
  KO_2: data/bams/KO_2.bam
  KO_ctl: data/bams/KO_ctl.bam
```

When `bam2bakr: True:`
+ Adapter trimming is skipped
+ Alignment is skipped
+ Index generation is skipped

---

## 5. Paired-end vs single-end
Set globally in the config file.

`PE: True`
or
`PE: False`

+ All samples are assumed to follow the same layout

---

## 6. Genome and annotation
The following files must exist before running the pipeline.
```yaml
genome: data/genome/genome.fasta
annotation: data/annotation/genome.gtf
```
+ Alignment indices will be created automatically if not present.
```
indices: data/indices/star_index
```

---

## 7. Adding or changing samples
To add new samples:

Define new sample IDs

Add them to:

+ `pnews`
+ `polds`
+ `samples`
+ Update `control_samples` if applicable

Place FASTQs or BAMs in the specified locations
No other pipeline files need to be edited

---

## 8. Environment modules (modified from original workflow)
All software dependencies are defined in:
```yaml
config/modules.yaml
```
Snakemake rules load modules using:
```yaml
envmodules:
  - tool_name
```

This allows:
+ Centralized software version control
+ Reproducible execution on shared HPC systems
+ Compatibility with snakemake --use-envmodules

---

## 9. Running the pipeline (modified from original workflow)
### Instructions to run on Slurm managed HPC  
9A. Download version controlled repository
```
git clone https://github.com/DonczewLab/fastq2EZbakR_EnvMods.git
```
9B. Load modules
```
module purge
module load slurm python/3.10
```
9C. Modify config file
```
vim config/config.yml
```
9D. Dry Run
```
snakemake -npr
```
9E. Run on HPC with config.yml options
```
sbatch submit_fastq2ezbark.sh
```

---

## 10. Licensing

This repository contains workflow infrastructure, configuration adaptations,
and HPC environment-module integration developed by the Donczew Lab.

These additions are licensed under the Apache License 2.0 — see the [LICENSE](LICENSE) file for details

This repository is designed to execute the fastq2EZbakR workflow developed by
Vock et al. The original fastq2EZbakR source code is not relicensed here and
remains subject to its original terms and authorship.

Users should obtain the original fastq2EZbakR software from:
https://github.com/isaacvock/fastq2EZbakR
[fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR)
