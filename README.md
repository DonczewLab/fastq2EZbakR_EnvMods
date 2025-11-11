# fastq2EZbakR_EnvMods

## Description
A Snakemake pipeline designed to process nucleotide recoding RNA-seq data (NR-seq, e.g., TimeLapse-seq, SLAM-seq, TUC-seq, etc.). fastq2EZbakR provides output readily compatible with EZbakR, but is also designed to provide processed NR-seq data in a convenient form. fastq2EZbakR_EnvMods takes the same original workflow but converts to env mods 

## Attribution
This workflow is adapted from [fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR),
originally developed by Isaac Vock under the Apache License 2.0.
Modifications by Kevin A. Boyd (OMRF) include conversion from Conda to
HPC environment modules, refactoring of the Snakemake workflow,
and OMRF-specific cluster configuration files.

## License
This project is licensed under the Apache License 2.0 — see the [LICENSE](LICENSE) file for details.

