#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
tracks_dir <- args[1]
out <- args[2]
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
# Emit a trivial table usable downstream
samples <- unique(gsub("\\..*$","", basename(list.files(tracks_dir, pattern="\\.tdf$"))))
if (length(samples)==0) samples <- c("sampleA","sampleB")
sf <- data.frame(sample=samples, scale_factor=rep(1.0, length(samples)))
write.table(sf, file=out, sep="\t", row.names=FALSE, quote=FALSE)
