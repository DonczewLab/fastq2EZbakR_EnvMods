#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
tracks_dir <- args[1]
out <- args[2]

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

# Stub: set all scale factors to 1.0
files <- list.files(tracks_dir, pattern = "\\.tdf$", full.names = FALSE)

if (length(files) == 0) {
  samples <- c("sampleA", "sampleB")
} else {
  # sample names = everything before first "."
  samples <- unique(sub("\\..*$", "", files))
}

sf <- data.frame(sample = samples,
                 scale_factor = rep(1.0, length(samples)))

write.table(sf, file = out, sep = "\t",
            row.names = FALSE, quote = FALSE)
