# Extract the vignette subset from the processed baseline data.
# Usage: Rscript inst/scripts/make_real_data.R betas.rds meta.rds [output_dir]
# Source: Ammala et al. (2025), https://doi.org/10.1111/jsr.14438
# Original data: https://osf.io/qx7w3/
# Inputs are the preprocessed baseline beta matrix and sample metadata.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L))
    stop("Supply betas.rds, meta.rds, and optionally an output directory.")
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
output_dir <- if (length(args) == 3L) args[3] else
    file.path(script_dir, "..", "extdata")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

beta <- readRDS(args[1])
samples <- readRDS(args[2])
sample_index <- match(colnames(beta), as.character(samples$Sample_Name))
if (anyNA(sample_index) || anyDuplicated(samples$Sample_Name))
    stop("Sample metadata must uniquely match the beta-matrix columns.")
samples <- samples[sample_index, , drop = FALSE]

# Select loci before fitting, retaining all observed baseline samples.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(2026)
selected <- sort(sample.int(nrow(beta), 2000L))
sample_info <- data.frame(
    time = samples$time - min(samples$time),
    subject = samples$id,
    source_time = samples$time,
    row.names = colnames(beta))
beta_file <- gzfile(file.path(output_dir, "real_beta.csv.gz"), "wt")
write.csv(beta[selected, , drop = FALSE], beta_file)
close(beta_file)
write.csv(sample_info, file.path(output_dir, "real_samples.csv"))
