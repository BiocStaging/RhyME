# Generate synthetic observations using real-data parameter estimates.
# Usage: Rscript inst/scripts/make_example_data.R /path/to/result_PQL_ML.rds
# Reference: full-model, unpenalized PQL fits from the 2026-09-16 real-data
# analysis (811,808 CpGs, 112 samples, 15 subjects; reference MD5
# ac74dc64d109a3682979daa973533c9f). Parameter rows are sampled jointly.
# Rhythmic rows use BH-significant Wald reference fits (padj_pql_ml < 0.05);
# null rows use the complete usable pool with sine/cosine coefficients set to 0.
# The 20% rhythmic fraction and balanced 8-subject design are example settings.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
package_root <- if (length(file_arg)) {
    normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."))
} else if (file.exists("DESCRIPTION")) {
    normalizePath(".")
} else {
    stop("Run this script from the package directory or with Rscript.")
}
out_dir <- file.path(package_root, "inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reference_file <- commandArgs(trailingOnly = TRUE)
if (length(reference_file) != 1L)
    stop("Supply the reference result_PQL_ML.rds file as the only argument.")
reference <- readRDS(reference_file)
fields <- c("b0", "b_sin", "b_cos", "log_phi", "tau")
valid <- reference$converged & reference$tau >= 0 &
    Reduce(`&`, lapply(reference[fields], is.finite))
pool <- reference[which(valid), ]
rhythmic <- which(is.finite(pool$padj_pql_ml) & pool$padj_pql_ml < 0.05)

set.seed(42)
n_cpg <- 10000
n_subjects <- 8
n_timepoints <- 6
n_samples <- n_subjects * n_timepoints
n_rhythmic <- 2000

index <- c(sample(rhythmic, n_rhythmic, replace = TRUE),
    sample(seq_len(nrow(pool)), n_cpg - n_rhythmic, replace = TRUE))
parameters <- pool[index, fields]
parameters[(n_rhythmic + 1):n_cpg, c("b_sin", "b_cos")] <- 0

times <- rep(seq(0, 20, length.out = n_timepoints), n_subjects)
subjects <- rep(paste0("S", seq_len(n_subjects)), each = n_timepoints)
omega <- 2 * pi / 24

beta_mat <- matrix(NA_real_, nrow = n_cpg, ncol = n_samples)
cpg_names <- paste0("cg", sprintf("%05d", seq_len(n_cpg)))
sample_names <- paste0(subjects, "_T", sprintf("%.0f", times))
rownames(beta_mat) <- cpg_names
colnames(beta_mat) <- sample_names

for (i in seq_len(n_cpg)) {
    phi <- exp(parameters$log_phi[i])
    u <- rnorm(n_subjects, 0, parameters$tau[i])
    u_expanded <- rep(u, each = n_timepoints)

    eta <- parameters$b0[i] + parameters$b_sin[i] * sin(omega * times) +
        parameters$b_cos[i] * cos(omega * times) + u_expanded
    mu <- plogis(eta)
    a <- mu * phi
    b <- (1 - mu) * phi

    y <- rbeta(n_samples, a, b)
    beta_mat[i, ] <- pmin(pmax(y, 1e-4), 1 - 1e-4)
}

meta <- data.frame(
    time = times,
    subject = subjects,
    Sample_Name = sample_names,
    stringsAsFactors = FALSE
)
rownames(meta) <- sample_names

se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(beta = beta_mat),
    colData = S4Vectors::DataFrame(meta)
)

example_data <- list(
    se = se,
    sample_meta = meta
)

saveRDS(example_data, file.path(out_dir, "example_data.rds"), compress = "xz")

cat("Example data saved to inst/extdata/example_data.rds\n")
cat(sprintf("  %d CpGs, %d samples (%d subjects x %d timepoints)\n",
            n_cpg, n_samples, n_subjects, n_timepoints))
cat(sprintf("  %d rhythmic, %d null\n", n_rhythmic, n_cpg - n_rhythmic))
cat(sprintf("  Reference: %d usable fits, %d BH-significant fits\n",
    nrow(pool), length(rhythmic)))
cat("  Reference MD5:", unname(tools::md5sum(reference_file)), "\n")
