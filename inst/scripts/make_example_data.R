# Script to generate example data for RhyME package.
# Uses simulation rather than subsetting real data for portability.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
package_root <- if (length(file_arg)) {
    normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."))
} else if (file.exists("DESCRIPTION")) {
    normalizePath(".")
} else {
    normalizePath(file.path("packages", "RhyME"))
}
out_dir <- file.path(package_root, "inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
n_cpg <- 3000
n_subjects <- 8
n_timepoints <- 6
n_samples <- n_subjects * n_timepoints
n_rhythmic <- 600

times <- rep(seq(0, 20, length.out = n_timepoints), n_subjects)
subjects <- rep(paste0("S", seq_len(n_subjects)), each = n_timepoints)
omega <- 2 * pi / 24

beta_mat <- matrix(NA_real_, nrow = n_cpg, ncol = n_samples)
cpg_names <- paste0("cg", sprintf("%05d", seq_len(n_cpg)))
sample_names <- paste0(subjects, "_T", sprintf("%.0f", times))
rownames(beta_mat) <- cpg_names
colnames(beta_mat) <- sample_names

for (i in seq_len(n_cpg)) {
    b0 <- rnorm(1, mean = 0, sd = 0.5)
    if (i <= n_rhythmic) {
        b_sin <- rnorm(1, mean = 0.3, sd = 0.08)
        b_cos <- rnorm(1, mean = 0.3, sd = 0.08)
    } else {
        b_sin <- 0
        b_cos <- 0
    }
    tau <- 0.2
    phi <- exp(rnorm(1, mean = 3, sd = 0.3))

    u <- rnorm(n_subjects, 0, tau)
    u_expanded <- rep(u, each = n_timepoints)

    eta <- b0 + b_sin * sin(omega * times) + b_cos * cos(omega * times) +
        u_expanded
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

saveRDS(example_data, file.path(out_dir, "example_data.rds"))

cat("Example data saved to inst/extdata/example_data.rds\n")
cat(sprintf("  %d CpGs, %d samples (%d subjects x %d timepoints)\n",
            n_cpg, n_samples, n_subjects, n_timepoints))
cat(sprintf("  %d rhythmic, %d null\n", n_rhythmic, n_cpg - n_rhythmic))
