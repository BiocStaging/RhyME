# Regenerate the vignette data from the bundled parameters and sampling design.
# Usage: Rscript inst/scripts/make_example_data.R [output.rds]
# See the adjacent README for the input provenance and simulation model.

library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply at most one output path.")
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
data_dir <- file.path(script_dir, "..", "extdata")
input <- readRDS(file.path(data_dir, "example_parameters.rds"))
parameters <- input$parameters
design <- input$design
provenance <- input$provenance
design$sin24 <- sin(2 * pi * design$time / provenance$period)
design$cos24 <- cos(2 * pi * design$time / provenance$period)
do.call(RNGkind, as.list(provenance$rng))
output <- if (length(args)) args[1] else
    file.path(data_dir, "example_data.rds")

# Repeated measurements: 15 subjects, 112 samples, 2,000 CpGs.
set.seed(provenance$observation_seed)
beta <- do.call(rbind, lapply(seq_len(nrow(parameters)), function(g) {
    u <- rnorm(nlevels(design$id), sd = parameters$tau[g])
    eta <- parameters$b0[g] + parameters$b_sin[g] * design$sin24 +
        parameters$b_cos[g] * design$cos24 + u[as.integer(design$id)]
    mu <- plogis(eta)
    y <- rbeta(nrow(design), mu * parameters$phi[g],
        (1 - mu) * parameters$phi[g])
    pmin(pmax(y, 1e-4), 1 - 1e-4)
}))
dimnames(beta) <- list(parameters$cpg, design$Sample_Name)
se <- SummarizedExperiment(assays = list(beta = beta),
    rowData = DataFrame(true_status = ifelse(parameters$true,
        "Rhythmic", "Non-rhythmic"),
        true_logit_amplitude = parameters$amplitude_link,
        true_phase = parameters$phase_true),
    colData = DataFrame(time = design$time, subject = design$id,
        row.names = design$Sample_Name))
metadata(se)$simulation <- provenance

# Independent observations for the default fixed-effects example.
p <- parameters[seq_len(1200), ]
time <- rep(seq(0, 20, by = 4), 8)
mu <- plogis(p$b0 + outer(p$b_sin, sin(2 * pi * time / 24)) +
    outer(p$b_cos, cos(2 * pi * time / 24)))
set.seed(123)
beta_independent <- matrix(rbeta(length(mu), mu * p$phi,
    (1 - mu) * p$phi), nrow = nrow(p),
    dimnames = list(p$cpg, paste0("sample", seq_along(time))))
se_independent <- SummarizedExperiment(
    assays = list(beta = beta_independent),
    rowData = rowData(se)[seq_len(1200), ],
    colData = DataFrame(time = time))
metadata(se_independent)$simulation <- list(
    source = "First 1,200 parameter rows; independent observations without subject effects",
    observation_seed = 123L, n_samples = 48L, period = 24)

saveRDS(list(se = se, se_independent = se_independent,
    sample_meta = as.data.frame(colData(se))),
    output, compress = "xz", version = 3)
