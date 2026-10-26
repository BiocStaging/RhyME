.simulate_auto_prior_data <- function() {
    withr::local_seed(123)
    n_cpg <- 1200L
    rhythmic <- seq_len(n_cpg) <= 120L
    time <- rep(seq(0, 22, by = 2), 4)
    phase <- rep(c(0, 0.25, 6, 12, 18, 23.75), length.out = n_cpg)
    baseline <- rep(qlogis(c(0.2, 0.5, 0.8)), length.out = n_cpg)
    amplitude <- ifelse(rhythmic, 0.6, 0)
    phi <- exp(rnorm(n_cpg, log(50), 0.4))
    # Generate from peak time directly, independently of the fitted coefficients.
    eta <- baseline + amplitude * cos(2 * pi * outer(phase, time, "-") / 24)
    mu <- plogis(eta)
    beta <- matrix(rbeta(length(mu), mu * phi, (1 - mu) * phi),
        nrow = n_cpg,
        dimnames = list(paste0("cg", seq_len(n_cpg)),
            paste0("sample", seq_along(time))))
    SummarizedExperiment(
        assays = list(beta = beta),
        rowData = S4Vectors::DataFrame(rhythmic = rhythmic,
            true_phase = ifelse(rhythmic, phase, NA_real_)),
        colData = S4Vectors::DataFrame(time = time))
}

test_auto_prior_fixed_effects <- function() {
    se <- .simulate_auto_prior_data()
    # Omitting prior and subject_col exercises the public defaults.
    fit <- fitCosinor(se, BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE)
    result <- rowData(fit)
    prior <- S4Vectors::metadata(fit)$RhyME$ebPrior
    checkEquals(S4Vectors::metadata(fit)$RhyME$design$mode, "fixed")
    checkEquals(prior$n, nrow(se))
    checkTrue(is.finite(prior$mean))
    checkTrue(is.finite(prior$sd) && prior$sd > 0)
    checkTrue(all(result$RhyME_converged))
    checkTrue(all(is.finite(result$RhyME_pvalue)))
    checkTrue(all(result$RhyME_pvalue >= 0 & result$RhyME_pvalue <= 1))
    checkEquals(result$RhyME_padj, p.adjust(result$RhyME_pvalue, "BH"))

    # Strong rhythmic signals should separate from the simulated nulls.
    rhythmic <- result$rhythmic
    top <- order(result$RhyME_pvalue)[seq_len(sum(rhythmic))]
    checkTrue(mean(rhythmic[top]) > 0.9)
    checkTrue(mean(result$RhyME_padj[rhythmic] < 0.05) > 0.9)
    checkTrue(median(result$RhyME_pvalue[!rhythmic]) > 0.25)

    phase <- result$RhyME_phase[rhythmic]
    checkTrue(all(is.finite(phase) & phase >= 0 & phase < 24))
    # Circular distance includes signals whose peaks are close to 0/24 hours.
    error <- abs((phase - result$true_phase[rhythmic] + 12) %% 24 - 12)
    checkTrue(median(error) < 0.5)
    checkTrue(mean(error < 1.5) > 0.95)
}

test_auto_prior_repeated_measures <- function() {
    example <- readRDS(system.file("extdata", "example_data.rds",
        package = "RhyME"))
    fit <- fitCosinor(example$se, subject_col = "subject", prior = "auto",
        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    info <- S4Vectors::metadata(fit)$RhyME
    checkEquals(info$design$mode, "mixed")
    checkEquals(info$design$n_subjects,
        length(unique(example$se$subject)))
    checkTrue(info$ebPrior$n >= 1000L && info$ebPrior$n <= nrow(fit))
    checkTrue(is.finite(info$ebPrior$mean))
    checkTrue(is.finite(info$ebPrior$sd) && info$ebPrior$sd > 0)
    checkTrue(all(is.finite(rowData(fit)$RhyME_pvalue)))
    checkTrue(all(rowData(fit)$RhyME_converged))
}

test_auto_prior_requires_usable_cpgs <- function() {
    se <- .simulate_auto_prior_data()
    error <- tryCatch(fitCosinor(se[seq_len(999), ], prior = "auto",
        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE),
        error = conditionMessage)
    checkTrue(is.character(error))
    checkTrue(grepl("at least 1,000 CpGs", error, fixed = TRUE))

    # Enough input rows, but fewer than 1,000 usable pilot fits.
    assay(se)[seq_len(201), ] <- NA_real_
    error <- tryCatch(fitCosinor(se, prior = "auto",
        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE),
        error = conditionMessage)
    checkTrue(is.character(error))
    checkTrue(grepl("at least 1,000 usable null pilot fits", error,
        fixed = TRUE))
}
