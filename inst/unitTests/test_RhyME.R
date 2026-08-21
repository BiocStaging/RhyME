make_test_se <- function(n_cpg = 50, n_subjects = 5, n_timepoints = 6,
                         n_rhythmic = 10, seed = 42) {
    set.seed(seed)

    times <- seq(0, 20, length.out = n_timepoints)
    n_samples <- n_subjects * n_timepoints
    omega <- 2 * pi / 24

    meta <- data.frame(
        time = rep(times, n_subjects),
        subject = rep(paste0("S", seq_len(n_subjects)), each = n_timepoints),
        stringsAsFactors = FALSE
    )
    meta$sample_id <- paste0(meta$subject, "_T", meta$time)

    sin_t <- sin(omega * meta$time)
    cos_t <- cos(omega * meta$time)

    beta_mat <- matrix(NA_real_, nrow = n_cpg, ncol = n_samples)
    rownames(beta_mat) <- paste0("cg", sprintf("%04d", seq_len(n_cpg)))
    colnames(beta_mat) <- meta$sample_id

    for (i in seq_len(n_cpg)) {
        b0 <- rnorm(1, mean = 0, sd = 0.5)
        if (i <= n_rhythmic) {
            b_sin <- rnorm(1, mean = 0.3, sd = 0.1)
            b_cos <- rnorm(1, mean = 0.3, sd = 0.1)
        } else {
            b_sin <- 0
            b_cos <- 0
        }
        tau <- 0.2
        phi <- exp(rnorm(1, mean = 3, sd = 0.3))

        u <- rnorm(n_subjects, 0, tau)
        eta <- b0 + b_sin * sin_t + b_cos * cos_t +
            rep(u, each = n_timepoints)
        mu <- plogis(eta)

        y <- rbeta(n_samples, mu * phi, (1 - mu) * phi)
        beta_mat[i, ] <- pmin(pmax(y, 1e-4), 1 - 1e-4)
    }

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(beta = beta_mat),
        colData = S4Vectors::DataFrame(meta)
    )
    list(se = se)
}

make_test_se_fixed <- function(n_cpg = 30, n_samples = 20,
                               n_rhythmic = 5, seed = 123) {
    set.seed(seed)

    times <- runif(n_samples, 0, 24)
    omega <- 2 * pi / 24

    meta <- data.frame(
        time = times,
        sample_id = paste0("S", seq_len(n_samples)),
        stringsAsFactors = FALSE
    )

    sin_t <- sin(omega * times)
    cos_t <- cos(omega * times)

    beta_mat <- matrix(NA_real_, nrow = n_cpg, ncol = n_samples)
    rownames(beta_mat) <- paste0("cg", sprintf("%04d", seq_len(n_cpg)))
    colnames(beta_mat) <- meta$sample_id

    for (i in seq_len(n_cpg)) {
        b0 <- rnorm(1, mean = 0, sd = 0.5)
        if (i <= n_rhythmic) {
            b_sin <- rnorm(1, mean = 0.4, sd = 0.1)
            b_cos <- rnorm(1, mean = 0.4, sd = 0.1)
        } else {
            b_sin <- 0
            b_cos <- 0
        }
        phi <- exp(rnorm(1, mean = 3, sd = 0.3))

        eta <- b0 + b_sin * sin_t + b_cos * cos_t
        mu <- plogis(eta)

        y <- rbeta(n_samples, mu * phi, (1 - mu) * phi)
        beta_mat[i, ] <- pmin(pmax(y, 1e-4), 1 - 1e-4)
    }

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(beta = beta_mat),
        colData = S4Vectors::DataFrame(meta)
    )
    list(se = se)
}

test_RhyME_mixed_workflow <- function() {
    dat <- make_test_se(n_cpg = 60, n_subjects = 5, n_timepoints = 6,
                        n_rhythmic = 10, seed = 42)

    set.seed(1)
    se <- fitCosinor(
        dat$se,
        time_col = "time",
        subject_col = "subject",
        prior = "none",
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )

    res <- SummarizedExperiment::rowData(se)
    meta <- S4Vectors::metadata(se)$RhyME

    checkTrue(inherits(se, "SummarizedExperiment"))
    checkEquals(nrow(res), 60)
    checkEquals(meta$design$mode, "mixed")
    checkEquals(meta$prior$method, "none")
    checkEquals(meta$test$method, "Wald")
    checkEquals(meta$test$p.adjust.method, "BH")
    checkEquals(length(meta$ebPrior), 0L)
    checkEquals(meta$fitSummary$total, 60)
    checkEquals(meta$fitSummary$converged, sum(res$converged %in% TRUE))
    checkEquals(meta$fitSummary$failed, 60 - sum(res$converged %in% TRUE))
    checkTrue(!"relative_amplitude" %in% names(res))
    hidden_cols <- c(
        "phi", "log_phi", "log_phi_pilot", "tau", "error",
        "amplitude_link", "mesor_link", "statistic", "df",
        "se_b0", "se_sin", "se_cos", "cov_sin_cos"
    )
    checkTrue(!any(hidden_cols %in% names(res)))

    top <- topCpGs(se, n = 5)
    checkTrue(inherits(top, "DataFrame"))
    checkTrue(nrow(top) <= 5)
    checkEquals(names(top), c("cpg", "amplitude_response", "phase",
                              "mesor_response", "pvalue", "padj"))
    if (nrow(top) > 1) {
        checkTrue(all(diff(top$pvalue) >= 0, na.rm = TRUE))
    }

    p <- plotCosinorFit(se, cpgs = top$cpg[1], n_grid = 20)
    checkTrue(inherits(p, "ggplot"))
}

test_RhyME_fixed_tests <- function() {
    dat <- make_test_se_fixed(n_cpg = 25, n_samples = 20,
                              n_rhythmic = 5, seed = 123)

    for (method in c("Wald", "LRT")) {
        se <- fitCosinor(
            dat$se,
            time_col = "time",
            prior = "none",
            test = method,
            BPPARAM = BiocParallel::SerialParam(),
            verbose = FALSE
        )

        res <- SummarizedExperiment::rowData(se)
        meta <- S4Vectors::metadata(se)$RhyME
        checkEquals(meta$design$mode, "fixed")
        checkEquals(meta$test$method, method)
        checkEquals(meta$test$p.adjust.method, "BH")
        checkTrue(all(res$pvalue >= 0 & res$pvalue <= 1, na.rm = TRUE))
        checkTrue(!any(c(
            "statistic", "df", "se_b0", "se_sin", "se_cos",
            "cov_sin_cos"
        ) %in% names(res)))
    }
}

test_RhyME_covariate_adjustment <- function() {
    dat <- make_test_se(n_cpg = 12, n_subjects = 8, n_timepoints = 6,
                        n_rhythmic = 4, seed = 314)
    cd <- as.data.frame(SummarizedExperiment::colData(dat$se))
    cd$sex <- factor(rep(rep(c("female", "male"), 4), each = 6))
    cd$age <- rep(seq(30, 65, length.out = 8), each = 6)
    SummarizedExperiment::colData(dat$se) <- S4Vectors::DataFrame(cd)

    fit <- fitCosinor(
        dat$se,
        time_col = "time",
        subject_col = "subject",
        covariates = c("sex", "age"),
        prior = "none",
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )

    res <- SummarizedExperiment::rowData(fit)
    meta <- S4Vectors::metadata(fit)$RhyME
    checkEquals(meta$design$covariates, c("sex", "age"))
    checkTrue(any(is.finite(res$pvalue)))
    checkTrue(!any(grepl("sex|age", names(res))))

    X_cov <- RhyME:::.makeCovariateMatrix(cd, c("sex", "age"))
    terms <- RhyME:::.makeCosinorTerms(cd$time)
    design <- RhyME:::.makeDesignMatrices(terms, X_cov)
    checkEquals(ncol(design$full), ncol(design$null) + 2L)
    checkEquals(
        unname(design$full[, seq_len(ncol(design$null)), drop = FALSE]),
        unname(design$null)
    )
}

test_RhyME_mixed_wald_joint <- function() {
    dat <- make_test_se(n_cpg = 12, n_subjects = 5, n_timepoints = 6,
                        n_rhythmic = 4, seed = 321)

    se <- fitCosinor(
        dat$se,
        time_col = "time",
        subject_col = "subject",
        prior = "none",
        test = "Wald",
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )
    res <- SummarizedExperiment::rowData(se)
    usable <- which(is.finite(res$pvalue))

    checkEquals(S4Vectors::metadata(se)$RhyME$design$mode, "mixed")
    checkTrue(length(usable) > 0)
    checkTrue(all(res$pvalue[usable] >= 0 & res$pvalue[usable] <= 1))
}

test_RhyME_default_assay <- function() {
    dat <- make_test_se_fixed(n_cpg = 20, n_samples = 18,
                              n_rhythmic = 4, seed = 99)
    SummarizedExperiment::assayNames(dat$se) <- "methylation"

    se <- fitCosinor(
        dat$se,
        time_col = "time",
        prior = "none",
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )

    meta <- S4Vectors::metadata(se)$RhyME
    checkEquals(meta$design$assay_name, "methylation")
    checkEquals(meta$design$mode, "fixed")
}

test_RhyME_mode_inference <- function() {
    dat <- make_test_se(n_cpg = 12, n_subjects = 4, n_timepoints = 4,
                        n_rhythmic = 3, seed = 7)

    fixed <- fitCosinor(
        dat$se,
        time_col = "time",
        prior = "none",
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )
    checkEquals(S4Vectors::metadata(fixed)$RhyME$design$mode, "fixed")
    checkTrue(is.null(S4Vectors::metadata(fixed)$RhyME$design$subject_col))

    mixed <- fitCosinor(
        dat$se,
        time_col = "time",
        subject_col = "subject",
        prior = "none",
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )
    checkEquals(S4Vectors::metadata(mixed)$RhyME$design$mode, "mixed")
    checkEquals(S4Vectors::metadata(mixed)$RhyME$design$subject_col,
                "subject")
}

test_RhyME_internal_summaries <- function() {
    terms <- RhyME:::.makeCosinorTerms(c(0, 6, 12), period = 24)
    checkEquals(terms$sin_term, c(0, 1, 0), tolerance = 1e-10)
    checkEquals(terms$cos_term, c(1, 0, -1), tolerance = 1e-10)

    params <- RhyME:::.computeCircadianParams(
        b0 = 0, b_sin = 0.3, b_cos = 0.4, period = 24
    )
    checkEquals(names(params), c("amplitude_link", "amplitude_response",
                                 "phase", "mesor_link", "mesor_response"))
    checkEquals(params$amplitude_link, 0.5, tolerance = 1e-10)

    # With eta(t) = sin(2*pi*t/24), the fitted peak is at 6 h.
    sine_peak <- RhyME:::.computeCircadianParams(
        b0 = 0, b_sin = 1, b_cos = 0, period = 24
    )
    checkEquals(sine_peak$phase, 6, tolerance = 1e-10)

    # The numerical period average must not count both 0 h and 24 h.
    grid_params <- RhyME:::.computeCircadianParams(
        b0 = -0.4, b_sin = 0.7, b_cos = 0.2,
        period = 24, grid_by = 1
    )
    tgrid <- 0:23
    expected_mesor <- mean(plogis(
        -0.4 + 0.7 * sin(2 * pi * tgrid / 24) +
            0.2 * cos(2 * pi * tgrid / 24)
    ))
    checkEquals(grid_params$mesor_response, expected_mesor,
                tolerance = 1e-12)

    marginal <- RhyME:::.computeCircadianParams(
        b0 = -1, b_sin = 0.5, b_cos = 0.3, period = 24, tau = 1
    )
    checkTrue(!isTRUE(all.equal(params$mesor_response,
                                marginal$mesor_response)))

    set.seed(42)
    log_phi <- rnorm(2000, mean = 3.0, sd = 0.5)
    prior <- RhyME:::.estimateLogPhiPrior(log_phi)
    checkEquals(names(prior), c("mean", "sd"))
    checkEquals(prior$mean, 3.0, tolerance = 0.15)
    checkEquals(prior$sd, 0.5, tolerance = 0.15)

    contaminated <- c(log_phi, -100, 100)
    robust_prior <- RhyME:::.estimateLogPhiPrior(contaminated)
    checkEquals(robust_prior$mean, prior$mean, tolerance = 0.02)
    checkEquals(robust_prior$sd, prior$sd, tolerance = 0.02)

    checkException(RhyME:::.estimateLogPhiPrior(rep(NA_real_, 10)))
}
