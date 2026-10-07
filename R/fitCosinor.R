#' Fit Cosinor Models Across CpG Sites
#'
#' Fits beta cosinor models using PQL, with optional random intercepts,
#' adjustment covariates, and an empirical Bayes (EB) prior on log precision.
#' A null-model score test assesses rhythmicity; a full model estimates
#' rhythm parameters.
#'
#' @param se A `SummarizedExperiment` with CpGs in rows and samples in
#' columns. The assay contains numeric beta values in \eqn{[0,1]}, with
#' missing values allowed.
#' @param time_col Name of the collection-time column in `colData(se)`,
#' in hours.
#' @param subject_col Name of the subject column for a random-intercept
#' model, or `NULL` for a fixed-effects model of independent samples.
#' @param assay_name Assay name, or `NULL` to use the first assay.
#' @param period Oscillation period in hours. Default is 24.
#' @param prior `"auto"` estimates an EB prior; `"none"` uses no prior.
#' A named list with `mean` and positive `sd` supplies a normal prior on
#' log precision. EB estimation requires at least 1,000 usable pilot fits.
#' @param BPPARAM A `BiocParallelParam` object. Default is
#' `BiocParallel::MulticoreParam()` (serial on Windows).
#' @param verbose Whether to print fitting progress.
#' @param covariates Names of numeric or categorical `colData(se)` columns
#' for adjustment. Covariate coefficients are not returned.
#'
#' @return The input `SummarizedExperiment` with these `rowData(se)` columns:
#' - `b0`, `b_sin`, `b_cos`: logit-scale intercept and sine/cosine coefficients.
#' - `log_phi`, `tau`: log precision and random-intercept SD.
#' - `score_chisq`, `pvalue`, `padj`: score statistic, p-value, and
#'   BH-adjusted p-value across all input CpGs.
#' - `logit_amplitude`, `amplitude`: logit-scale amplitude and half the
#'   marginal mean curve's peak-to-trough difference on the methylation scale.
#' - `phase`, `mesor`: peak time in hours in \eqn{[0,\mathrm{period})} and
#'   the marginal mean curve averaged over one period.
#' - `converged`: whether both null and full fits converged.
#'
#' Unusable null fits have `NA` p-values; unusable full fits have `NA`
#' parameter estimates. Usable results from either fit are retained.
#' Phase is undefined at zero amplitude. Existing row annotations are retained.
#' `metadata(se)$RhyME` stores `design`, `ebPrior`, and `fitSummary`.
#'
#' @seealso [plotCosinorFit()]
#' @export
#'
#' @examples
#' ex <- readRDS(system.file("extdata", "example_data.rds", package = "RhyME"))
#' fit <- fitCosinor(ex$se[seq_len(5), ], subject_col = "subject",
#'     prior = "none", verbose = FALSE)
#' head(SummarizedExperiment::rowData(fit), 3)
fitCosinor <- function(se, time_col = "time", subject_col = NULL,
                       assay_name = NULL, period = 24,
                       prior = c("auto", "none"),
                       BPPARAM = BiocParallel::MulticoreParam(), verbose = TRUE,
                       covariates = NULL) {
    if (!is.null(subject_col) && subject_col %in% covariates)
        stop("'subject_col' cannot also be an adjustment covariate.", call. = FALSE)
    if (!is.numeric(period) || length(period) != 1L ||
        !is.finite(period) || period <= 0)
        stop("'period' must be a positive finite number.", call. = FALSE)
    if (is.list(prior)) {
        .validatePQLPrior(prior)
        prior_spec <- prior
        prior_method <- "supplied"
    } else {
        prior_method <- match.arg(prior)
        prior_spec <- NULL
    }
    input <- .getRhyMEInput(se, assay_name, time_col, subject_col, covariates)
    beta <- input$beta
    cd <- input$col_data
    # Center once, before omitting missing beta values within each CpG.
    X_cov <- .makeCovariateMatrix(cd, covariates)
    X <- cbind(`(Intercept)` = rep(1, nrow(cd)), X_cov,
        sin_term = sin(2 * pi * cd[[time_col]] / period),
        cos_term = cos(2 * pi * cd[[time_col]] / period)
    )
    subject <- if (is.null(subject_col)) NULL else as.integer(factor(cd[[subject_col]])) - 1L
    .checkPQLDesign(X, subject)
    if (prior_method == "auto") {
        if (nrow(beta) < 1000L)
            stop("EB prior estimation requires at least 1,000 CpGs; use prior = 'none' for smaller datasets.",
                 call. = FALSE)
        pilot_index <- .pqlPilotIndex(nrow(beta))
        if (verbose) message("Estimating EB prior on log precision from null pilot fits...")
        pilot <- BiocParallel::bplapply(
            .pqlBatches(beta[pilot_index, , drop = FALSE], BiocParallel::bpworkers(BPPARAM)),
            .fitPQLPilotBatch, design = X, subject = subject, BPPARAM = BPPARAM)
        prior_spec <- .estimateLogPhiPrior(unlist(pilot, use.names = FALSE))
    }
    if (verbose) message("Fitting null score tests and full-model parameters...")
    final_results <- BiocParallel::bplapply(
        .pqlBatches(beta, BiocParallel::bpworkers(BPPARAM)), .fitPQLBatch,
        design = X, subject = subject, prior_phi = prior_spec, BPPARAM = BPPARAM)
    if (verbose) message("Combining batches and computing rhythm summaries and BH-adjusted p-values...")
    se <- .assembleResults(se, final_results, period)
    fit_summary <- list(total = nrow(se),
                        usable_tests = sum(!is.na(rowData(se)$pvalue)),
                        usable_parameters = sum(!is.na(rowData(se)$b0)))
    metadata(se)$RhyME <- list(
        design = list(period = period,
                      mode = if (is.null(subject_col)) "fixed" else "mixed",
                      assay_name = input$assay_name,
                      time_col = time_col, subject_col = subject_col,
                      covariates = if (is.null(covariates)) character() else covariates,
                      n_samples = ncol(beta),
                      n_subjects = if (is.null(subject_col)) NA_integer_ else
                          length(unique(subject))),
        ebPrior = prior_spec,
        fitSummary = fit_summary
    )
    if (verbose) message(sprintf(
        "RhyME PQL: %d/%d usable score tests; %d/%d usable parameter fits.",
        fit_summary$usable_tests, nrow(se), fit_summary$usable_parameters, nrow(se)))
    se
}

.assembleResults <- function(se, results, period) {
    res <- do.call(rbind, results)
    circ <- .computeCircadianParams(res$b0, res$b_sin, res$b_cos,
                                    period, tau = res$tau)
    parameters <- c("b0", "b_sin", "b_cos", "log_phi", "tau", "score_chisq", "pvalue")
    result <- DataFrame(res[, parameters],
        padj = stats::p.adjust(ifelse(is.na(res$pvalue), 1, res$pvalue), method = "BH"),
        as.data.frame(circ),
        converged = res$converged)
    old <- rowData(se)
    keep <- setdiff(names(old), c("RhyME", "converged", names(res), names(result),
        "status_null", "status_full", "error_null", "error_full",
        "phi", "amplitude_link", "amplitude_response", "mesor_link", "mesor_response"))
    rowData(se) <- cbind(result, old[, keep, drop = FALSE])
    se
}
