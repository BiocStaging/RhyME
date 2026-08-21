#' Fit RhyME across CpG sites
#'
#' Fits a RhyME to each CpG site in a methylation array
#' dataset. Supports both fixed-effects (cross-sectional) and mixed-effects
#' (repeated measures with random intercepts) designs, and optionally
#' applies empirical Bayes shrinkage to the precision parameter (phi).
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] object. The
#'   specified assay should contain methylation beta values in [0, 1]. Values
#'   exactly at the boundaries are moved slightly inward for beta regression.
#' @param time_col Character, column name in `colData(se)` containing
#'   collection times in hours (default `"time"`).
#' @param subject_col Character or NULL. Column name in `colData(se)` for
#'   subject identifiers. If non-NULL, random intercepts per subject are
#'   included (mixed-effects mode). If NULL, fixed-effects mode is used.
#' @param covariates Character vector or NULL. Optional `colData(se)` columns
#'   included as adjustment covariates in both the rhythmic and null models.
#'   Covariate coefficients are treated as nuisance parameters and are not
#'   tested or returned.
#' @param assay_name Character or NULL, name of the assay to use. The default
#'   NULL uses the first assay in `se`.
#' @param period Numeric, the oscillation period in hours (default 24 for
#'   circadian rhythms). Can be set to other values for ultradian or
#'   infradian rhythms.
#' @param prior Character, either `"auto"` (default) to estimate a robust,
#'   bootstrap-corrected empirical-Bayes prior for log(phi), or `"none"` for
#'   ordinary maximum likelihood.
#' @param test Rhythmicity test: a covariance-aware joint Wald test of the sine
#'   and cosine coefficients (`"Wald"`, default), or a likelihood-ratio test
#'   (`"LRT"`).
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object for execution.
#'   Defaults to [BiocParallel::MulticoreParam()] for parallel fitting on
#'   Unix-like systems. Pass [BiocParallel::SerialParam()] for serial fitting or
#'   [BiocParallel::SnowParam()] for socket-based parallelism.
#' @param verbose Logical, print messages and show supported progress bars.
#'
#' @return The input `SummarizedExperiment` with per-CpG results appended to
#'   `rowData(se)` and design info, EB prior, test settings, and fit summary
#'   stored in `metadata(se)$RhyME`.
#'
#' @details When `covariates` are supplied, their centered design columns are
#' included in both the full and null models. Consequently, the reported
#' p-value still tests only the joint sine/cosine rhythmicity effect, while
#' response-scale summaries describe the mean covariate profile.
#'
#' @export
#'
#' @examples
#' ex <- readRDS(system.file("extdata", "example_data.rds",
#'                          package = "RhyME"))
#' se <- ex$se[seq_len(5), ]
#'
#' fit <- fitCosinor(se, time_col = "time", subject_col = "subject",
#'                   prior = "none", verbose = FALSE)
#' head(SummarizedExperiment::rowData(fit))
fitCosinor <- function(se,
                       time_col = "time",
                       subject_col = NULL,
                       assay_name = NULL,
                       period = 24,
                       prior = c("auto", "none"),
                       test = c("Wald", "LRT"),
                       BPPARAM = BiocParallel::MulticoreParam(),
                       verbose = TRUE,
                       covariates = NULL) {

    prior <- match.arg(prior)
    test_method <- match.arg(test)
    if (!is.null(subject_col) && subject_col %in% covariates) {
        stop("'subject_col' cannot also be an adjustment covariate.",
             call. = FALSE)
    }
    if (!is.numeric(period) || length(period) != 1L ||
        !is.finite(period) || period <= 0) {
        stop("'period' must be a positive number.", call. = FALSE)
    }

    squeeze_bounds <- .RHYME_SQUEEZE_BOUNDS
    # Keep the pilot and bootstrap passes quiet. A single progress bar is
    # enabled immediately before the final genome-wide fit.
    BPPARAM <- .setBPProgressBar(BPPARAM, FALSE)

    # A subject column is the single switch between fixed and mixed modes.
    mode <- if (is.null(subject_col)) "fixed" else "mixed"

    input <- .getRhyMEInput(
        se, assay_name, time_col, subject_col, covariates
    )
    beta_mat <- input$beta
    cd <- input$col_data
    assay_name <- input$assay_name
    cpg_ids <- rownames(beta_mat)
    n_cpg <- nrow(beta_mat)

    # Build per-sample data frame with cosinor terms
    cos_terms <- .makeCosinorTerms(cd[[time_col]], period = period)
    sample_df <- data.frame(
        sin_term = cos_terms$sin_term,
        cos_term = cos_terms$cos_term
    )
    if (mode == "mixed") {
        sample_df$subject <- cd[[subject_col]]
    }
    X_cov <- .makeCovariateMatrix(cd, covariates)

    prior_spec <- NULL

    if (prior == "auto") {
        # ---- Pass 1: Pilot fit on a random subset (no prior) ----
        n_pilot <- min(.RHYME_PILOT_N, n_cpg)
        if (n_pilot < .RHYME_MIN_PILOT_CPGS) {
            stop(sprintf(
                paste0("EB prior needs at least %d pilot CpGs, but only %d ",
                       "are available. Supply more CpGs or set ",
                       "prior = 'none'."),
                .RHYME_MIN_PILOT_CPGS, n_pilot
            ), call. = FALSE)
        }
        pilot_idx <- if (n_pilot < n_cpg) {
            sample(seq_len(n_cpg), n_pilot)
        } else {
            seq_len(n_cpg)
        }

        pilot_results <- BiocParallel::bplapply(
            pilot_idx,
            function(i) {
                .fitOneCpG(
                    cpg_id = cpg_ids[i],
                    y = as.numeric(beta_mat[i, ]),
                    df = sample_df,
                    X_cov = X_cov,
                    mode = mode,
                    prior_phi = NULL,
                    eps = squeeze_bounds[1],
                    test = "none"
                )
            },
            BPPARAM = BPPARAM
        )

        pilot_log_phi_subset <- vapply(pilot_results, function(r) {
            if (isTRUE(r$converged)) r$log_phi else NA_real_
        }, numeric(1))

        raw_prior <- .estimateLogPhiPrior(pilot_log_phi_subset)
        prior_spec <- .bootstrapPhiPrior(
            pilot_results = pilot_results,
            beta_pilot = beta_mat[pilot_idx, , drop = FALSE],
            sample_df = sample_df,
            X_cov = X_cov,
            mode = mode,
            raw_prior = raw_prior,
            reps = .RHYME_BOOTSTRAP_REPS,
            variance_floor = .RHYME_BOOTSTRAP_VARIANCE_FLOOR,
            squeeze_bounds = squeeze_bounds,
            BPPARAM = BPPARAM
        )
        if (verbose) {
            message(sprintf(
                "RhyME prior: log(phi) mean=%.3f, sd=%.3f (bootstrap B=%d)",
                prior_spec$mean, prior_spec$sd,
                .RHYME_BOOTSTRAP_REPS
            ))
        }
    }

    BPPARAM <- .setBPProgressBar(BPPARAM, verbose)
    final_results <- BiocParallel::bplapply(
        seq_len(n_cpg),
        function(i) {
            .fitOneCpG(
                cpg_id = cpg_ids[i],
                y = as.numeric(beta_mat[i, ]),
                df = sample_df,
                X_cov = X_cov,
                mode = mode,
                prior_phi = prior_spec,
                eps = squeeze_bounds[1],
                test = test_method
            )
        },
        BPPARAM = BPPARAM
    )

    se <- .assembleResults(se, final_results, period)

    converged <- rowData(se)$converged %in% TRUE
    metadata(se)$RhyME <- list(
        design = list(
            period = period,
            mode = mode,
            assay_name = assay_name,
            time_col = time_col,
            subject_col = subject_col,
            covariates = if (is.null(covariates)) character() else covariates,
            n_samples = ncol(beta_mat),
            n_subjects = if (mode == "mixed") {
                length(unique(cd[[subject_col]]))
            } else {
                NA_integer_
            }
        ),
        prior = list(method = prior),
        ebPrior = prior_spec,
        test = list(method = test_method, p.adjust.method = "BH"),
        fitSummary = list(
            total = n_cpg,
            converged = sum(converged),
            failed = n_cpg - sum(converged)
        )
    )

    if (verbose) .messageFitSummary(se)
    se
}

#' Assemble per-CpG results into the SummarizedExperiment
#' @noRd
.assembleResults <- function(se, results, period) {

    num <- function(field) vapply(results, `[[`, numeric(1), field)
    chr <- function(field) vapply(results, `[[`, character(1), field)
    flag <- function(field) vapply(results, `[[`, logical(1), field)

    b0 <- num("b0")
    b_sin <- num("b_sin")
    b_cos <- num("b_cos")
    tau <- num("tau")

    # Circadian parameters (computed from coefficients; response-scale
    # summaries are marginal over random intercepts when tau is available).
    circ <- .computeCircadianParams(b0, b_sin, b_cos, period, tau = tau)

    pvalue <- num("pvalue")
    res_df <- DataFrame(
        b0 = b0,
        b_sin = b_sin,
        b_cos = b_cos,
        amplitude_link = circ$amplitude_link,
        amplitude_response = circ$amplitude_response,
        phase = circ$phase,
        mesor_link = circ$mesor_link,
        mesor_response = circ$mesor_response,
        tau = tau,
        pvalue = pvalue,
        padj = stats::p.adjust(pvalue, method = "BH"),
        converged = flag("converged"),
        row.names = chr("cpg")
    )

    rd <- rowData(se)
    keep_cols <- setdiff(names(rd), c("RhyME", names(res_df)))
    rowData(se) <- cbind(rd[, keep_cols, drop = FALSE], res_df)
    se
}

#' Enable BiocParallel progress bars when supported
#' @noRd
.setBPProgressBar <- function(BPPARAM, progressbar) {
    tryCatch({
        BiocParallel::bpprogressbar(BPPARAM) <- isTRUE(progressbar)
        BPPARAM
    }, error = function(e) BPPARAM)
}

#' Print a compact fit summary
#' @noRd
.messageFitSummary <- function(se) {
    fs <- metadata(se)$RhyME$fitSummary
    message(sprintf(
        "RhyME fit: %d/%d CpGs converged (%.1f%%); %d failed.",
        fs$converged, fs$total, 100 * fs$converged / fs$total, fs$failed
    ))
    invisible(NULL)
}
