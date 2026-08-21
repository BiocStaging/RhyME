# Internal utility functions

.RHYME_MIN_PILOT_CPGS <- 1000L
.RHYME_PILOT_N <- 10000L
.RHYME_BOOTSTRAP_REPS <- 5L
.RHYME_BOOTSTRAP_VARIANCE_FLOOR <- 0.01
.RHYME_SQUEEZE_BOUNDS <- c(1e-4, 1 - 1e-4)

.invlogit <- function(x) {
    1 / (1 + exp(-x))
}

#' Build cosinor terms from a time vector
#'
#' @param time Numeric vector of collection times (in hours).
#' @param period Numeric scalar, the period in hours (default 24).
#' @return A data.frame with columns `sin_term` and `cos_term`.
#' @noRd
.makeCosinorTerms <- function(time, period = 24) {
    omega <- 2 * pi / period
    data.frame(
        sin_term = sin(omega * time),
        cos_term = cos(omega * time)
    )
}

#' Construct a fixed-effect matrix for adjustment covariates
#' @noRd
.makeCovariateMatrix <- function(col_data, covariates = NULL) {
    if (is.null(covariates) || !length(covariates)) {
        return(matrix(numeric(), nrow = nrow(col_data), ncol = 0L))
    }
    if (!is.character(covariates) || anyNA(covariates) ||
        any(!nzchar(covariates)) || anyDuplicated(covariates)) {
        stop("'covariates' must contain unique colData column names.",
             call. = FALSE)
    }

    cov_data <- col_data[, covariates, drop = FALSE]
    if (any(!stats::complete.cases(cov_data))) {
        stop("Adjustment covariates cannot contain missing values.",
             call. = FALSE)
    }

    X <- stats::model.matrix(~ ., data = cov_data)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    X <- sweep(X, 2L, colMeans(X), FUN = "-")
    colnames(X) <- paste0("cov_", make.names(colnames(X), unique = TRUE))
    X
}

#' Construct full and null model matrices
#' @noRd
.makeDesignMatrices <- function(sample_df, X_cov = NULL) {
    if (is.null(X_cov)) {
        X_cov <- matrix(numeric(), nrow = nrow(sample_df), ncol = 0L)
    }
    base <- cbind(
        `(Intercept)` = rep(1, nrow(sample_df)),
        X_cov
    )
    list(
        full = cbind(
            base,
            sin_term = sample_df$sin_term,
            cos_term = sample_df$cos_term
        ),
        null = base
    )
}

#' Estimate an empirical-Bayes prior for beta precision
#'
#' Estimates a Normal prior on log(phi) from per-CpG pilot estimates. The
#' center is the median and the scale is IQR / 1.349, making the estimate less
#' sensitive to unstable per-CpG precision fits. The prior is used in the
#' second-pass TMB fit to shrink CpG-specific precision estimates toward the
#' genome-wide distribution.
#'
#' @param log_phi_hat Numeric vector of per-CpG log(phi) estimates from the
#'   pilot fit.
#' @return A list containing `mean` and `sd`.
#' @noRd
.estimateLogPhiPrior <- function(log_phi_hat) {
    valid <- log_phi_hat[is.finite(log_phi_hat)]
    n_ok <- length(valid)

    if (n_ok < .RHYME_MIN_PILOT_CPGS) {
        stop(sprintf(
            paste0("Automatic prior requires at least %d finite log(phi) ",
                   "estimates; got %d."),
            .RHYME_MIN_PILOT_CPGS, n_ok
        ), call. = FALSE)
    }

    center <- stats::median(valid)
    robust_sd <- stats::IQR(valid) / 1.349

    list(
        mean = center,
        sd = robust_sd
    )
}

#' Parametric-bootstrap correction for the log(phi) prior

#' @noRd
.bootstrapPhiPrior <- function(pilot_results, beta_pilot, sample_df, X_cov,
                               mode,
                               raw_prior, reps, variance_floor,
                               squeeze_bounds, BPPARAM) {
    X_full <- .makeDesignMatrices(sample_df, X_cov)$full
    usable <- vapply(pilot_results, function(fit) {
        parameter_ok <- length(fit$beta) == ncol(X_full) &&
            all(is.finite(fit$beta)) && is.finite(fit$log_phi)
        tau_ok <- mode != "mixed" ||
            (length(fit$tau) == 1L && is.finite(fit$tau) && fit$tau >= 0)
        isTRUE(fit$converged) && parameter_ok && tau_ok
    }, logical(1))

    if (sum(usable) < .RHYME_MIN_PILOT_CPGS) {
        stop(sprintf(
            paste0("Automatic-prior bootstrap has only %d usable pilot ",
                   "fits; at least %d are required."),
            sum(usable), .RHYME_MIN_PILOT_CPGS
        ), call. = FALSE)
    }

    fits <- pilot_results[usable]
    beta_pilot <- as.matrix(beta_pilot[usable, , drop = FALSE])
    n_fit <- length(fits)
    n_sample <- ncol(beta_pilot)
    common_phi <- exp(raw_prior$mean)

    if (mode == "mixed") {
        subject_index <- as.integer(factor(sample_df$subject))
        n_subject <- max(subject_index)
    }

    center_biases <- numeric(reps)
    noise_variances <- numeric(reps)
    for (b in seq_len(reps)) {
        pseudo_beta <- matrix(NA_real_, nrow = n_fit, ncol = n_sample)
        for (j in seq_len(n_fit)) {
            fit <- fits[[j]]
            eta <- drop(X_full %*% fit$beta)
            if (mode == "mixed") {
                random_intercept <- stats::rnorm(
                    n_subject, mean = 0, sd = fit$tau
                )
                eta <- eta + random_intercept[subject_index]
            }

            mu <- .squeezeValues(stats::plogis(eta), squeeze_bounds)
            observed <- is.finite(beta_pilot[j, ])
            pseudo_beta[j, observed] <- stats::rbeta(
                sum(observed),
                shape1 = mu[observed] * common_phi,
                shape2 = (1 - mu[observed]) * common_phi
            )
        }

        pseudo_results <- BiocParallel::bplapply(
            seq_len(n_fit),
            function(j) {
                .fitOneCpG(
                    cpg_id = paste0("bootstrap_", b, "_", j),
                    y = pseudo_beta[j, ],
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

        pseudo_log_phi <- vapply(
            pseudo_results,
            function(x) {
                if (isTRUE(x$converged)) x$log_phi else NA_real_
            },
            numeric(1)
        )
        pseudo_prior <- .estimateLogPhiPrior(pseudo_log_phi)
        center_biases[b] <- pseudo_prior$mean - raw_prior$mean
        noise_variances[b] <- pseudo_prior$sd^2
    }

    center_bias <- mean(center_biases)
    noise_variance <- mean(noise_variances)
    corrected_variance <- max(raw_prior$sd^2 - noise_variance,
                              variance_floor)

    list(
        mean = raw_prior$mean - center_bias,
        sd = sqrt(corrected_variance)
    )
}

#' Squeeze beta values away from 0 and 1
#'
#' @param x Numeric vector of beta values.
#' @param bounds Length-2 numeric vector with lower and upper bounds.
#' @return Numeric vector with values clamped to `(bounds[1], bounds[2])`.
#' @noRd
.squeezeValues <- function(x, bounds = c(1e-5, 1 - 1e-5)) {
    pmin(pmax(x, bounds[1]), bounds[2])
}

#' Extract the assay and sample metadata used by RhyME
#' @noRd
.getRhyMEInput <- function(se, assay_name, time_col, subject_col = NULL,
                           covariates = NULL) {
    if (!inherits(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
    }

    cd <- as.data.frame(colData(se))
    missing_cols <- setdiff(c(time_col, subject_col, covariates), names(cd))
    if (length(missing_cols)) {
        stop(sprintf("Missing colData column(s): %s.",
                     paste(missing_cols, collapse = ", ")), call. = FALSE)
    }
    if (!is.numeric(cd[[time_col]])) {
        stop(sprintf("Column '%s' in colData must be numeric (hours).",
                     time_col), call. = FALSE)
    }
    if (anyNA(cd[[time_col]]) || any(!is.finite(cd[[time_col]]))) {
        stop(sprintf("Column '%s' in colData must contain finite times.",
                     time_col), call. = FALSE)
    }
    if (!is.null(subject_col)) {
        subject <- cd[[subject_col]]
        invalid_subject <- anyNA(subject) ||
            (is.numeric(subject) && any(!is.finite(subject))) ||
            any(!nzchar(as.character(subject)))
        if (invalid_subject) {
            stop(sprintf(
                "Column '%s' in colData contains missing or invalid subjects.",
                subject_col
            ), call. = FALSE)
        }
    }

    if (is.null(assay_name)) {
        beta <- assay(se, 1L)
        available <- assayNames(se)
        assay_name <- if (length(available) && !is.na(available[1]) &&
                          nzchar(available[1])) available[1] else NULL
    } else {
        beta <- assay(se, assay_name)
    }

    cpg_ids <- rownames(beta)
    if (is.null(cpg_ids) || length(cpg_ids) != nrow(beta) ||
        anyNA(cpg_ids) || any(!nzchar(cpg_ids)) || anyDuplicated(cpg_ids)) {
        stop("The fitted assay must have unique, non-missing CpG row names.",
             call. = FALSE)
    }
    if (!is.numeric(beta)) {
        stop("The fitted assay must contain numeric beta values.",
             call. = FALSE)
    }
    if (any(is.infinite(beta), na.rm = TRUE) ||
        any(beta < 0 | beta > 1, na.rm = TRUE)) {
        stop("Observed beta values must be finite and within [0, 1].",
             call. = FALSE)
    }

    list(beta = beta, col_data = cd, assay_name = assay_name)
}
