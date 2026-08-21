# Internal utility functions

.RHYME_MIN_PILOT_CPGS <- 1000L
.RHYME_PILOT_N <- 10000L
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

    # Samples without an observed collection time cannot contribute to the
    # cosinor design and are omitted from the fit.
    keep_samples <- is.finite(cd[[time_col]])
    beta <- beta[, keep_samples, drop = FALSE]
    cd <- cd[keep_samples, , drop = FALSE]

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
